#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2026
# (c) University of Strathclyde 2026
# (c) James Hutton Institute 2026
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# ehobbs@ebi.ac.uk
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Retrieve Pfam domain matches for UniProt accessions from the InterPro API.

InterPro's ``entry/pfam/protein/uniprot`` endpoint takes exactly one accession per request
(confirmed live: a comma-separated list of accessions returns 204 No Content rather than
being treated as a batch), so unlike the UniProt ID-mapping API there is no server-side batch
call to make here. "Batching" in this module therefore means what it means for the PDB and
NCBI Entrez retrieval modules: the accession list is chunked purely so progress is reported
and data is persisted batch by batch, while each accession within a batch is still fetched
with its own request.
"""


import logging
import time

from argparse import Namespace
from pathlib import Path

import requests

from requests.exceptions import RequestException
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm

from src.sql.interface.add_data.add_pfam_data import add_pfam_matches, merge_temp_pfam_relationships
from src.sql.interface.connect import get_sqlite3_connection
from src.sql.interface.temp_tables import create_temp_pfam_protein_table, drop_temp_pfam_protein_table

logger = logging.getLogger(__name__)


# upper bound on the exponential retry backoff, matching the PDB/GTDB/NCBI genome modules
MAX_BACKOFF_SECONDS = 60

# sanity cap on cursor-pagination per accession - comfortably above any protein's real Pfam
# match count, just there to stop a malformed 'next' link from looping forever
MAX_PAGES_PER_ACCESSION = 20

INTERPRO_MATCHES_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot"


class PfamMatchRecord:
    """One Pfam domain match location on a protein."""
    def __init__(
        self,
        uniprot_acc: str,
        pfam_accession: str,
        name: str,
        annotation_type: str,
        interpro_accession: str,
        match_start: int,
        match_end: int,
        release: str,
    ):
        self.uniprot_acc = uniprot_acc
        self.pfam_accession = pfam_accession
        self.name = name
        self.annotation_type = annotation_type
        self.interpro_accession = interpro_accession
        self.match_start = match_start
        self.match_end = match_end
        self.release = release
        self.protein_id = None  # filled in once mapped to a local db record


def fetch_pfam_matches(uniprot_acc: str, args: Namespace) -> tuple[list[dict], str, bool]:
    """Retrieve the raw Pfam match results for one UniProt accession, following pagination.

    The API returns 204 No Content both for a protein with no Pfam matches and for an
    accession it does not recognise - it does not distinguish the two - so that case is
    reported back as "no results" rather than as a connection error.

    Returns (raw result entries, InterPro release version, connection_error_flag).
    """
    url = f"{INTERPRO_MATCHES_URL}/{uniprot_acc}/"
    all_results = []
    release = None
    max_attempts = max(1, args.retries)
    pages = 0

    while url and pages < MAX_PAGES_PER_ACCESSION:
        payload = None

        for attempt in range(1, max_attempts + 1):
            try:
                response = requests.get(url, timeout=args.timeout)
                if response.status_code == 204:
                    return all_results, release, False
                response.raise_for_status()
                payload = response.json()
                release = response.headers.get("InterPro-Version", release)
                break
            except (RequestException, ValueError) as err:
                logger.warning(
                    "Attempt %s/%s failed retrieving InterPro Pfam matches for %s: %s",
                    attempt, max_attempts, uniprot_acc, err,
                )
                if attempt < max_attempts:
                    time.sleep(min(2 ** attempt, MAX_BACKOFF_SECONDS))

        if payload is None:
            return all_results, release, True

        all_results.extend(payload.get("results") or [])
        url = payload.get("next")
        pages += 1

    return all_results, release, False


def parse_pfam_results(uniprot_acc: str, results: list[dict], release: str) -> list[PfamMatchRecord]:
    """Flatten raw InterPro 'results' entries into one PfamMatchRecord per match location."""
    records = []

    for entry in results:
        metadata = entry.get("metadata") or {}
        if metadata.get("source_database") != "pfam":
            continue  # the endpoint queried is pfam-only, but guard against surprises anyway

        pfam_accession = metadata.get("accession")
        name = metadata.get("name")
        annotation_type = metadata.get("type")
        interpro_accession = metadata.get("integrated")

        for protein in entry.get("proteins") or []:
            for location in protein.get("entry_protein_locations") or []:
                for fragment in location.get("fragments") or []:
                    records.append(PfamMatchRecord(
                        uniprot_acc,
                        pfam_accession,
                        name,
                        annotation_type,
                        interpro_accession,
                        fragment.get("start"),
                        fragment.get("end"),
                        release,
                    ))

    return records


def get_pfam_data(
    uniprot_acc_to_protein_id: dict[str, int],
    cache_dir: Path,
    time_stamp: str,
    args: Namespace,
) -> dict[str, int]:
    """Get Pfam domain matches from InterPro and persist them batch by batch.

    :param uniprot_acc_to_protein_id: maps each UniProt accession to query to its local
        protein_id, used to attach retrieved matches to the right protein
    :param cache_dir: directory for caching accessions that failed or had no matches
    :param time_stamp: timestamp for the operation
    :param args: command line arguments

    Return a dict of statistics about the operation.
    """
    accessions = sorted(uniprot_acc_to_protein_id)
    batches = get_chunks_list(accessions, args.batch_size)

    stats = {
        "accessions queried": 0,
        "accessions with no matches": 0,
        "new pfam families": 0,
        "new matches": 0,
        "batches processed": 0,
        "failed accessions": 0,
    }

    conn = get_sqlite3_connection(args.database)
    create_temp_pfam_protein_table(conn)

    connection_err_cache = []
    no_match_cache = []

    for batch in tqdm(batches, desc="Retrieving Pfam data from InterPro"):
        batch_records = []

        for uniprot_acc in batch:
            results, release, connection_err = fetch_pfam_matches(uniprot_acc, args)
            stats["accessions queried"] += 1

            if connection_err:
                stats["failed accessions"] += 1
                connection_err_cache.append(uniprot_acc)
                continue

            if not results:
                stats["accessions with no matches"] += 1
                no_match_cache.append(uniprot_acc)
                continue

            records = parse_pfam_results(uniprot_acc, results, release)
            for record in records:
                record.protein_id = uniprot_acc_to_protein_id[uniprot_acc]
            batch_records.extend(records)

        if batch_records:
            new_family_count, new_match_count = add_pfam_matches(batch_records, conn)
            stats["new pfam families"] += new_family_count
            stats["new matches"] += new_match_count

        stats["batches processed"] += 1

    conn.close()

    merge_temp_pfam_relationships(args.database)

    cleanup_conn = get_sqlite3_connection(args.database)
    drop_temp_pfam_protein_table(cleanup_conn)
    cleanup_conn.close()

    if connection_err_cache:
        with open(cache_dir / f"pfam_connection_errors_{time_stamp}.txt", "w") as fh:
            for accession in connection_err_cache:
                fh.write(f"{accession}\n")

    if no_match_cache:
        with open(cache_dir / f"pfam_no_matches_{time_stamp}.txt", "w") as fh:
            for accession in no_match_cache:
                fh.write(f"{accession}\n")

    logger.info("Pfam retrieval from InterPro completed:")
    logger.info("  - Accessions queried: %d", stats["accessions queried"])
    logger.info("  - Accessions with no Pfam matches: %d", stats["accessions with no matches"])
    logger.info("  - New Pfam families: %d", stats["new pfam families"])
    logger.info("  - New protein-Pfam matches: %d", stats["new matches"])
    logger.info("  - Batches processed: %d", stats["batches processed"])
    logger.info("  - Accessions failed after retries: %d", stats["failed accessions"])

    return stats
