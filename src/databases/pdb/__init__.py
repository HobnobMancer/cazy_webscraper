#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# Author:
# Emma E. M. Hobbs

# Contact
# ehobbbs@ebi.ac.uk

# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Get genome genome data from NCBI Assembly database"""

import gzip
import logging
import re
import urllib.error
import requests
import time
import sqlite3

#
# Bio.PDB reference:
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class
# implemented in Python. Bioinformatics 19: 2308-2310
"""Retrieve protein structure data and structure files from the RCSB PDB database."""


import logging
import sqlite3
import time

from argparse import Namespace
from pathlib import Path

import requests

from requests.exceptions import RequestException
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm

import Bio.PDB

from src.sql.interface.connect import get_sqlite3_connection
from src.sql.interface.temp_tables import create_temp_pdb_structure_table


logger = logging.getLogger(__name__)


# upper bound on the exponential retry backoff
MAX_BACKOFF_SECONDS = 60


RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"

# Bio.PDB's PDBList defaulted to ftp://ftp.wwpdb.org until biopython 1.85, and wwPDB has
# retired anonymous FTP, so on older biopython every download fails silently ("Desired
# structure doesn't exist") and writes no file. setup.py now requires biopython>=1.85, whose
# default is already https, but the server is still set explicitly here so the download does
# not depend on a library default that has changed once already.
PDB_SERVER = "https://files.wwpdb.org"

# One request returns every entry in the batch, so the number of round trips scales with
# the batch count rather than the accession count.
ENTRY_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    exptl { method }
    rcsb_entry_info { resolution_combined }
  }
}
"""


def fetch_pdb_metadata(batch: list[str], args: Namespace) -> tuple[dict[str, tuple], bool]:
    """Retrieve experimental method and resolution for a batch of PDB accessions.

    RCSB simply omits accessions it does not recognise from the response rather than
    erroring, so the caller identifies invalid ids by diffing the request against the result.

    Returns ({pdb_accession: (method, resolution)}, connection_error_flag).
    """
    metadata = {}
    max_attempts = max(1, args.retries)

    for attempt in range(1, max_attempts + 1):
        try:
            response = requests.post(
                RCSB_GRAPHQL_URL,
                json={"query": ENTRY_QUERY, "variables": {"ids": batch}},
                timeout=args.timeout,
            )
            response.raise_for_status()
            payload = response.json()

            if payload.get("errors"):
                logger.warning("RCSB returned errors for batch: %s", payload["errors"][:1])
                return metadata, True

            for entry in payload.get("data", {}).get("entries") or []:
                accession = entry.get("rcsb_id")
                if not accession:
                    continue

                methods = [m.get("method") for m in (entry.get("exptl") or []) if m.get("method")]
                method = ", ".join(methods) if methods else None

                # resolution_combined is a list, and is absent for methods that do not
                # produce one (e.g. solution NMR)
                resolutions = (entry.get("rcsb_entry_info") or {}).get("resolution_combined")
                resolution = resolutions[0] if resolutions else None

                metadata[accession] = (method, resolution)

            return metadata, False

        except (RequestException, ValueError) as err:
            logger.warning(
                "Attempt %s/%s failed retrieving PDB metadata: %s", attempt, max_attempts, err
            )
            if attempt < max_attempts:
                # capped: --retries defaults to 10, and an uncapped 2**attempt
                # would sleep ~17 min on the last try alone
                time.sleep(min(2 ** attempt, MAX_BACKOFF_SECONDS))

    return metadata, True


def dump_pdb_metadata(metadata: dict[str, tuple], conn: sqlite3.Connection) -> None:
    """Write a batch of retrieved PDB metadata straight into the temp table."""
    if not metadata:
        return

    cursor = conn.cursor()
    cursor.executemany(
        """INSERT OR REPLACE INTO TEMP_PDB_STRUCTURE (pdb_accession, method, resolution)
           VALUES (?, ?, ?)""",
        [(acc, method, resolution) for acc, (method, resolution) in metadata.items()],
    )
    conn.commit()
    cursor.close()


def download_structure_files(accessions: list[str], args: Namespace) -> int:
    """Download structure files for a batch of accessions from RCSB PDB.

    Returns the number of files that are present on disk after the call.
    """
    if args.skip_download or not args.file_formats:
        return 0

    pdbl = Bio.PDB.PDBList(server=PDB_SERVER)
    outdir = args.outdir if args.outdir else Path.cwd()
    downloaded = 0

    for file_format in args.file_formats:
        for accession in accessions:
            try:
                # retrieve_pdb_file is used rather than download_pdb_files so that a single
                # unavailable structure does not take the rest of the batch with it, and so
                # the count reflects files actually written rather than files requested
                written = pdbl.retrieve_pdb_file(
                    accession,
                    file_format=file_format,
                    overwrite=args.overwrite,
                    pdir=str(outdir),
                )
                if written and Path(written).exists():
                    downloaded += 1
                else:
                    logger.warning(
                        "No %s structure file available from PDB for %s",
                        file_format, accession,
                    )
            except Exception as err:
                logger.warning(
                    "Failed downloading %s structure file for %s: %s",
                    file_format, accession, err,
                )

    return downloaded


def get_pdb_data(
    pdb_accessions: list[str],
    cache_dir: Path,
    time_stamp: str,
    args: Namespace,
) -> dict[str, int]:
    """Retrieve PDB entry metadata and structure files, persisting metadata batch by batch.

    Metadata is written into TEMP_PDB_STRUCTURE as each batch comes back, so an interrupted
    run keeps everything it had already retrieved and can be resumed.

    Returns a dict of statistics about the operation.
    """
    conn = get_sqlite3_connection(args.database)
    create_temp_pdb_structure_table(conn)

    stats = {
        "metadata retrieved": 0,
        "invalid accessions": 0,
        "structure files downloaded": 0,
        "batches processed": 0,
        "failed batches": 0,
    }

    connection_err_cache = []
    invalid_id_cache = []

    batches = get_chunks_list(pdb_accessions, args.batch_size)

    for batch in tqdm(batches, desc="Retrieving PDB data"):
        metadata, connection_err = fetch_pdb_metadata(batch, args)

        if connection_err:
            stats["failed batches"] += 1
            connection_err_cache.extend(batch)
            continue

        dump_pdb_metadata(metadata, conn)
        stats["metadata retrieved"] += len(metadata)

        missing = [acc for acc in batch if acc not in metadata]
        invalid_id_cache.extend(missing)
        stats["invalid accessions"] += len(missing)

        # only ask PDB for files it told us exist
        retrievable = [acc for acc in batch if acc in metadata]
        if retrievable:
            stats["structure files downloaded"] += download_structure_files(retrievable, args)

        stats["batches processed"] += 1

    conn.close()

    if connection_err_cache:
        cache_path = cache_dir / f"pdb_connection_errors_{time_stamp}.txt"
        with open(cache_path, "a") as fh:
            for accession in connection_err_cache:
                fh.write(f"{accession}\n")

    if invalid_id_cache:
        cache_path = cache_dir / f"pdb_invalid_ids_{time_stamp}.txt"
        with open(cache_path, "a") as fh:
            for accession in invalid_id_cache:
                fh.write(f"{accession}\n")

    logger.warning("PDB retrieval completed:")
    logger.warning("  - Entry metadata retrieved: %d", stats["metadata retrieved"])
    logger.warning("  - Structure files downloaded: %d", stats["structure files downloaded"])
    logger.warning("  - Accessions not found in PDB: %d", stats["invalid accessions"])
    logger.warning("  - Batches processed: %d", stats["batches processed"])
    logger.warning("  - Failed batches: %d", stats["failed batches"])

    return stats
