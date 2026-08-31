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

"""Retrieve taxonomic classifications from the Genome Taxonomy Database (GTDB)."""


import gzip
import logging
import sqlite3
import time

from argparse import Namespace
from pathlib import Path

import mechanicalsoup
import requests

from requests.exceptions import RequestException
from tqdm import tqdm

from src.sql.interface.connect import get_sqlite3_connection
from src.sql.interface.temp_tables import create_temp_gtdb_table


logger = logging.getLogger(__name__)


# upper bound on the exponential retry backoff
MAX_BACKOFF_SECONDS = 60

GTDB_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"

# GTDB publishes one taxonomy file per domain, named after the marker gene set used
DOMAIN_FILES = {
    "archaea": "ar53_taxonomy.tsv.gz",
    "bacteria": "bac120_taxonomy.tsv.gz",
}

# the seven ranks in a GTDB lineage string, in the order they appear, mapped to the
# GtdbTaxs columns they populate ("d__" is stored in the kingdom column)
LINEAGE_COLUMNS = (
    "kingdom", "phylum", "tax_class", "tax_order", "family", "genus", "species",
)


def get_gtdb_release() -> str:
    """Read the GTDB release version (e.g. 'v232') from the release's VERSION.txt.

    The v2 code derived this from the data file name, which yields the marker gene set
    ('ar53'/'bac120') rather than the release, so every lineage was stored with the wrong
    release string. VERSION.txt is the authoritative source.

    Returns the release string, or None if it could not be read.
    """
    try:
        response = requests.get(f"{GTDB_URL}VERSION.txt", timeout=30)
        response.raise_for_status()
        # the file starts with the version on its own line, e.g. "v232"
        release = response.text.strip().split("\n")[0].strip()
        return release or None
    except RequestException as err:
        logger.warning("Could not retrieve the GTDB release version: %s", err)
        return None


def download_gtdb_file(file_name: str, out_path: Path, args: Namespace) -> Path | None:
    """Download one GTDB taxonomy file, streaming it to disk.

    Returns the path written, or None if the download failed.
    """
    if out_path.exists():
        logger.info("Using already downloaded GTDB file %s", out_path)
        return out_path

    url = f"{GTDB_URL}{file_name}"
    max_attempts = max(1, args.retries)

    for attempt in range(1, max_attempts + 1):
        try:
            with requests.get(url, stream=True, timeout=args.timeout) as response:
                response.raise_for_status()
                total = int(response.headers.get("Content-Length", 0)) or None

                with open(out_path, "wb") as fh, tqdm(
                    total=total, unit="B", unit_scale=True,
                    desc=f"Downloading {file_name}",
                ) as pbar:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            fh.write(chunk)
                            pbar.update(len(chunk))

            return out_path

        except RequestException as err:
            logger.warning(
                "Attempt %s/%s failed downloading %s: %s", attempt, max_attempts, file_name, err
            )
            out_path.unlink(missing_ok=True)  # do not leave a truncated file behind
            if attempt < max_attempts:
                time.sleep(min(2 ** attempt, MAX_BACKOFF_SECONDS))

    logger.error("Could not download GTDB file %s", file_name)
    return None


def parse_lineages_to_db(
    gtdb_file: Path,
    accession_to_genome_id: dict,
    release: str,
    conn: sqlite3.Connection,
    batch_size: int = 1000,
) -> int:
    """Stream a GTDB taxonomy file, dumping the lineages of genomes of interest into the db.

    The bacterial file is ~120 MB and several hundred thousand rows, so it is read line by
    line and only matching rows are kept; nothing is loaded into a dataframe. Rows are written
    to TEMP_GTDB in batches as they are found, so an interrupted run keeps its progress.

    Returns the number of lineages written.
    """
    written = 0
    batch = []
    cursor = conn.cursor()

    def flush():
        nonlocal written, batch
        if not batch:
            return
        cursor.executemany(
            """INSERT OR REPLACE INTO TEMP_GTDB
               (assembly_accession, genome_id, kingdom, phylum, tax_class, tax_order,
                family, genus, species, release)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            batch,
        )
        conn.commit()
        written += len(batch)
        batch = []

    with gzip.open(gtdb_file, "rt", encoding="utf-8") as fh:
        for line in tqdm(fh, desc=f"Parsing {gtdb_file.name}"):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            # accessions are prefixed by their source archive: RS_ (RefSeq) or GB_ (GenBank)
            accession = parts[0].replace("RS_", "", 1).replace("GB_", "", 1)
            genome_id = accession_to_genome_id.get(accession)
            if genome_id is None:
                continue

            # "d__Archaea;p__Methanobacteriota;...;s__Methanocatella smithii"
            ranks = {}
            for rank in parts[1].split(";"):
                prefix, _, name = rank.partition("__")
                ranks[prefix.strip()] = name.strip() or None

            lineage = [
                ranks.get(prefix)
                for prefix in ("d", "p", "c", "o", "f", "g", "s")
            ]

            batch.append((accession, genome_id, *lineage, release))
            if len(batch) >= batch_size:
                flush()

    flush()
    cursor.close()

    return written


def get_gtdb_taxonomies(
    accession_to_genome_id: dict,
    cache_dir: Path,
    args: Namespace,
) -> dict[str, int]:
    """Download the GTDB release files and dump the lineages of interest into the local db.

    Returns a dict of statistics about the operation.
    """
    conn = get_sqlite3_connection(args.database)
    create_temp_gtdb_table(conn)

    stats = {"lineages retrieved": 0, "files parsed": 0, "files failed": 0}

    release = get_gtdb_release()
    if release is None:
        logger.warning(
            "Could not determine the GTDB release version; lineages will be stored without one"
        )

    for domain in args.taxs:
        file_name = DOMAIN_FILES[domain]
        out_path = cache_dir / file_name

        gtdb_file = download_gtdb_file(file_name, out_path, args)
        if gtdb_file is None:
            stats["files failed"] += 1
            continue

        stats["lineages retrieved"] += parse_lineages_to_db(
            gtdb_file, accession_to_genome_id, release, conn
        )
        stats["files parsed"] += 1

        if not args.nodelete_cache:
            gtdb_file.unlink(missing_ok=True)

    conn.close()

    logger.warning("GTDB retrieval completed:")
    logger.warning("  - GTDB release: %s", release)
    logger.warning("  - Lineages retrieved for genomes of interest: %d", stats["lineages retrieved"])
    logger.warning("  - Release files parsed: %d", stats["files parsed"])
    if stats["files failed"]:
        logger.warning("  - Release files that could not be downloaded: %d", stats["files failed"])

    return stats
