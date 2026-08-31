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

from argparse import Namespace
from http.client import IncompleteRead
from pathlib import Path
from requests.exceptions import RequestException

from Bio import Entrez
from Bio.Entrez.Parser import NotXMLError, CorruptedXMLError
from saintBioutils.genbank import entrez_retry
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm

from src.cache.ncbi import cache_connection_errors, cache_invalid_ids
from src.databases.ncbi import post_acc_to_entrez, validate_batch, handle_entrez_errors
from src.sql.interface.connect import get_sqlite3_connection
from src.sql.interface.temp_tables import create_temp_genome_tables
from src.sql.interface.get_data.get_genomes import get_cached_proteins


logger = logging.getLogger(__name__)


# upper bound on the exponential retry backoff
MAX_BACKOFF_SECONDS = 60


# RefSeq protein accessions are two letters + underscore (NP_, XP_, WP_, YP_, ...).
# GenBank/EMBL/DDBJ protein accessions never carry an underscore in that position.
REFSEQ_PROTEIN_RE = re.compile(r"^[A-Z]{2}_\d")


def is_refseq_protein(accession: str) -> bool:
    """Return True if the accession is a RefSeq protein accession."""
    return bool(REFSEQ_PROTEIN_RE.match(accession))


class NcbiGenome:
    """Represent one genome assembly record from the NCBI Assembly database"""

    def __init__(self):
        self.genome_id = None         # NCBI Assembly UID
        self.assembly_name = None
        self.gbk_accession = None     # GCA_...
        self.refseq_accession = None  # GCF_...

    def parse_genome_record(self, record) -> None:
        """Parse one DocumentSummary returned by esummary(db="assembly").

        The summary carries both the GenBank (GCA_) and the RefSeq (GCF_) accession of the
        same assembly under 'Synonym'. Both are kept because they are not interchangeable
        downstream: a GenBank protein accession is only listed in the GCA feature table and
        a RefSeq protein accession only in the GCF one.
        """
        try:
            self.genome_id = record.attributes.get('uid')
            self.assembly_name = (record.get('AssemblyName') or '').strip() or None

            synonyms = record.get('Synonym', {})
            self.gbk_accession = (synonyms.get('Genbank') or '').strip() or None
            self.refseq_accession = (synonyms.get('RefSeq') or '').strip() or None

            # older/partial summaries may omit Synonym - fall back to the primary accession
            if not self.gbk_accession and not self.refseq_accession:
                accession = (record.get('AssemblyAccession') or '').strip() or None
                if accession and accession.startswith('GCF_'):
                    self.refseq_accession = accession
                elif accession:
                    self.gbk_accession = accession
        except (KeyError, AttributeError) as err:
            logger.warning("Error parsing assembly record: %s", err)

    def compile_url(self, accession: str) -> str:
        """Compile the NCBI FTP URL of the feature table for one assembly accession."""
        prefix, digits = accession.split("_")
        return (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
            f"{prefix}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/"
            f"{accession}_{self.assembly_name}/"
            f"{accession}_{self.assembly_name}_feature_table.txt.gz"
        )

    def save_to_temp_table(self, db_connection: sqlite3.Connection) -> None:
        """Save assembly data to the temp table"""
        cursor = db_connection.cursor()
        cursor.execute(
            """INSERT OR REPLACE INTO TEMP_GENOME
               (ncbi_genome_id, gbk_accession, refseq_accession, assembly_name)
               VALUES (?, ?, ?, ?)""",
            (self.genome_id, self.gbk_accession, self.refseq_accession, self.assembly_name)
        )
        db_connection.commit()

    def parse_feature_table_to_db(
        self,
        feature_table_path: Path,
        protein_accs: set[str],
        db_connection: sqlite3.Connection
    ) -> int:
        """Parse a feature table and store the protein->assembly mappings it evidences.

        Args:
            feature_table_path: Path to the (gzipped) feature table
            protein_accs: Protein accessions to look for in this table
            db_connection: Database connection

        Returns the number of mappings written.
        """
        mappings_written = 0
        protein_accs = list(protein_accs)
        if not protein_accs:
            return 0

        try:
            cursor = db_connection.cursor()
            placeholders = ','.join('?' * len(protein_accs))
            cursor.execute(
                f"SELECT protein_accession, protein_id FROM Proteins "
                f"WHERE protein_accession IN ({placeholders})",
                protein_accs
            )
            acc_to_id = dict(cursor.fetchall())
            wanted = set(protein_accs)

            protein_mappings = []
            batch_size = 500

            def flush():
                nonlocal mappings_written
                if not protein_mappings:
                    return
                cursor.executemany(
                    """INSERT OR REPLACE INTO TEMP_GENOME2PROTEIN
                       (ncbi_genome_id, protein_id)
                       VALUES (?, ?)""",
                    protein_mappings
                )
                db_connection.commit()
                mappings_written += len(protein_mappings)
                protein_mappings.clear()

            with gzip.open(feature_table_path, "rt", encoding='utf-8') as f:
                for line in f:
                    if not line.startswith("CDS"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 12:
                        continue
                    # col 10 is product_accession (GenBank ids in a GCA table, RefSeq ids
                    # in a GCF table); col 11 is non-redundant_refseq, which some GCA
                    # tables populate with the equivalent WP_ accession
                    for prot_acc in (parts[10].strip(), parts[11].strip()):
                        if prot_acc and prot_acc in wanted:
                            protein_mappings.append((self.genome_id, acc_to_id[prot_acc]))
                            if len(protein_mappings) >= batch_size:
                                flush()
                            break

            flush()

            logger.info(
                "Mapped %d proteins to assembly %s from its feature table",
                mappings_written, self.gbk_accession or self.refseq_accession
            )

            try:
                if feature_table_path.exists():
                    feature_table_path.unlink()
                    logger.debug("Deleted feature table %s", feature_table_path)
            except OSError as del_err:
                logger.warning("Failed to delete feature table %s: %s", feature_table_path, del_err)

        except Exception as e:
            logger.error("Error parsing feature table %s: %s", feature_table_path, e)

        return mappings_written


def get_ncbi_genome_data(
    protein_accs: list[str],
    cache_dir: Path,
    args: Namespace
) -> None:
    """Get genome data for proteins from NCBI Assembly database and store in temp tables

    Args:
        protein_accs: List of protein accessions
        cache_dir: Directory for caching failed requests
        db_connection: Database connection
        args: Command line arguments
    """
    conn = get_sqlite3_connection(args.database)

    create_temp_genome_tables(conn)

    # Check which proteins are already cached
    cached_proteins = get_cached_proteins(protein_accs, conn)
    uncached_proteins = [acc for acc in protein_accs if acc not in cached_proteins]

    if cached_proteins:
        logger.info("Found %d proteins already cached, processing %d new proteins",
                   len(cached_proteins), len(uncached_proteins))

    if not uncached_proteins:
        logger.info("All proteins already cached, skipping NCBI retrieval")
        return

    batches = get_chunks_list(uncached_proteins, args.batch_size)
    all_genome_ids = set()
    connection_err_cache = []
    invalid_id_cache = []

    def process_batch(batch: list[str]) -> None:
        nonlocal connection_err_cache, invalid_id_cache

        validated_batch, connection_err_cache = validate_batch(batch, connection_err_cache, args)
        if not validated_batch:
            return

        epost_webenv, epost_query_key, invalid_err, connection_err = post_acc_to_entrez(validated_batch, args)
        if connection_err:
            connection_err_cache.append(validated_batch)
            return
        if invalid_err:
            # If batch size is 1, cache as invalid ID
            if len(validated_batch) == 1:
                invalid_id_cache.append(validated_batch[0])
            else:
                # Split batch in half and retry
                mid = len(validated_batch) // 2
                process_batch(validated_batch[:mid])
                process_batch(validated_batch[mid:])
            return

        # Link proteins to assemblies
        genome_links, connection_err = link_proteins_to_assemblies(
            epost_webenv, epost_query_key, args
        )
        if connection_err:
            connection_err_cache.append(validated_batch)
            return

        all_genome_ids.update(genome_links)

    for batch in tqdm(batches, desc="Fetching genome links from NCBI"):
        process_batch(batch)

    # Fetch and process genome data
    if all_genome_ids:
        # split the targets by accession type once: a GenBank protein accession is only
        # listed in the GCA feature table and a RefSeq one only in the GCF table, so each
        # assembly only needs the table(s) matching the accessions we are trying to map
        refseq_targets = {acc for acc in uncached_proteins if is_refseq_protein(acc)}
        gbk_targets = set(uncached_proteins) - refseq_targets

        genome_batches = get_chunks_list(list(all_genome_ids), args.batch_size)
        for batch in tqdm(genome_batches, desc="Fetching genome data from NCBI"):
            batch_genomes, connection_err = fetch_and_store_genome_records(batch, conn, args)
            if connection_err:
                connection_err_cache.append(batch)
                continue

            # Download feature tables and parse proteins for each genome
            for genome in batch_genomes:
                for accession, targets in (
                    (genome.gbk_accession, gbk_targets),
                    (genome.refseq_accession, refseq_targets),
                ):
                    if not accession or not targets:
                        continue
                    feature_table_path = download_feature_table(genome, accession, cache_dir, args)
                    if feature_table_path:
                        genome.parse_feature_table_to_db(feature_table_path, targets, conn)

    # Cache any errors
    if connection_err_cache:
        cache_connection_errors(connection_err_cache, cache_dir)

    if invalid_id_cache:
        cache_invalid_ids(invalid_id_cache, cache_dir)

    conn.commit()
    conn.close()

    return len(all_genome_ids)


def fetch_and_store_genome_records(
    genome_ids: list[str],
    db_connection: sqlite3.Connection,
    args: Namespace
) -> tuple[list[NcbiGenome], bool]:
    """Fetch genome records from NCBI and store in temp table

    Args:
        genome_ids: List of genome IDs to fetch
        db_connection: Database connection
        args: Command line arguments

    Returns:
        Tuple of (list of NcbiGenome objects, error occurred flag)
    """
    genomes = []
    connection_err = False

    if not genome_ids:
        return genomes, connection_err

    try:
        with entrez_retry(
            args.retries,
            Entrez.esummary,
            db="assembly",
            id=",".join(genome_ids),
            retmode="xml"
        ) as handle:
            records = Entrez.read(handle, validate=False)

            # esummary against the assembly db nests the records one level down
            for record in records['DocumentSummarySet']['DocumentSummary']:
                genome = NcbiGenome()
                genome.parse_genome_record(record)
                if genome.assembly_name and (genome.gbk_accession or genome.refseq_accession):
                    genome.save_to_temp_table(db_connection)
                    genomes.append(genome)

    except urllib.error.HTTPError as err:
        if err.code == 500:
            logger.warning("NCBI server error (500). Server is experiencing issues.")
        else:
            logger.warning("HTTP error %s: %s", err.code, err)
        connection_err = True
    except (IncompleteRead, CorruptedXMLError, NotXMLError) as err:
        logger.warning("Entrez connection error: %s", err)
        connection_err = True
    except Exception as err:
        handle_entrez_errors(err)
        connection_err = True

    return genomes, connection_err


def link_proteins_to_assemblies(
    epost_webenv: str,
    epost_query_key: str,
    args: Namespace
) -> tuple[list[str], bool]:
    """Link protein records to the genome assemblies they belong to.

    NCBI offers no direct protein->assembly link: `acheck` against the live API confirms
    that a protein UID exposes `protein_nuccore` (always) and `protein_genome` (RefSeq
    proteins only). The previous implementation used `protein_genome`, which meant GenBank
    protein accessions - the bulk of CAZy - resolved to nothing at all, and RefSeq ones
    resolved to a species-level Genome entry that points at the species' representative
    assembly rather than the assembly the protein actually came from.

    So the link is made in two hops instead: protein -> nuccore (the specific nucleotide
    record encoding the protein) -> assembly.

    Returns:
        Tuple of (list of NCBI Assembly UIDs, error occurred flag)
    """
    assembly_ids = []
    connection_err = False

    try:
        # hop 1: protein -> nuccore
        with entrez_retry(
            args.retries,
            Entrez.elink,
            query_key=epost_query_key,
            WebEnv=epost_webenv,
            dbfrom="protein",
            db="nuccore",
            linkname="protein_nuccore"
        ) as handle:
            link_results = Entrez.read(handle, validate=False)

        nuccore_ids = [
            link['Id']
            for result in link_results
            for hit in result.get('LinkSetDb', [])
            for link in hit.get('Link', [])
        ]

        if not nuccore_ids:
            return assembly_ids, connection_err

        # hop 2: nuccore -> assembly
        with entrez_retry(
            args.retries,
            Entrez.elink,
            id=nuccore_ids,
            dbfrom="nuccore",
            db="assembly",
            linkname="nuccore_assembly"
        ) as handle:
            link_results = Entrez.read(handle, validate=False)

        assembly_ids = [
            link['Id']
            for result in link_results
            for hit in result.get('LinkSetDb', [])
            for link in hit.get('Link', [])
        ]

    except urllib.error.HTTPError as err:
        if err.code == 500:
            logger.warning("NCBI server error (500). Server is experiencing issues.")
        else:
            logger.warning("HTTP error %s: %s", err.code, err)
        connection_err = True
    except (IncompleteRead, CorruptedXMLError, NotXMLError) as err:
        logger.warning("Entrez connection error: %s", err)
        connection_err = True
    except Exception as err:
        handle_entrez_errors(err)
        connection_err = True

    return assembly_ids, connection_err


def download_feature_table(
    genome: NcbiGenome,
    accession: str,
    download_dir: Path,
    args: Namespace
) -> Path | None:
    """Download the feature table for one assembly accession.

    Args:
        genome: NcbiGenome object
        accession: the assembly accession (GCA_/GCF_) whose feature table to fetch
        download_dir: Directory to save downloaded file
        args: Command line arguments

    Returns:
        Path to the downloaded feature table file or None if download failed
    """
    url = genome.compile_url(accession)
    # First check if the URL is valid
    try:
        head_response = requests.head(url, timeout=10)
        if head_response.status_code == 404:
            logger.warning(
                "Feature table not available for %s at %s. "
                "Cannot map genome back to proteins.",
                accession, url
            )
            return None
        elif head_response.status_code != 200:
            logger.warning(
                "Feature table URL returned status %s for "
                "%s. Cannot map genome back to proteins.",
                head_response.status_code, accession
            )
            return None
    except RequestException as e:
        logger.warning(
            "Cannot access feature table for %s at %s: %s. "
            "Cannot map genome back to proteins.",
            accession, url, e
        )
        return None

    feature_table_path = download_dir / f"{accession}.feature_table.txt.gz"

    # Check if already downloaded
    if feature_table_path.exists():
        return feature_table_path

    max_attempts = max(1, args.retries)

    for attempt in range(1, max_attempts + 1):
        try:
            response = requests.get(url, stream=True, timeout=30)
            response.raise_for_status()

            with open(feature_table_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

            return feature_table_path

        except RequestException as e:
            logger.warning(
                "Attempt %s/%s failed downloading feature table for "
                "%s: %s",
                attempt, max_attempts, accession, e
            )
            if attempt < max_attempts:
                # exponential backoff before next try
                # capped: --retries defaults to 10, and an uncapped 2**attempt
                # would sleep ~17 min on the last try alone
                time.sleep(min(2 ** attempt, MAX_BACKOFF_SECONDS))
            else:
                logger.error(
                    "Error downloading feature table for %s after "
                    "%s attempts: %s. Cannot map genome back to proteins.",
                    accession, max_attempts, e
                )
                return None
        except Exception as e:
            logger.error("Unexpected error downloading feature table for %s: %s", accession, e)
            return None

    return None
