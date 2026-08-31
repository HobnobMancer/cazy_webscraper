#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2025
# (c) University of Strathclyde 2025
# (c) James Hutton Institute 2025
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
"""Retrieve and parse data from UniProt KB"""


import sqlite3
import logging
import time


from argparse import Namespace
from pathlib import Path

from bioservices import UniProt as UniProtService
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm

from src.sql.interface.add_data.add_uniprot_data import (
    add_uniprot_records, add_ec_numbers, add_pdbs, add_go_terms,
    merge_temp_ec_relationships, merge_temp_pdb_relationships, merge_temp_go_relationships
)
from src.sql.interface.get_data.get_proteins import get_ncbi_acc_to_id
from src.sql.interface.connect import get_sqlite3_connection
from src.sql.interface.temp_tables import (
    create_temp_ec_protein_table,
    create_temp_pdb_protein_table,
    create_temp_go_protein_table,
    drop_temp_ec_protein_table,
    drop_temp_pdb_protein_table,
    drop_temp_go_protein_table
)

logger = logging.getLogger(__name__)


class UniProtRecord:
    def __init__(
        self,
        ncbi_acc: str, uniprot_id: str, uniparc_id: str,
        name: str, swissprot: bool, sequence: str, md5: str,
        mol_weight: float, sequence_date: str,
        ec_nums: set[str],
        pdbs: set[tuple[str]],
        go_terms: dict[str, str]
    ):
        self.ncbi_acc = ncbi_acc
        self.uniprot_id = uniprot_id
        self.uniparc_id = uniparc_id
        self.name = name
        self.swissprot = swissprot
        self.sequence = sequence
        self.md5 = md5
        self.mol_weight = mol_weight
        self.sequence_date = sequence_date
        self.ec_nums = ec_nums
        self.pdbs = pdbs
        self.go_terms = go_terms
        self.protein_id = None  # to be filled when mapped to local db record


class UniProtKB:
    """
    UniProtKB service class for retrieving and parsing UniProt records.

    Bioservices requests batch queries no larger than 200.

    Note that according to Uniprot (June 2022), there are various limits on ID Mapping Job Submission:

    ========= =====================================================================================
    Limit	  Details
    ========= =====================================================================================
    100,000	  Total number of ids allowed in comma separated param ids in /idmapping/run api
    500,000	  Total number of "mapped to" ids allowed
    100,000	  Total number of "mapped to" ids allowed to be enriched by UniProt data
    10,000	  Total number of "mapped to" ids allowed with filtering
    ========= =====================================================================================
    """
    def __init__(self, retries: int = 3):
        self.retries = retries
        self.service = UniProtService()

    def _map_batch(self, batch: list[str], fr: str) -> tuple[list[dict], list[str]]:
        """Submit one UniProt ID mapping job for a batch of ids using the given 'from' id type.

        bioservices swallows all connection-level errors itself and surfaces a connection
        failure to us in one of two ways depending on which step failed: mapping() returns a
        falsy value if the job submission itself failed, or raises TypeError if the job was
        submitted but checking on its status failed (confirmed against the live API/bioservices
        source - in neither case does a raw requests exception ever reach this code). Both are
        connection failures and warrant a retry, not just the TypeError case.

        Returns a tuple of (raw "results" entries, failed ids).
        """
        mappings = None
        attempt = 0
        while not mappings and attempt <= self.retries:
            try:
                mappings = self.service.mapping(
                    fr=fr,
                    to="UniProtKB",
                    query=",".join(batch),  # str of ids, separated by commas
                    progress=False
                )
            except TypeError:
                # raised by bioservices if the job was submitted but its status check failed
                mappings = None

            if not mappings:
                if attempt == self.retries:
                    break
                logger.warning("Could not connect to UniProt KB, retrying (%d/%d)", attempt + 1, self.retries)
                time.sleep(10)
                attempt += 1

        if mappings:
            return mappings.get("results", []), mappings.get("failedIds", [])
        return [], batch  # return all ids as failed if no mappings

    def get_uniprot_records(self, batch: list[str], swissprot_only: bool, get_ecs: bool, get_pdbs: bool, get_gos: bool) -> tuple[list, list[str]]:
        """Map NCBI protein accessions to UniProt IDs.

        Accessions are first mapped as classic GenBank/EMBL/DDBJ CDS ids. UniProt indexes
        RefSeq accessions (e.g. NP_/XP_/WP_ prefixed) under a separate id type, so anything
        that fails the first pass is retried as a RefSeq id before being treated as a genuine
        miss - confirmed against the live API that EMBL-GenBank-DDBJ_CDS alone never matches
        RefSeq accessions, even when the record exists in UniProt under RefSeq_Protein.

        Returns a tuple with (parsed_records, failed_ids)
        """
        results, failed_ids = self._map_batch(batch, "EMBL-GenBank-DDBJ_CDS")

        if failed_ids:
            refseq_results, failed_ids = self._map_batch(failed_ids, "RefSeq_Protein")
            results += refseq_results

        parsed_records = self.parse_mappings(results, swissprot_only, get_ecs, get_pdbs, get_gos)
        return parsed_records, failed_ids

    def parse_mappings(
        self,
        mappings: list[dict],
        swissprot_only: bool = False,
        get_ecs: bool = True,
        get_pdbs: bool = True,
        get_gos: bool = True
    ) -> list[UniProtRecord]:
        """Parse UniProt mapping results into list of UniProtKB instances"""
        parsed_mappings = []
        for entry in mappings:
            uniprot_id = uniparc_id = name = swissprot = sequence = md5 = mol_weight = sequence_date = None
            ec_nums = set()
            pdbs = set()
            go_terms = {}

            ncbi_acc = entry.get("from")
            record = entry.get("to", {})

            uniprot_id = record.get("uniProtkbId")
            uniparc_id = record.get("extraAttributes", {}).get("uniParcId")

            protein_desc = record.get("proteinDescription", {})
            if protein_desc:
                name = protein_desc.get("recommendedName", {}).get("fullName", {}).get("value")
                if not name:
                    for submitted_name in protein_desc.get("submittedNames", []):
                        name = submitted_name.get("fullName", {}).get("value")
                        if name:
                            break
                if get_ecs:
                    ec_nums = set([value['value'] for value in protein_desc.get("recommendedName", {}).get("ecNumbers", [])])

            if get_ecs:
                for comment in record.get("comments", []):
                    if comment.get("commentType") == "CATALYTIC ACTIVITY":
                        ec_nums.add(comment.get("reaction", {}).get("ecNumber"))

            swissprot = True if (record.get("entryType")).lower() == "uniprotkb reviewed (swiss-prot)" else False
            if swissprot_only and not swissprot:
                continue  # skip non-swissprot entries

            sequence = record.get("sequence", {}).get("value")
            md5 = record.get("sequence", {}).get("md5")
            mol_weight = record.get("sequence", {}).get("molWeight")
            sequence_date = record.get("entryAudit", {}).get("lastSequenceUpdateDate")

            if get_pdbs:
                for value in record.get("uniProtKBCrossReferences", []):
                    if value.get("database") == "PDB":
                        method = resolution = None
                        for prop in value.get("properties", []):
                            if prop.get("key") == "Method":
                                method = prop.get("value")
                            if prop.get("key") == "Resolution":
                                try:
                                    resolution = float(prop.get("value").rstrip(" A").strip())
                                except ValueError:
                                    resolution = None
                        pdbs.add((value.get("id"), method, resolution))

            if get_gos:
                for value in record.get("uniProtKBCrossReferences", []):
                    if value.get("database") == "GO":
                        go_id = value.get("id")
                        for prop in value.get("properties", []):
                            if prop.get("key") == "GoTerm":
                                go_term = prop.get("value")
                                go_terms[go_id] = go_term

            parsed_mappings.append(UniProtRecord(
                ncbi_acc,
                uniprot_id,
                uniparc_id,
                name,
                swissprot,
                sequence,
                md5,
                mol_weight,
                sequence_date,
                ec_nums,
                pdbs,
                go_terms
            ))

        return parsed_mappings


def get_uniprot_data(
    protein_accs: list[str],
    cache_dir: Path,
    time_stamp: str,
    args: Namespace
) -> dict[str, int]:
    """Get protein data from UniProt KB and persist them batch by batch.

    Args:
        protein_accs: List of protein accessions to retrieve
        cache_dir: Directory for caching failed requests
        time_stamp: Timestamp for the operation
        args: Command line arguments

    Returns:
        Dictionary with statistics about the operation
    """
    batches = get_chunks_list(protein_accs, args.batch_size)
    stats = {
        "proteins not in uniprot": 0,
        "uniprot ids retrieved": 0,
        "new uniprot records": 0,
        "sequences updated": 0,
        "new ecs": 0,
        "new pdbs": 0,
        "new go terms": 0,
        "batches processed": 0,
        "failed batches": 0
    }
    # create these up front rather than relying on the add_* helpers to create them
    # lazily: if no batch yields any EC/PDB/GO the merge step below would otherwise hit a
    # table that was never created
    setup_conn = get_sqlite3_connection(args.database)
    if args.ec:
        create_temp_ec_protein_table(setup_conn)
    if args.pdb:
        create_temp_pdb_protein_table(setup_conn)
    if args.go:
        create_temp_go_protein_table(setup_conn)
    setup_conn.close()

    connection_err_cache = open(cache_dir / f"uniprot_connection_errors_{time_stamp}.txt", "w")
    failed_ids_cache = open(cache_dir / f"uniprot_failed_ids_{time_stamp}.txt", "w")
    uniprotkb_service = UniProtKB(args.retries)

    def process_batch(batch: list[str]) -> None:
        """Process a single batch of accessions."""
        nonlocal stats
        nonlocal connection_err_cache
        nonlocal failed_ids_cache
        nonlocal uniprotkb_service

        records, failed_ids = uniprotkb_service.get_uniprot_records(batch, args.swissprot, args.ec, args.pdb, args.go)

        for failed_id in failed_ids:
            failed_ids_cache.write(f"{failed_id}\n")

        if not records:  # could not connect to uniprot
            stats["failed batches"] += 1
            connection_err_cache.write(f"Failed to retrieve batch: {','.join(batch)}\n")
            return

        stats["uniprot ids retrieved"] += len(records)

        new_record_count, updated_seq_count = add_uniprot_records(records, args.database, args.update)
        stats["new uniprot records"] += new_record_count
        stats["sequences updated"] += updated_seq_count

        if args.ec or args.pdb or args.go:
            # map protein ids to local db records
            conn = sqlite3.connect(args.database)
            ncbi2id = get_ncbi_acc_to_id(conn, set([record.ncbi_acc for record in records]))
            for record in records:
                if record.ncbi_acc in ncbi2id:
                    record.protein_id = ncbi2id[record.ncbi_acc]

            if args.ec:
                stats["new ecs"] += add_ec_numbers(records, conn)
            if args.pdb:
                stats["new pdbs"] += add_pdbs(records, conn)
            if args.go:
                stats["new go terms"] += add_go_terms(records, conn)

            conn.commit()
            conn.close()

    for batch in tqdm(batches, desc="Retrieving UniProt data"):
        process_batch(batch)
        stats["batches processed"] += 1

    # add all new ec, pdb and go relationships
    if args.ec:
        merge_temp_ec_relationships(args.database)
    if args.pdb:
        merge_temp_pdb_relationships(args.database)
    if args.go:
        merge_temp_go_relationships(args.database)

    # the drop_temp_* helpers take an open connection, not a db path
    cleanup_conn = get_sqlite3_connection(args.database)
    if args.ec:
        drop_temp_ec_protein_table(cleanup_conn)
    if args.pdb:
        drop_temp_pdb_protein_table(cleanup_conn)
    if args.go:
        drop_temp_go_protein_table(cleanup_conn)
    cleanup_conn.close()

    connection_err_cache.close()
    failed_ids_cache.close()

    logger.info("UniProt data retrieval completed:")
    logger.info("  - Proteins not in UniProt: %d", stats['proteins not in uniprot'])
    logger.info("  - UniProt ids retrieved: %d", stats['uniprot ids retrieved'])
    logger.info("  - New UniProt records: %d", stats['new uniprot records'])
    logger.info("  - Sequences updated: %d", stats['sequences updated'])
    logger.info("  - New EC numbers: %d", stats['new ecs'])
    logger.info("  - New PDB accessions: %d", stats['new pdbs'])
    logger.info("  - New GO terms: %d", stats['new go terms'])
    logger.info("  - Batches processed: %d", stats['batches processed'])
    logger.info("  - Failed batches: %d", stats['failed batches'])

    return stats
