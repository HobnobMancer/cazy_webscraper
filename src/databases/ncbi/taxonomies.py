#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) James Hutton Institute 2024
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
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

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Get taxonomies from the NCBI Taxonomy database"""


import sqlite3
import logging

from argparse import Namespace
from http.client import IncompleteRead
from pathlib import Path

from Bio import Entrez
from Bio.Entrez.Parser import NotXMLError, CorruptedXMLError
from saintBioutils.genbank import entrez_retry
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm

from src.cache.ncbi import cache_connection_errors, cache_invalid_ids
from src.databases.ncbi import handle_entrez_errors, post_acc_to_entrez, validate_batch
from src.sql.interface.connect import get_sqlite3_connection
from src.sql.interface.temp_tables import create_temp_tax_tables, drop_temp_tax_tables


logger = logging.getLogger(__name__)


class NcbiProtein:
    """Represent the data extract from a NCBI protein record.
    Used when scrapping CAZy and dealing with multiple Taxa."""
    def __init__(self):
        self.protein_id = None
        self.genus = None
        self.species = None
        self.kingdom = None

    def get_ncbi_data(self, record):
        """extract the data from the record"""
        self.protein_id = record['GBSeq_accession-version']
        self.genus = record['GBSeq_organism'].split(maxsplit=1)[0]
        self.species = " ".join(record['GBSeq_organism'].split()[1:])
        self.kingdom = record['GBSeq_taxonomy'].split(';')[0]


def get_taxonomy(epost_results, args: Namespace) -> None | list[NcbiProtein]:
    """Parse the ePost output from Entrez

    Used when scrapping CAZy and dealing with multiple Taxa.

    :param epost_results: Entrez ePost output
    :param args: cmd-line args parser
    """

    # Retrieve web environment and query key from Entrez epost
    epost_webenv = epost_results["WebEnv"]
    epost_query_key = epost_results["QueryKey"]

    try:
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            db="Protein",
            query_key=epost_query_key,
            WebEnv=epost_webenv,
            retmode="xml",
        ) as record_handle:
            protein_records = Entrez.read(record_handle, validate=False)

    # if no record is returned from call to Entrez
    except (TypeError, AttributeError, RuntimeError,  NotXMLError, IncompleteRead, CorruptedXMLError) as error:
        logger.error(
            "Entrez failed to retireve accession numbers.\n%s\n"
            "Will use the last taxonomy to returned from CAZy", error
        )
        return

    ncbi_data = []
    for record in tqdm(protein_records, desc="Retrieving organism from NCBI"):
        protein = NcbiProtein()
        protein.get_ncbi_data(record)
        ncbi_data.append(protein)

    return ncbi_data


def get_ncbi_taxa(ncbi_acc: list[str], cache_dir: Path, args: Namespace) -> dict[str, int]:
    """Get taxonomy record ids from NCBI"""
    conn = get_sqlite3_connection(args.database)
    create_temp_tax_tables(conn)
    connection_err_cache = []
    invalid_id_cache = []
    tax_stats = {'added': 0, 'updated': 0}

    def process_batch(batch):
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

        # Get the NCBI ids for each protein {id: acc} and dump them all into TEMP_PROTEIN
        err_occurred = fetch_prot_ids_from_entrez(epost_webenv, epost_query_key, conn, args)
        if err_occurred:
            connection_err_cache.append(validated_batch)
            return

        # Get the NCBI ids for the taxonomies, and dump them all into TEMP_TAXONOMY
        err_occurred = get_ncbi_tax_ids(epost_webenv, epost_query_key, conn, args)
        if err_occurred:
            connection_err_cache.append(validated_batch)
            return

    batches = get_chunks_list(ncbi_acc, args.batch_size)
    for batch in tqdm(batches, desc="Querying NCBI"):
        process_batch(batch)

    if args.update: # update existing taxonomy records
        query = """
            SELECT TT.ncbi_tax_id
            FROM TEMP_TAXONOMY TT
            LEFT JOIN NcbiTaxs NT ON TT.ncbi_tax_id = NT.ncbi_tax_id
            WHERE NT.ncbi_tax_id IS NOT NULL
        """
        cur = conn.cursor()
        tax_ids = [row[0] for row in cur.execute(query)]
        tax_batches = get_chunks_list(tax_ids, args.batch_size)
        for tax_batch in tax_batches:
            # ncbi_tax_id is the primary key, so REPLACE overwrites the existing row
            insert_query = """
                INSERT OR REPLACE INTO NcbiTaxs
                (domain, kingdom, phylum, tax_class, tax_order, family, genus, species, strain, ncbi_tax_id)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            err_occurred = get_ncbi_lineages(tax_batch, insert_query, conn, args)
            if err_occurred:
                connection_err_cache.extend(tax_batch)
                continue

    # Add new records
    query = """
        SELECT TT.ncbi_tax_id
        FROM TEMP_TAXONOMY TT
        LEFT JOIN NcbiTaxs NT ON TT.ncbi_tax_id = NT.ncbi_tax_id
        WHERE NT.ncbi_tax_id IS NULL
    """
    cur = conn.cursor()
    tax_ids = [row[0] for row in cur.execute(query)]
    tax_batches = get_chunks_list(tax_ids, args.batch_size)
    for tax_batch in tax_batches:
        insert_query = """
            INSERT INTO NcbiTaxs
            (domain, kingdom, phylum, tax_class, tax_order, family, genus, species, strain, ncbi_tax_id)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        err_occurred = get_ncbi_lineages(tax_batch, insert_query, conn, args)
        if err_occurred:
            connection_err_cache.extend(tax_batch)
            continue

    # update protein links
    # Build the {ncbi_protein_id: protein_accession} mapping needed to translate the numeric
    # protein UIDs that NCBI's elink call returns back into accessions for the UPDATE below.
    cur = conn.cursor()
    cur.execute("SELECT ncbi_protein_id, protein_acccession FROM TEMP_PROTEIN")
    # str()-cast the key: sqlite's INTEGER affinity returns ncbi_protein_id as a Python int,
    # but NCBI's elink response returns protein UIDs as strings - without normalising, the
    # lookup in get_linked_proteins_batch would silently match nothing.
    protein_id_to_acc = {str(row[0]): row[1] for row in cur}

    # Query TEMP_TAXONOMY directly for the full set of tax ids relevant to this run, rather than
    # re-running the "Add new records" query above: by this point those tax ids have already been
    # inserted and committed into NcbiTaxs, so re-running a "WHERE ncbi_tax_id IS NULL" query would
    # now match nothing.
    cur.execute("SELECT ncbi_tax_id FROM TEMP_TAXONOMY")
    all_tax_ids = [row[0] for row in cur]

    tax_id_batches = get_chunks_list(all_tax_ids, args.batch_size)
    for tax_id_batch in tax_id_batches:
        linked_proteins, err_occurred = get_linked_proteins_batch(tax_id_batch, protein_id_to_acc, args)
        if err_occurred:
            connection_err_cache.extend(tax_id_batch)
            continue
        for prot_batch in get_chunks_list(linked_proteins, 500):
            cur.executemany(
                """
                UPDATE Proteins
                SET ncbi_tax_id = ?
                WHERE protein_accession = ?
                """,
                prot_batch
            )

    cur.close()
    conn.commit()

    if connection_err_cache:
        cache_connection_errors(connection_err_cache, cache_dir)

    if invalid_id_cache:
        cache_invalid_ids(invalid_id_cache, cache_dir)

    drop_temp_tax_tables(conn)

    conn.close()

    return tax_stats


def fetch_prot_ids_from_entrez(
    epost_webenv,
    epost_query_key,
    conn,
    args: Namespace
) -> bool:
    """Retrieve protein Ids from NCBI using ePost results
    {protein ncbi id : protein version accession}"""
    err_occurred = False
    records = []

    try:
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            query_key=epost_query_key,
            WebEnv=epost_webenv,
            db="Protein",
            rettype="docsum"
        ) as handle:
            prot_records = Entrez.read(handle, validate=False)
            for record in prot_records:
                records.append((record['Id'], record['AccessionVersion']))
    except Exception as err:
        handle_entrez_errors(err)
        err_occurred = True

    if records:
        cur = conn.cursor()
        batches = get_chunks_list(records, 500)
        for batch in batches:
            cur.executemany(
                """
                INSERT OR IGNORE INTO TEMP_PROTEIN (ncbi_protein_id, protein_acccession)
                VALUES (?, ?)
                """,
                batch
            )
        cur.close()
        conn.commit()

    return err_occurred


def get_ncbi_tax_ids(epost_webenv, epost_query_key, conn, args) -> bool:
    """Link ncbi protein ids to taxonomy ids
    set(ncbi tax ids)
    """
    tax_ids = set()
    err_occurred = False

    try:
        with entrez_retry(
            args.retries,
            Entrez.elink,
            query_key=epost_query_key,
            WebEnv=epost_webenv,
            dbfrom="Protein",
            db="Taxonomy",
            linkname="protein_taxonomy"
        ) as handle:
            tax_links = Entrez.read(handle, validate=False)
            tax_ids = tax_ids.union([
                (link['Id'],)
                for result in tax_links
                for hit in result['LinkSetDb']
                for link in hit['Link']
            ])
    except Exception as err:
        handle_entrez_errors(err)
        err_occurred = True

    if tax_ids:
        cur = conn.cursor()
        batches = get_chunks_list(list(tax_ids), 500)
        for batch in batches:
            cur.executemany(
                """
                INSERT OR IGNORE INTO TEMP_TAXONOMY (ncbi_tax_id)
                VALUES (?)
                """,
                batch
            )
        cur.close()
        conn.commit()

    return  err_occurred


def get_linked_proteins_batch(
    tax_id_batch: list[str],
    prot_id_to_acc: dict[str, str],
    args: Namespace,
) -> tuple[list[tuple], bool]:
    """Retrieve the proteins linked to a batch of taxonomy ids from NCBI.

    Passing multiple ids directly as a list (i.e. as repeated &id= params, not joined into a
    single comma-separated string, and not via epost/WebEnv) makes NCBI's elink return one
    LinkSet per input id, each carrying its own IdList - so the tax_id <-> protein mapping
    survives batching. Joining ids into one string or posting them via epost/WebEnv instead
    merges every input id's links into a single combined LinkSet, losing that mapping.

    Return (list of (ncbi_tax_id, protein_accession) tuples, error occurred flag)
    """
    linked_proteins = []
    err_occurred = False
    try:
        with entrez_retry(
            args.retries,
            Entrez.elink,
            id=tax_id_batch,
            db="Protein",
            dbfrom="Taxonomy",
            linkname="taxonomy_protein",
        ) as handle:
            tax_links = Entrez.read(handle, validate=False)
            for result in tax_links:
                result_tax_ids = result.get('IdList') or []
                if not result_tax_ids:
                    continue
                ncbi_tax_id = result_tax_ids[0]
                for item in result['LinkSetDb']:
                    for link in item['Link']:
                        # link['Id'] is a single id string, not a list - iterating over it
                        # directly would walk individual characters, not the id itself
                        prot_id = link['Id']
                        if prot_id in prot_id_to_acc:
                            linked_proteins.append((ncbi_tax_id, prot_id_to_acc[prot_id]))
    except Exception as err:
        handle_entrez_errors(err)
        err_occurred = True

    return linked_proteins, err_occurred


def get_ncbi_lineages(
    tax_ids: list[str],
    query: str,
    conn: sqlite3.Connection,
    args: Namespace,
    batch_size: int = 500
) -> bool:
    """Retrieve the taxonomy data for a set of taxonomy ids from NCBI

    :param tax_ids: list of taxonomy ids to retrieve taxonomy data for
    :param args: cmd-line args parser
    """
    cur = conn.cursor()
    err_occurred = False
    epost_webenv, epost_query_key = None, None
    records = []
    tax_ids = [str(tax_id) for tax_id in tax_ids]
    try:
        with entrez_retry(
            args.retries,
            Entrez.epost,
            db="Taxonomy",
            id=",".join(tax_ids),
        ) as epost_handle:
            epost_results = Entrez.read(epost_handle, validate=False)
        epost_webenv = epost_results["WebEnv"]
        epost_query_key = epost_results["QueryKey"]
    except (TypeError, AttributeError, RuntimeError,  NotXMLError, IncompleteRead, CorruptedXMLError) as error:
        logger.error(
            "Entrez failed to post NCBI taxonomy accessions.\n%s", error
        )
        err_occurred = True
        return err_occurred

    try:
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            db="Taxonomy",
            query_key=epost_query_key,
            WebEnv=epost_webenv,
            retmode="xml",
        ) as record_handle:
            tax_records = Entrez.read(record_handle, validate=False)
            for record in tax_records:
                domain, kingdom, phylum, tax_class, tax_order, family, genus, species, strain = (None,)*9
                for tax in record['LineageEx']:
                    rank = tax['Rank']
                    name = tax['ScientificName']
                    if rank == 'domain':
                        domain = name
                    elif rank == 'kingdom':
                        kingdom = name
                    elif rank == 'phylum':
                        phylum = name
                    elif rank == 'class':
                        tax_class = name
                    elif rank == 'order':
                        tax_order = name
                    elif rank == 'family':
                        family = name
                    elif rank == 'genus':
                        genus = name
                sci_name = record['ScientificName']
                name_parts = sci_name.split(" ")
                # some taxa (e.g. genus-only or unclassified entries) have a single-word name
                species = name_parts[1] if len(name_parts) > 1 else None
                strain = " ".join(name_parts[2:]) if len(name_parts) > 2 else None                    
                records.append((
                    domain,
                    kingdom,
                    phylum,
                    tax_class,
                    tax_order,
                    family,
                    genus,
                    species,
                    strain,
                    int(record['TaxId'])
                ))
                if len(records) == batch_size:
                    cur.executemany(query, records)
                    records.clear()

        if records:
            cur.executemany(query, records)

    # if no record is returned from call to Entrez
    except (TypeError, AttributeError, RuntimeError,  NotXMLError, IncompleteRead, CorruptedXMLError) as error:
        logger.error(
            "Entrez failed to retrieve Taxonomy data from NCBI\n%s", error
        )
        err_occurred = True
        return err_occurred

    cur.close()
    conn.commit()
