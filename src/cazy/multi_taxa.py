#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
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
"""Handle proteins associated with multiple taxa."""


import logging
import sqlite3

from argparse import Namespace
from http.client import IncompleteRead

from Bio.Entrez.Parser import NotXMLError, CorruptedXMLError
from saintBioutils.genbank import entrez_retry
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm

from Bio import Entrez

from cazy_webscraper.ncbi.taxonomy.multiple_taxa import get_ncbi_tax
from cazy_webscraper.sql.sql_interface.add_data.add_ncbi_tax_data import replace_ncbi_taxonomy


logger = logging.getLogger(__name__)


def process_multi_taxa(
    db_path: str,
    multi_taxa_log: str,
    replaced_taxa_log: str,
    args: Namespace
) -> None:
    """Identify cases where a protein ID is associated with multiple taxonomies
    (specific distance genus-species).

    This is not applied to the proteins 
    from the JGI () which are assinged an ID number (int) by CAZy, which 
    does not allow the protein sequence to be tracked back to the entry in 
    the JGI database.

    If skip_ncbi is true, the latest taxonomies are not retrieved from their remote repositories
    and the first taxonomy retrieved from CAZy is used.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Execute the query
    cursor.execute("""
    SELECT protein_id
    FROM TempTable
    WHERE protein_id IN (
        SELECT protein_id
        FROM TempTable
        GROUP BY protein_id
        HAVING COUNT(DISTINCT genus || species) > 1
    ) AND source != 'jgi'
    ORDER BY protein_id;
    """)

    if args.skip_ncbi_tax:
        # delete all rows that have a different taxonomy to
        # first row retrieved from the table
        for row in cursor:
            protein_id = row[0]
            with open(multi_taxa_log, "a") as fh:
                fh.write(f"{protein_id}\n")
            keep_first_taxa(protein_id, conn)

    else:
        protein_ids = [row[0] for row in cursor]
        with open(multi_taxa_log, "a") as fh:
            for protein_id in protein_ids:
                fh.write(f"{protein_id}\n")
        invalid_ids = False  # default to presume all IDs are valid
        success = use_latest_taxa(protein_ids, db_path, args, replaced_taxa_log, conn, invalid_ids)

    # Close the connection
    cursor.close()
    conn.commit()
    conn.close()


def keep_first_taxa(protein_id: str, conn: sqlite3.Connection) -> None:
    """Delete all rows that have a different taxonomy to 
    first row retrieved from the table"""
    cur = conn.cursor()
    cur.execute("""
        SELECT *
        FROM TempTable
        WHERE protein_id = ?;
    """, (protein_id,)
    )

    first_row = cur.fetchone()
    cur.execute("""
        DELETE FROM TempTable
        WHERE protein_id = ? 
            AND (genus != ? OR species != ?);
    """, (protein_id, first_row[1], first_row[2]))

    conn.commit()


def use_latest_taxa(
    protein_ids: list[str],
    db_path: str,
    args: Namespace,
    replaced_taxa_log: str,
    conn: sqlite3.Connection,
    invalid_ids: bool,
) -> None:
    """Retrieve the latest taxonomy from NCBI and use this for
    the taxonomy"""
    batches = get_chunks_list(protein_ids, args.ncbi_batch_size)
    for batch in tqdm(batches, desc=f"Batch retrieving tax info from NCBI. Batch size:{args.ncbi_batch_size}"):
        id_post_list = str(",".join(batch))
        success = False

        try:
            epost_results = Entrez.read(
                entrez_retry(
                    args.retries,
                    Entrez.epost,
                    "Protein",
                    id=id_post_list,
                )
            )
            success = True

        except (TypeError, AttributeError, RuntimeError,  NotXMLError, IncompleteRead, CorruptedXMLError):  # if no record is returned from call to Entrez
            # error not due to the presence of invalid IDs
            logger.error(
                "Entrez failed to post assembly IDs for this batch.\n"
                "Not retrieving tax data from NCBI for these proteins"
                "Selecting the first organism retrieved from CAZy as the source organism\nProtein accessions:\n"
                "%s", batch
            )
            # cazy_data, gbk_accessions, replaced_taxa_logger
            for protein_id in batch:
                keep_first_taxa(protein_id, conn)
            success = True
            continue

        except RuntimeError:
            logger.warning("Found GenBank accessions in CAZy data that are no longer in NCBI")

            if invalid_ids:
                # replace_multiple_tax was called by replace_multiple_tax_with_invalid_ids
                # return results, don't use recursive programming
                continue

            else:
                # first time replace_multiple_tax was called
                replace_multiple_tax_with_invalid_ids(
                    batch,
                    db_path,
                    replaced_taxa_log,
                    args,
                    conn,
                )
                success = True

        if not success:
            logger.error(
                "Could not retrieve taxonomy data from NCBI for this batch,\n"
                "Using the first source organism retrieved from CAZy for each GenBank accession\n"
                "Protein accessions:\n"
                "%s", batch
            )
            for protein_id in batch:
                keep_first_taxa(protein_id, conn)

        else:
            ncbi_data = get_ncbi_tax(epost_results, args)
            if not ncbi_data:
                for protein_id in batch:
                    keep_first_taxa(protein_id, conn)

            for protein in ncbi_data:
                replace_ncbi_taxonomy(protein, conn, replaced_taxa_log)

    return success


def replace_multiple_tax_with_invalid_ids(
    batch: list[str],
    db_path: str,
    replaced_taxa_log: str,
    args: Namespace,
    conn: sqlite3.Connection,
) -> bool:
    """Split up a batch where there may be invalid protein ids, then
    try retrieving data"""
    # retrieve the first half of the list
    mid_point = int((len(batch)/2))
    half_batch = batch[:mid_point]

    success = use_latest_taxa(
        half_batch,
        db_path,
        args,
        replaced_taxa_log,
        conn,
        invalid_ids=True
    )

    if success:
        # invalid IDs are stored in the second half of the accessions list
        half_batch = batch[mid_point:]
        for protein_id in half_batch:
            success = use_latest_taxa(
                [protein_id],
                db_path,
                args,
                replaced_taxa_log,
                conn,
                invalid_ids=True
            )

            if not success:
                keep_first_taxa(protein_id, conn)

    else:
        # invalid IDs are stored in the first half of the accessions list
        for protein_id in half_batch:
            success = use_latest_taxa(
                [protein_id],
                db_path,
                args,
                replaced_taxa_log,
                conn,
                invalid_ids=True
            )

            if not success:
                keep_first_taxa(protein_id, conn)

        # parse the second half of the list
        half_batch = batch[mid_point:]
        success = use_latest_taxa(
            half_batch,
            db_path,
            args,
            replaced_taxa_log,
            conn,
            invalid_ids=True
        )

        if not success:
            # invalid gbk ID present in the second half of the accessions list
            for protein_id in half_batch:
                success = use_latest_taxa(
                    [protein_id],
                    db_path,
                    args,
                    replaced_taxa_log,
                    conn,
                    invalid_ids=True
                )

                if not success:
                    keep_first_taxa(protein_id, conn)
