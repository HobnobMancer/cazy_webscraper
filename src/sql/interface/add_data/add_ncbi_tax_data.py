#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
#
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
"""Add NCBI Taxonomy data to a local SQLite database"""


import logging
import sqlite3

from sqlalchemy import text
from tqdm import tqdm

from src.sql.interface.get_data.get_table_dicts import get_ncbi_tax_table
from src.databases.ncbi.taxonomies import NcbiProtein


logger = logging.getLogger(__name__)


def update_ncbi_taxs(tax_dict, db, args):
    """Add NCBI taxonomy data to the local db using sqlite3 commands

    :param tax_dict: dict, keyed by ncbi tax db id, and valued by lineage data
    :param connection: open sqlite3.Connection to a db
    :param args: cmd-line args parser

    Return nothing
    """
    conn = sqlite3.connect(db)
    cur = conn.cursor()

    # Load ncbiTax table into dict
    db_ncbi_tax_table = get_ncbi_tax_table(conn)  # {ncbi_tax_id: local db id}

    records_to_add = set()
    records_to_update = set()

    for ncbi_tax_id in tax_dict:
        tax_data = (
            ncbi_tax_id,
            tax_dict[ncbi_tax_id]['kingdom'],
            tax_dict[ncbi_tax_id]['phylum'],
            tax_dict[ncbi_tax_id]['class'],
            tax_dict[ncbi_tax_id]['order'],
            tax_dict[ncbi_tax_id]['family'],
            tax_dict[ncbi_tax_id]['genus'],
            tax_dict[ncbi_tax_id]['species'],
            tax_dict[ncbi_tax_id]['strain'],
        )

        try:
            db_ncbi_tax_table[int(ncbi_tax_id)]
            if args.update_taxs:
                records_to_update.add(tax_data)
        except KeyError:
            records_to_add.add(tax_data)

    if args.update_taxs and len(records_to_update) != 0:
        # update tax records using sqlite3
        for record in records_to_update:
            cur.execute(
                """
                UPDATE NcbiTaxs
                SET kingdom = ?,
                    phylum = ?,
                    tax_class = ?,
                    tax_order = ?,
                    family = ?,
                    genus = ?,
                    species = ?,
                    strain = ?
                WHERE ncbi_tax_id = ?
                """,
                (
                    record[1], record[2], record[3], record[4],
                    record[5], record[6], record[7], record[8], record[0]
                )
            )
        conn.commit()

    if len(records_to_add) != 0:
        # add records using sqlite3
        cur.executemany(
            """
            INSERT INTO NcbiTaxs (
                ncbi_tax_id,
                kingdom,
                phylum,
                tax_class,
                tax_order,
                family,
                genus,
                species,
                strain
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            list(records_to_add)
        )
        conn.commit()

    cur.close()
    return


def replace_ncbi_taxonomy(
    protein: NcbiProtein,
    conn: sqlite3.Connection,
    replaced_taxa_log: str,
) -> None:
    """Replace all taxa for the protein id with the taxa retrieved from NCBI"""
    cur = conn.cursor()

    cur.execute("""
        UPDATE TempTable
        SET genus = ?,
            species = ?,
            kingdom = ?
        WHERE protein_id = ?;
    """, (protein.genus, protein.species, protein.kingdom, protein.protein_id))

    message = (
        f"{protein.protein_id}\t"
        f"\tREPLACED TAXA WITH:\t{protein.kingdom}\t{protein.genus}:\t{protein.species}"
    )
    with open(replaced_taxa_log, "a") as fh:
        fh.write(f"{message}\n")
        logger.warning(message)

    cur.close()
    conn.commit()
