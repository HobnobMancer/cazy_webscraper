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
"""Add CAZyme data to a local SQLite database"""


import logging
import sqlite3

from pathlib import Path

from sqlalchemy import delete, text
from tqdm import tqdm

from cazy_webscraper.sql.sql_interface import insert_data
from cazy_webscraper.sql.sql_interface.get_data.get_table_dicts import (
    get_kingdom_table_dict,
    get_taxs_table_dict,
    get_fams_table_dict,
    get_protein_table_dict,
    get_prot_fam_table_dict,
)
from cazy_webscraper.sql.sql_orm import proteins_families


logger = logging.getLogger(__name__)


def add_kingdoms(db: Path) -> None:
    """Add new Kingdoms objects to database from the TempTable

    Check existing kingdom objects in the db against kingdoms retrieved from the 
    CAZy txt file, so as to only add new kingdoms.
    """
    conn = sqlite3.connect(db)
    kingdom_table_dict = get_kingdom_table_dict(conn)  # dict {kingdom: kingdom_id}

    # retrieve the Kingdoms retrieved from the CAZy txt file
    cur = conn.cursor()
    cur.execute("""SELECT DISTINCT(kingdom) FROM TempTable""")
    # create list of tuples for db insert
    kingdoms_db_insert_values = [(row[0],) for row in cur if row[0] not in kingdom_table_dict]
    cur.close()

    if len(kingdoms_db_insert_values) != 0:
        insert_data(conn, 'Kingdoms', ['kingdom'], kingdoms_db_insert_values)

    conn.commit()
    conn.close()


def add_source_organisms(db: Path) -> None:
    """Add taxonomy (source organism) data from TempTable to the Tax table in
    the local CAZyme database

    Check existing kingdom objects in the db against kingdoms retrieved from the 
    CAZy txt file, so as to only add new kingdoms.
    """
    conn = sqlite3.connect(db)
    # get the kingdom_id for the Taxs table
    kingdom_table_dict = get_kingdom_table_dict(conn)  # dict {kingdom: kingdom_id}
    # {genus species: {'tax_id': db_tax_id, 'kingdom_id': kingdom_id}
    tax_table_dict = get_taxs_table_dict(conn)
    # {genus species: {'tax_id': int(db_tax_id), 'kingdom_id': int(kingdom_id)}

    tax_insert_values = set()  # new taxa to add to db (kingdom_id, genus, species)
    records_to_update = set()  # used incase kingdom has changed for a species
    cur = conn.cursor()
    cur.execute("""SELECT DISTINCT kingdom, genus, species FROM TempTable""")
    for row in cur:
        if f"{row[1]} {row[2]}" in tax_table_dict:
            if int(row[0]) != int(tax_table_dict[f"{row[1]} {row[2]}"]['kingdom_id']):
                records_to_update.add((row[0], row[1], row[2],))
        else:
            tax_insert_values.add((row[1], row[2], kingdom_table_dict[row[0]]))

    if len(tax_insert_values) != 0:
        logger.warning("Inserting %s new taxs into the local db", len(tax_insert_values))
        insert_data(conn, 'Taxs', ['genus', 'species', 'kingdom_id'], list(tax_insert_values))

    if len(records_to_update) != 0:
        logger.warning("Updating %s tax records in th local CAZyme db", len(records_to_update))
        for record in records_to_update:
            conn.execute(
                """UPDATE Taxs SET kingdom_id = ? WHERE genus = ? AND species = ?""",
                record
            )

    conn.commit()
    conn.close()


def add_cazy_families(db: Path) -> None:
    """Add CAZy families and subfamilies to local CAZyme database"""
    conn = sqlite3.connect(db)
    fam_table_dict = get_fams_table_dict(conn)  # {'<fam> <subfam|_>': db fam id}

    families_db_insert_values = set()  # new fam records to add to db
    cur = conn.cursor()
    cur.execute("""SELECT DISTINCT family FROM TempTable""")
    for row in cur:
        fam, subfam = (row[0].split('_')[0], row[0]) if '_' in row[0] else (row[0], None)

        if f"{fam} {subfam}" not in fam_table_dict:
            families_db_insert_values.add((fam, subfam if subfam != '_' else None))

    if len(families_db_insert_values) != 0:
        logger.warning(
            "Inserting %s new families-subfamilies into the local db",
            len(families_db_insert_values)
        )
        insert_data(conn, 'CazyFamilies', ['family', 'subfamily'], list(families_db_insert_values))

    conn.commit()
    conn.close()


def add_proteins(db: Path) -> None:
    """Add GenBank accs from TempTable to Proteins table"""
    conn = sqlite3.connect(db)
    protein_table_dict = get_protein_table_dict(conn)
    # {prot acc: {'taxa_id': str, 'id': int}}
    taxa_table_dict = get_taxs_table_dict(conn)
    # {genus species: {'tax_id': int(db_tax_id), 'kingdom_id': int(kingdom_id)}

    prot_record_updates = set()  # {prot_accession: 'taxa_id': (new taxa_id) int, 'protein_id': int}
    prot_db_insert_values = set()

    cur = conn.cursor()
    cur.execute("""SELECT DISTINCT protein_id, source, genus, species FROM TempTable""")
    for row in cur:
        protein_acc = row[0]
        tax_id = taxa_table_dict[f"{row[2]} {row[3]}"]['tax_id']
        if protein_acc not in protein_table_dict:
            prot_db_insert_values.add((protein_acc, tax_id, row[1]))
        elif tax_id != protein_table_dict[protein_acc]['taxa_id']:
            prot_record_updates.add((tax_id, protein_table_dict[protein_acc]['protein_id']))

    if len(prot_db_insert_values) != 0:
        logger.warning("Inserting %s Protein accessions into the db", len(prot_db_insert_values))
        insert_data(
            conn,
            'Proteins',
            ['protein_accession', 'taxonomy_id', 'source'],
            list(prot_db_insert_values)
        )

    if len(prot_record_updates) != 0:
        logger.warning(
            "Updating %s Protein table records with new taxonomy IDs",
            len(prot_record_updates)
        )
        for record in prot_record_updates:
            conn.execute(
                """UPDATE Proteins SET taxonomy_id = ? WHERE protein_id = ?""",
                record
            )

    conn.commit()
    conn.close()


def add_protein_fam_relationships(db: Path) -> None:
    """Add Protein accession - CAZy family relationships to db"""
    insert_values = set()  # new records to add
    records_to_del = set()  # records/relationships to delete

    conn = sqlite3.connect(db)
    fam_table_dict = get_fams_table_dict(conn)  # {'fam subfam': fam db id}
    prot_fam_table_dict = get_prot_fam_table_dict(conn)  # {protein_id: {family_id}}

    cur = conn.cursor()
    cur.execute("""SELECT DISTINCT protein_id, family FROM TempTable""")

    for row in cur:
        protein_acc = row[0]
        fam, subfam = (row[1].split('_')[0], row[1]) if '_' in row[1] else (row[1], 'None')
        fam_id = fam_table_dict[f"{fam} {subfam}"]

        prot_cur = conn.cursor()
        prot_cur.execute(
            "SELECT protein_id FROM Proteins WHERE protein_accession = ?",
            (protein_acc,)
        )

        result = prot_cur.fetchone()
        if result:
            protein_id = result[0]
        else:
            logger.error((
                "Could not find protein id for accession %s in the Proteins table\n"
                "Not adding fams for this protein"
            ))
            continue
        prot_cur.close()

        if protein_id not in prot_fam_table_dict:
            insert_values.add((protein_id, fam_id))
        elif fam_id not in prot_fam_table_dict[protein_id]:
            records_to_del.add((protein_id, fam_id))

    cur.close()

    if len(insert_values) != 0:
        logger.warning("Inserting %s protein-family relationships into the db", len(insert_values))
        insert_data(
            conn,
            'Proteins_CazyFamilies',
            ['protein_id', 'family_id'],
            list(insert_values),
        )

    if (len(records_to_del) != 0):
        logger.warning("Deleting %s defunct protein-family relationships", len(records_to_del))
        del_cur = conn.connect()
        for record in records_to_del:
            del_cur.execute(
                "DELETE FROM Proteins_CazyFamilies WHERE protein_id = ? AND family_id = ?",
                (record[0], record[1])
            )
            conn.commit()
        del_cur.close()

    conn.commit()
    conn.close()

    return
