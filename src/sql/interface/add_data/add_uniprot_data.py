#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2025
# (c) University of Strathclyde 2025
# (c) James Hutton Institute 2025
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
"""Add data retrieved from UniProt to a local SQLite database"""


import logging
import sqlite3

from pathlib import Path

from src.sql.interface import insert_data
from src.sql.interface.get_data.get_uniprot_data import (
    get_uniprot_summary, get_uniprot_to_id,
    get_db_ecs, get_ec_to_id,
    get_db_pdbs, get_pdb_to_id,
    get_db_gos, get_go_to_id
)
from src.sql.interface.temp_tables import (
    create_temp_ec_protein_table,
    create_temp_pdb_protein_table,
    create_temp_go_protein_table
)

logger = logging.getLogger(__name__)


def add_uniprot_records(records: list, db: Path, update_seq: bool, batch_size: int = 500) -> tuple[int, int]:
    """Add UniProt data to the local CAZyme database

    :param records: list of UniProtRecord objects
    :param db: path to the local CAZyme database
    :param update_seq: whether to update existing sequences in the db
    :param batch_size: number of records to add in each batch
    """
    conn = sqlite3.connect(db)
    records_to_update = []
    records_to_add = set()  # There are a few cases of multiple ncbis --> 1 uniprot
    db_table = get_uniprot_summary(conn)  # uniprot_acc: bool if it has a seq

    for record in records:
        if record.uniprot_id not in db_table:
            records_to_add.add(
                (
                    record.uniprot_id,
                    record.uniparc_id,
                    record.name,
                    record.swissprot,
                    record.sequence,
                    record.md5,
                    record.mol_weight,
                    record.sequence_date,
                )
            )
        elif not db_table[record.uniprot_id] or update_seq:
            records_to_update.append(
                (
                    record.sequence,
                    record.md5,
                    record.mol_weight,
                    record.sequence_date,
                    record.uniprot_id,
                )
            )

    if records_to_add:
        insert_data(conn, "Uniprots", ["uniprot_accession", "uniparc_id", "name", "swissprot", "sequence", "md5", "molecular_weight", "seq_update_date"], list(records_to_add))

    cursor = conn.cursor()
    for i in range(0, len(records_to_update), batch_size):
        batch = records_to_update[i:i + batch_size]
        cursor.executemany(
            "UPDATE Uniprots SET sequence = ?, md5 = ?, molecular_weight = ?, seq_update_date = ? WHERE uniprot_id = ?",
            batch
        )
        logger.debug("Processed batch %d-%d of %d records for uniprot seq update", i + 1, min(i + batch_size, len(records_to_update)), len(records_to_update))

    uniprot2dbid = get_uniprot_to_id(conn, set([record.uniprot_id for record in records]))
    protein_records_updates = [(uniprot2dbid[record.uniprot_id], record.ncbi_acc) for record in records]
    for i in range(0, len(protein_records_updates), batch_size):
        batch = protein_records_updates[i:i + batch_size]
        cursor.executemany(
            "UPDATE Proteins SET uniprot_id = ? WHERE protein_accession = ?",
            batch
        )
        logger.debug("Processed batch %d-%d of %d records for ncbi-uniprot update", i + 1, min(i + batch_size, len(protein_records_updates)), len(protein_records_updates))

    conn.commit()
    cursor.close()

    return len(records_to_add), len(records_to_update)


def add_ec_numbers(records: list, conn: sqlite3.Connection, batch_size: int = 500) -> int:
    db_ecs = get_db_ecs(conn)
    records_to_add = set()
    all_ecs = set()
    for record in records:
        for ec in record.ec_nums:
            if ec not in db_ecs:
                records_to_add.add( (ec,) )
            all_ecs.add(ec)

    records_to_add = list(records_to_add)
    for i in range(0, len(records_to_add), batch_size):
        batch = records_to_add[i:i + batch_size]
        insert_data(conn, "Ecs", ["ec_number"], batch)
        logger.debug("Processed batch %d-%d of %d EC numbers", i + 1, min(i + batch_size, len(records_to_add)), len(records_to_add))

    ec2id = get_ec_to_id(conn, all_ecs)
    if ec2id:
        all_relationships = set()
        for record in records:
            for ec in record.ec_nums:
                ec_id = ec2id[ec]
                all_relationships.add( (record.protein_id, ec_id) )

        create_temp_ec_protein_table(conn)
        insert_data(conn, "TEMP_EC_PROTEIN", ["protein_id", "ec_id"], list(all_relationships))

    conn.commit()

    return len(records_to_add)


def add_pdbs(records: list, conn: sqlite3.Connection, batch_size: int = 500) -> int:
    """Add PDB accessions to the database"""
    db_pdbs = get_db_pdbs(conn)
    pdbs_to_update = []  # if the resolution has changed
    pdbs_to_add = []
    all_pdbs = set()
    for record in records:
        for pdb_acc, method, resolution in record.pdbs:
            all_pdbs.add(pdb_acc)
            if pdb_acc not in db_pdbs:
                pdbs_to_add.append((pdb_acc, method, resolution))
            elif resolution != db_pdbs[pdb_acc]:
                pdbs_to_update.append((pdb_acc, method, resolution))

    cursor = conn.cursor()
    for i in range(0, len(pdbs_to_update), batch_size):
        batch = pdbs_to_update[i:i + batch_size]
        cursor.executemany(
            "UPDATE PDBs SET resolution = ? WHERE pdb_accession = ?",
            batch
        )
    conn.commit()
    cursor.close()

    for i in range(0, len(pdbs_to_add), batch_size):
        batch = pdbs_to_add[i:i + batch_size]
        insert_data(conn, "PDBs", ["pdb_accession", "method", "resolution"], batch)

    pdb2id = get_pdb_to_id(conn, all_pdbs)
    if pdb2id:
        all_relationships = set()
        for record in records:
            for pdb_ac in all_pdbs:
                pdb_id = pdb2id[pdb_ac]
                all_relationships.add( (record.protein_id, pdb_id) )

        create_temp_pdb_protein_table(conn)
        insert_data(conn, "TEMP_PDB_PROTEIN", ["protein_id", "pdb_id"], list(all_relationships))

    conn.commit()

    return len(pdbs_to_add)


def add_go_terms(records: list, conn: sqlite3.Connection, batch_size: int = 500) -> int:
    """Add GO terms to the database"""

    db_gos = get_db_gos(conn)
    gos_to_add = []
    all_gos = set()
    for record in records:
        for goterm_id, desc in record.go_terms.items():
            all_gos.add(goterm_id)
            if goterm_id not in db_gos:
                gos_to_add.append( (goterm_id, desc) )

    for i in range(0, len(gos_to_add), batch_size):
        batch = gos_to_add[i:i + batch_size]
        insert_data(conn, "GoTerms", ["goterm_id", "description"], batch)

    go2id = get_go_to_id(conn, all_gos)
    if go2id:
        all_relationships = set()
        for record in records:
            for go_id in all_gos:
                go_db_id = go2id[go_id]
                all_relationships.add( (record.protein_id, go_db_id) )

        create_temp_go_protein_table(conn)
        insert_data(conn, "TEMP_GO_PROTEIN", ["protein_id", "go_id"], list(all_relationships))

    conn.commit()

    return len(gos_to_add)


def merge_temp_ec_relationships(db_path: str) -> int:
    """Merge EC-protein relationships from TEMP_EC_PROTEIN into Proteins_Ecs

    Only inserts relationships that don't already exist in Proteins_Ecs.

    :param db_path: path to the database
    :return: number of new relationships added
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Insert relationships that don't already exist
    cursor.execute("""
        DELETE FROM TEMP_EC_PROTEIN
        WHERE rowid NOT IN (
            SELECT MIN(rowid)
            FROM TEMP_EC_PROTEIN
            GROUP BY protein_id, ec_id
        );
    """)
    cursor.execute("""
        INSERT OR IGNORE INTO Proteins_Ecs (protein_id, ec_id)
        SELECT protein_id, ec_id
        FROM TEMP_EC_PROTEIN
    """)

    conn.commit()
    cursor.close()
    conn.close()

    return


def merge_temp_pdb_relationships(db_path: str) -> int:
    """Merge PDB-protein relationships from TEMP_PDB_PROTEIN into Proteins_PDBs

    Only inserts relationships that don't already exist in Proteins_PDBs.

    :param db_path: path to the database
    :return: number of new relationships added
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Insert relationships that don't already exist
    cursor.execute("""
        DELETE FROM TEMP_PDB_PROTEIN
        WHERE rowid NOT IN (
            SELECT MIN(rowid)
            FROM TEMP_PDB_PROTEIN
            GROUP BY protein_id, pdb_id
        );
    """)
    cursor.execute("""
        INSERT OR IGNORE INTO Proteins_PDBs (protein_id, pdb_id)
        SELECT protein_id, pdb_id
        FROM TEMP_PDB_PROTEIN
    """)

    conn.commit()
    cursor.close()
    conn.close()

    return


def merge_temp_go_relationships(db_path: str) -> int:
    """Merge GO-protein relationships from TEMP_GO_PROTEIN into Proteins_GOs

    Only inserts relationships that don't already exist in Proteins_GOs.

    :param db_path: path to the database
    :return: number of new relationships added
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Insert relationships that don't already exist
    cursor.execute("""
        DELETE FROM TEMP_GO_PROTEIN
        WHERE rowid NOT IN (
            SELECT MIN(rowid)
            FROM TEMP_GO_PROTEIN
            GROUP BY protein_id, go_id
        );
    """)
    cursor.execute("""
        INSERT OR IGNORE INTO Proteins_GoTerms (protein_id, go_id)
        SELECT protein_id, go_id
        FROM TEMP_GO_PROTEIN
    """)

    conn.commit()
    cursor.close()
    conn.close()

    return
