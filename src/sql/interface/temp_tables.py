#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

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

import sqlite3

from pathlib import Path

from src.sql.interface.connect import get_sqlite3_connection


def build_temp_table(db: Path) -> None:
    """Create the TEMP_TABLE matching the SQLAlchemy model definition."""
    conn = get_sqlite3_connection(db)
    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS TEMP_TABLE (
            record_id INTEGER PRIMARY KEY,
            family VARCHAR(255),
            kingdom VARCHAR(255),
            genus VARCHAR(255),
            species VARCHAR(255),
            protein_id VARCHAR(255),
            source VARCHAR(255)
        )
    """)

    # Create the index as defined in the SQLAlchemy model
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS record_id ON TEMP_TABLE (protein_id)
    """)

    conn.commit()
    cursor.close()
    conn.close()


def drop_temp_table(db: Path) -> None:
    """Drop the TEMP_TABLE."""
    conn = get_sqlite3_connection(db)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS TEMP_TABLE")
    cursor.execute("DROP INDEX IF EXISTS record_id")

    conn.commit()
    cursor.close()
    conn.close()


def create_temp_tax_tables(conn: sqlite3.Connection) -> None:
    """Create tables for taxonomy data and mappings"""
    cursor = conn.cursor()

    # Create regular tables instead of temporary ones
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS TEMP_TAXONOMY (
            ncbi_tax_id INTEGER PRIMARY KEY
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS TEMP_PROTEIN (
            ncbi_protein_id INTEGER PRIMARY KEY,
            protein_acccession VARCHAR(50)
        )
    """)

    conn.commit()
    cursor.close()


def drop_temp_tax_tables(conn: sqlite3.Connection) -> None:
    """Drop tables for taxonomy data and mappings"""
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS TEMP_TAXONOMY")
    cursor.execute("DROP TABLE IF EXISTS TEMP_PROTEIN")

    conn.commit()
    cursor.close()


def create_temp_genome_tables(conn: sqlite3.Connection) -> None:
    """Create tables for genome data and mappings"""
    cursor = conn.cursor()

    # Create regular tables instead of temporary ones
    # one row per NCBI Assembly UID, carrying both accessions of that assembly: an
    # assembly has a GenBank (GCA_) and a RefSeq (GCF_) accession and both are needed,
    # since which one a protein is listed under depends on the protein's own accession type
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS TEMP_GENOME (
            ncbi_genome_id INTEGER PRIMARY KEY,
            gbk_accession VARCHAR(50),
            refseq_accession VARCHAR(50),
            assembly_name VARCHAR(255)
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS TEMP_GENOME2PROTEIN (
            ncbi_genome_id INTEGER,
            protein_id INTEGER,
            genome_id INTEGER NULL,
            PRIMARY KEY (ncbi_genome_id, protein_id)
        )
    """)

    conn.commit()
    cursor.close()


def rm_temp_genome_tables(db_path: str) -> None:
    """Remove tables for genome data and mappings"""
    conn = get_sqlite3_connection(db_path)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS TEMP_GENOME")
    cursor.execute("DROP TABLE IF EXISTS TEMP_GENOME2PROTEIN")

    conn.commit()
    cursor.close()
    conn.close()


def create_temp_ec_protein_table(conn: sqlite3.Connection) -> None:
    """Create temporary table for EC to protein mappings"""
    cursor = conn.cursor()

    # Create regular table instead of temporary one
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS TEMP_EC_PROTEIN (
            protein_id INTEGER,
            ec_id INTEGER
        )
    """)

    conn.commit()
    cursor.close()


def drop_temp_ec_protein_table(conn: sqlite3.Connection) -> None:
    """Drop temporary table for EC to protein mappings"""
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS TEMP_EC_PROTEIN")

    conn.commit()
    cursor.close()


def create_temp_pdb_protein_table(conn: sqlite3.Connection) -> None:
    """Create temporary table for PDB to protein mappings"""
    cursor = conn.cursor()

    # Create regular table instead of temporary one
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS TEMP_PDB_PROTEIN (
            protein_id INTEGER,
            pdb_id INTEGER
        )
    """)

    conn.commit()
    cursor.close()


def drop_temp_pdb_protein_table(conn: sqlite3.Connection) -> None:
    """Drop temporary table for PDB to protein mappings"""
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS TEMP_PDB_PROTEIN")

    conn.commit()
    cursor.close()


def create_temp_go_protein_table(conn: sqlite3.Connection) -> None:
    """Create temporary table for GO to protein mappings"""
    conn.commit()
    cursor = conn.cursor()

    # Create regular table instead of temporary one
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS TEMP_GO_PROTEIN (
            protein_id INTEGER,
            go_id INTEGER
        )
    """)

    conn.commit()
    cursor.close()


def drop_temp_go_protein_table(conn: sqlite3.Connection) -> None:
    """Drop temporary table for GO to protein mappings"""
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS TEMP_GO_PROTEIN")

    conn.commit()
    cursor.close()


def create_temp_pfam_protein_table(conn: sqlite3.Connection) -> None:
    """Create temporary table for Pfam-protein match data"""
    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS TEMP_PFAM_PROTEIN (
            protein_id INTEGER,
            pfam_id INTEGER,
            interpro_accession VARCHAR(50),
            match_start INTEGER,
            match_end INTEGER
        )
    """)

    conn.commit()
    cursor.close()


def drop_temp_pfam_protein_table(conn: sqlite3.Connection) -> None:
    """Drop temporary table for Pfam-protein match data"""
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS TEMP_PFAM_PROTEIN")

    conn.commit()
    cursor.close()


def create_temp_pdb_structure_table(conn: sqlite3.Connection) -> None:
    """Create the temp table PDB entry metadata is dumped into as it is retrieved.

    Named ..._STRUCTURE to keep it distinct from TEMP_PDB_PROTEIN, which the UniProt
    pipeline uses for PDB-to-protein relationships.
    """
    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS TEMP_PDB_STRUCTURE (
            pdb_accession VARCHAR(50) PRIMARY KEY,
            method VARCHAR(255),
            resolution FLOAT
        )
    """)

    conn.commit()
    cursor.close()


def drop_temp_pdb_structure_table(conn: sqlite3.Connection) -> None:
    """Drop the temporary PDB entry metadata table"""
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS TEMP_PDB_STRUCTURE")

    conn.commit()
    cursor.close()


def create_temp_gtdb_table(conn: sqlite3.Connection) -> None:
    """Create the temp table GTDB lineages are dumped into as the release files are parsed."""
    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS TEMP_GTDB (
            assembly_accession VARCHAR(50) PRIMARY KEY,
            genome_id INTEGER,
            kingdom VARCHAR(255),
            phylum VARCHAR(255),
            tax_class VARCHAR(255),
            tax_order VARCHAR(255),
            family VARCHAR(255),
            genus VARCHAR(255),
            species VARCHAR(255),
            release VARCHAR(50)
        )
    """)

    conn.commit()
    cursor.close()


def drop_temp_gtdb_table(conn: sqlite3.Connection) -> None:
    """Drop the temporary GTDB lineage table"""
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS TEMP_GTDB")

    conn.commit()
    cursor.close()
