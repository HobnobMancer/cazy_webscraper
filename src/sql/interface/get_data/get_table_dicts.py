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
"""Retrieve all objects from a db table and parse the data to build a dict, repr the current table state."""


import sqlite3

from tqdm import tqdm

from src.sql.sql_orm import (
    Genome,
    GtdbTax,
    Pdb,
    Protein,
    Session,
    Uniprot,
)


def get_gbk_table_seq_dict(connection):
    """Compile a dict of the data in the Proteins table
    
    :param connection: open connection to an SQLite3 database
    
    Return dict {genbank_accession: 'sequence': str, 'seq_date': str}
    """
    with Session(bind=connection) as session:
        all_genbank = session.query(Protein).all()

    db_protein_dict = {}  # {genbank_accession: 'sequence': str, 'seq_date': str}
    
    for gbk in all_genbank:
        db_protein_dict[f"{gbk.genbank_accession}"] = {
            'sequence': gbk.sequence,
            'seq_date': gbk.seq_update_date
        }
    
    return db_protein_dict


def get_pdb_table_dict(connection):
    """Create dict of objects present in the Pdbs table.
    
    :param connection: open sqlalchemy db engine connection
    
    Return dict {pdb_accession: pdb_db_id}
    """
    with Session(bind=connection) as session:
        db_pdb_records = session.query(Pdb).all()

    pdb_table_dict = {}  # {pdb_accession: pdb_db_id}

    for record in tqdm(db_pdb_records, desc="Loading existing PDB db records"):
        pdb_table_dict[record.pdb_accession] = record.pdb_id

    return pdb_table_dict 


def get_gbk_pdb_table_dict(connection):
    """Create dict of objects present in the Proteins_Pdbs table.
    
    :param connection: open sqlalchemy db engine connection
    
    Return dict {gbk_db_id: set(pdb_db_ids) }
    """
    with Session(bind=connection) as session:
        all_gbk_pdb_records = session.query(Protein, Pdb).\
            join(Pdb, Protein.pdbs).\
            all()

    gbk_pdb_table_dict = {}  # {pdb_accession: pdb_db_id}

    for record in tqdm(all_gbk_pdb_records, desc="Loading existing Protein_Pdbs db records"):
        genbank_id = record[0].genbank_id
        pdb_id = record[1].pdb_id

        try:
            gbk_pdb_table_dict[genbank_id].add(pdb_id)
        except KeyError:
            gbk_pdb_table_dict[genbank_id] = {pdb_id}

    return gbk_pdb_table_dict 








def get_uniprot_table_dict(connection):
    """Create dict of objects present in the Uniprots table.
    
    :param connection: open sqlalchemy db engine connection
    
    Return dict {acc: {db_id: int, name: str, seq: str, seq_date:str } }
    """
    with Session(bind=connection) as session:
        db_uniprot_records = session.query(Uniprot).all()

    uniprot_table_dict = {}  # {acc: {name: str, gbk_id: int, seq: str, seq_date:str } }

    for record in tqdm(db_uniprot_records, desc="Retrieving existing UniProt records from db"):
        uniprot_table_dict[record.uniprot_accession] = {
            "db_id": record.uniprot_id,
            "name": record.uniprot_name,
            "seq": record.sequence,
            "seq_date": record.seq_update_date,
        }
    
    return uniprot_table_dict




def get_gtdb_table_dict(connection):
    """Load GTDB table into dict

    :param connection: open connection to an sql db

    Return dict {db id: (tuple of lineage data)}
    """
    with Session(bind=connection) as session:
        query_results = session.query(GtdbTax).all()

    gtdb_dict = {}

    for record in tqdm(query_results, desc="Loading GtdbTax table into dict"):
        gtdb_dict[record.gtdb_tax_id] = (
            record.kingdom,
            record.phylum,
            record.tax_class,
            record.tax_order,
            record.family,
            record.genus,
            record.species,
        )

    return gtdb_dict


def get_genome_table(connection):
    """Load genome table into a dict

    :param connection: open sql db connection

    Return dict {genomic version accession: db id}
    """
    with Session(bind=connection) as session:
        genome_records = session.query(Genome).all()

    db_genome_dict = {}  # {genomic ver acc: {'db_id': db_id, 'gtdb_id': gtdb_id}}

    for record in tqdm(genome_records, desc="Retrieving genome records from the local db"):
        gbk_acc = record.gbk_version_accession
        ref_acc = record.refseq_version_accession
        db_id = record.genome_id
        gtdb_id = record.gtdb_tax_id

        if gbk_acc is not None:
            db_genome_dict[gbk_acc] = {'db_id': db_id, 'gtdb_id': gtdb_id}
        if ref_acc is not None:
            db_genome_dict[ref_acc] = {'db_id': db_id, 'gtdb_id': gtdb_id}

    return db_genome_dict
