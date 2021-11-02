#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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
"""Add data retrieved from UniProt to a local SQLite database"""


from tqdm import tqdm

from scraper.sql.sql_interface import insert_data
from scraper.sql.sql_orm import Session, Ec


def add_uniprot_accessions(uniprot_dict, gbk_dict, connection, args):
    """Add UniProt data to the local CAZyme database
    
    :param uniprot_dict: dict containing data retrieved from UniProt
    :param gbk_dict: dict representing data from the Genbanks table
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing.
    """
    uniprot_insert_values = set()

    if args.sequence:
        for uniprot_acc in tqdm(uniprot_dict, desc="Adding UniProt data per protein"):
            genbank_acc = uniprot_dict[uniprot_acc]["genbank_accession"]
            gbk_id = gbk_dict[genbank_acc]
            uniprot_name = uniprot_dict[uniprot_acc]["name"]
            seq = uniprot_dict[uniprot_acc]["sequence"]
            date = uniprot_dict[uniprot_acc]["seq_date"]

            uniprot_insert_values.add( (gbk_id, uniprot_acc, uniprot_name, seq, date) )
        
        columns = ['genbank_id', 'uniprot_accession', 'uniprot_name', 'sequence', 'seq_update_date']
    
    else:
        for uniprot_acc in tqdm(uniprot_dict, desc="Adding UniProt data per protein"):
            genbank_acc = uniprot_dict[uniprot_acc]["genbank_accession"]
            gbk_id = gbk_dict[genbank_acc]
            uniprot_name = uniprot_dict[uniprot_acc]["name"]

            uniprot_insert_values.add( (gbk_id, uniprot_acc, uniprot_name) )
        
        columns = ['genbakn_id', 'uniprot_accession', 'uniprot_name']
    
    insert_data(connection, "Uniprots", columns, list(uniprot_insert_values))

    return


def add_ec_numbers(uniprot_dict, all_ecs, gbk_dict, connection):
    """Add EC numbers to the local CAZyme database

    :param uniprot_dict: dict containing data retrieved from UniProt
    :param all_ecs: set of all EC numbers retrieved from UniProt
    :param gbk_dict: dict representing data from the Genbanks table
    :param connection: open sqlalchemy conenction to an SQLite db engine

    Return nothing.
    """
    insert_data(connection, "Ecs", ["ec_number"], list(all_ecs))

    # load Ecs table into memory to retrieve the EC number IDs
    with Session(bind=connection) as session:
        all_db_ecs = session.query(Ec).all()

    # parse Ecs table to create a dict
    ec_dict = {}
    for ec_object in all_db_ecs:
        ec_dict[ec_object.ec_number] = ec_object.ec_id

    gbk_ec_dict = {}
    for uniprot_acc in tqdm(uniprot_dict, desc="Adding EC numbers per protein"):
        genbank_acc = uniprot_dict[uniprot_acc]["genbank_accession"]

        ec_numbers = uniprot_dict[uniprot_acc]["ec"]

        gbk_id = gbk_dict[genbank_acc]

        try:
            gbk_ec_dict[gbk_id]
        except KeyError:
            gbk_ec_dict[gbk_id] = set()

        for ec in ec_numbers:
            ec_id = ec_dict[ec]
            gbk_ec_dict[gbk_id].add(ec_id)
        
    # compile list of tuples to insert into the Genbanks_Ecs table
    gbk_ec_insert_values = set()

    for gbk_id in gbk_ec_dict:
        ec_ids = gbk_ec_dict[gbk_id]
        for ec_id in ec_ids:
            gbk_ec_insert_values.add( (gbk_id, ec_id) )
    
    
    insert_data(connection, "Genbanks_Ecs", ["genbank_id", "ec_id"], list(gbk_ec_insert_values))

    return


def add_pdb_accessions(uniprot_dict, gbk_dict, connection):
    """Add PDB accessions to the local CAZyme database

    :param uniprot_dict: dict containing data retrieved from UniProt
    :param gbk_dict: dict representing data from the Genbanks table
    :param connection: open sqlalchemy conenction to an SQLite db engine

    Return nothing.
    """
    pdb_insert_values = set()

    for uniprot_acc in tqdm(uniprot_dict, desc="Adding PDB accessions per protein"):
        genbank_acc = uniprot_dict[uniprot_acc]["genbank_accession"]
        gbk_id = gbk_dict[genbank_acc]

        pdbs = uniprot_dict[uniprot_acc]["pdb"]
        if len(pdbs) == 0:
            continue
        print(pdbs)

        for pdb in pdbs:
            pdb_insert_values.add(pdb, gbk_id)

    if len(pdb_insert_values) == 0:
        return

    insert_data(connection, "Pdbs", ["pdb_accession", "genbank_id"], list(pdb_insert_values))

    return

