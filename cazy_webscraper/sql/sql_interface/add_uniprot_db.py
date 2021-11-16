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


from sqlalchemy import text
from tqdm import tqdm

from cazy_webscraper.sql.sql_interface import insert_data
from cazy_webscraper.sql.sql_orm import Session, Ec


def add_uniprot_accessions(uniprot_dict, gbk_dict, uniprot_table_dict, connection, args):
    """Add UniProt data to the local CAZyme database
    
    :param uniprot_dict: dict containing data retrieved from UniProt
    :param gbk_dict: dict representing data from the Genbanks table
    :param uniprot_table_dict: dict representing the Uniprots table
        {acc: {name: str, gbk_id: int, seq: str, seq_date:str } }
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing.
    """
    uniprot_insert_values = set()

    # the following variables containing objects in the local db that are to be updated
    update_record_gbk_id = set()  # ((uniprot_accession, gbk_id),)
    update_record_name = set()  # ((uniprot_accession, retrieved_name), )
    update_record_seq = set()  # ((uniprot_accession, retrieved_seq), )

    for uniprot_acc in tqdm(uniprot_dict, desc="Separating new and existing records"):
        # check if the UniProt accession is already in the local CAZyme db
        try:
            uniprot_table_dict[uniprot_acc]
            # Uniprot accession is already in the local CAZyme db

            # check if the GenBank id is the same
            existing_gbk_id = uniprot_table_dict[uniprot_acc]['genbank_id']
            retrieved_gbk_acc = uniprot_dict[uniprot_acc]["genbank_accession"]
            retrieved_gbk_id = gbk_dict[retrieved_gbk_acc]
            if existing_gbk_id != retrieved_gbk_id:
                update_record_gbk_id.add( (uniprot_acc, retrieved_gbk_id) )
                
            if args.update_name:
                # check if the name has changed
                existing_name = uniprot_table_dict[uniprot_acc]['name']
                retrieved_name = uniprot_dict[uniprot_acc]['name']
                
                if existing_name != retrieved_name:
                    update_record_name.add( (uniprot_acc, retrieved_name) )
            
            if args.sequence_new:  # only add seq if sequence is not there
                if uniprot_table_dict[uniprot_acc]['seq'] is None:
                    update_record_seq.add( (uniprot_acc, uniprot_dict[uniprot_acc]["sequence"]) )
                
            if args.sequence_update:  # update seq if a newer version is available or no seq present
                if uniprot_table_dict[uniprot_acc]['seq'] is None:
                        update_record_seq.add( (
                            uniprot_acc,
                            uniprot_dict[uniprot_acc]["sequence"],
                            uniprot_dict[uniprot_acc]["seq_date"],
                        ) )
                else:
                    existing_date = uniprot_table_dict[uniprot_acc]['seq_date']
                    retrieved_date = uniprot_dict[uniprot_acc]['seq_date']
                    
                    if existing_date < retrieved_date:  # existing date is older, update seq
                        update_record_seq.add( (
                            uniprot_acc,
                            uniprot_dict[uniprot_acc]["sequence"],
                            uniprot_dict[uniprot_acc]["seq_date"],
                        ) )
        
        except KeyError:  # new record to add to local CAZyme db
            
            if args.sequence_new or args.sequence_update:
                genbank_acc = uniprot_dict[uniprot_acc]["genbank_accession"]
                gbk_id = gbk_dict[genbank_acc]
                uniprot_name = uniprot_dict[uniprot_acc]["name"]
                seq = uniprot_dict[uniprot_acc]["sequence"]
                date = uniprot_dict[uniprot_acc]["seq_date"]

                uniprot_insert_values.add( (gbk_id, uniprot_acc, uniprot_name, seq, date) )

            else:  # not retrieving protein sequences
                genbank_acc = uniprot_dict[uniprot_acc]["genbank_accession"]
                gbk_id = gbk_dict[genbank_acc]
                uniprot_name = uniprot_dict[uniprot_acc]["name"]
                
                uniprot_insert_values.add( (gbk_id, uniprot_acc, uniprot_name) )
    
    if len(uniprot_insert_values) != 0:
        if args.sequence_new or args.sequence_update:
            columns = ['genbank_id', 'uniprot_accession', 'uniprot_name', 'sequence', 'seq_update_date']
        else:
            columns = ['genbakn_id', 'uniprot_accession', 'uniprot_name']
    
        insert_data(connection, "Uniprots", columns, list(uniprot_insert_values))

    if len(update_record_gbk_id) != 0:
        for record in tqdm(update_record_gbk_id, desc="Updating UniProt-Gbk relationships"):
            connection.execute(
                text(
                    "UPDATE Uniprots "
                    f"SET genbank_id = {record[1]} "
                    f"WHERE uniprot_accession = '{record[0]}'"
                )
            )
        
    if len(update_record_name) != 0:
        for record in tqdm(update_record_name, desc="Updating UniProt protein names"):
            connection.execute(
                text(
                    "UPDATE Uniprots "
                    f"SET uniprot_name = {record[1]} "
                    f"WHERE uniprot_accession = '{record[0]}'"
                )
            )

    if len(update_record_seq) != 0:
        for record in tqdm(update_record_seq, desc="Updating UniProt protein seqs"):
            connection.execute(
                text(
                    "UPDATE Uniprots "
                    f"SET sequence = {record[1]}, seq_update_date = {record[2]} "
                    f"WHERE uniprot_accession = '{record[0]}'"
                )
            )

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

        for pdb in pdbs:
            pdb_insert_values.add(pdb, gbk_id)

    if len(pdb_insert_values) == 0:
        return

    insert_data(connection, "Pdbs", ["pdb_accession", "genbank_id"], list(pdb_insert_values))

    return

