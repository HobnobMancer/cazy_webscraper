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


from sqlalchemy import delete, text
from tqdm import tqdm

from cazy_webscraper.sql.sql_interface import insert_data, get_table_dicts
from cazy_webscraper.sql.sql_orm import genbanks_ecs


def add_uniprot_accessions(uniprot_dict, gbk_dict, connection, args):
    """Add UniProt data to the local CAZyme database
    
    :param uniprot_dict: dict containing data retrieved from UniProt
    :param gbk_dict: dict representing data from the Genbanks table
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing.
    """
    # load the current Uniprot records in the local CAZyme db
    # {acc: {name: str, gbk_id: int, seq: str, seq_date:str } }
    uniprot_table_dict = get_table_dicts.get_uniprot_table_dict(connection)

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


def add_ec_numbers(uniprot_dict, all_ecs, gbk_dict, connection, args):
    """Add EC numbers to the local CAZyme database

    :param uniprot_dict: dict containing data retrieved from UniProt
    :param all_ecs: set of all EC numbers retrieved from UniProt
    :param gbk_dict: dict representing data from the Genbanks table
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing.
    """
    # load EC records in the local CAZyme db
    ec_table_dict = get_table_dicts.get_ec_table_dict(connection)

    # identify new EC numbers to add to the EC table
    existing_ecs = list(ec_table_dict.keys())
    ec_insert_values = set()

    for ec in all_ecs:
        if ec not in existing_ecs:
            ec_insert_values.add( (ec,) )

    if len(ec_insert_values) != 0:
        insert_data(connection, "Ecs", ["ec_number"], list(ec_insert_values))

        # load in the newly updated EC table from the local CAZyme db
        ec_table_dict = get_table_dicts.get_ec_table_dict(connection)
    
    # load in gbk_ec table, contains the gbk-ec number relationships
    ec_gbk_table_dict = get_table_dicts.get_ec_gbk_table_dict(connection)  # {ec_id: {gbk ids}}

    if args.delete_old_ec:
        gbk_ec_table_dict = get_table_dicts.get_gbk_ec_table_dict(connection)  # {gbk_id: {ec ids}}
        ecs_rel_to_delete = set()  # set of tuples (gbk_id, ec_id)

    # compile list of tuples to insert into the Genbanks_Ecs table
    gbk_ec_insert_values = set()

    for uniprot_acc in tqdm(uniprot_dict, desc="Updating EC numbers"):
        genbank_acc = uniprot_dict[uniprot_acc]["genbank_accession"]
        gbk_id = gbk_dict[genbank_acc]
        
        retrieved_ec_numbers = uniprot_dict[uniprot_acc]["ec"]  # EC#s retrieved from UniProt
        for ec in retrieved_ec_numbers:
            ec_id = ec_table_dict[ec]
            
            existing_gbk_relationships = ec_gbk_table_dict[ec_id]
            if gbk_id not in existing_gbk_relationships:
                gbk_ec_insert_values.add( (gbk_id, ec_id) )
        
        if args.delete_old_ec:
            existing_ec_relationships = gbk_ec_table_dict[gbk_id]
            for ec in retrieved_ec_numbers:
                ec_id = ec_table_dict[ec]
                if ec_id not in existing_ec_relationships:
                    ecs_rel_to_delete.add( (gbk_id, ec_id) )

    if len(gbk_ec_insert_values) != 0:
        insert_data(connection, "Genbanks_Ecs", ["genbank_id", "ec_id"], list(gbk_ec_insert_values))

    if args.delete_old_ec and len(ecs_rel_to_delete) != 0:
        for record in tqdm(ecs_rel_to_delete, desc="Deleteing old GenBank-EC relationships"):
            # record = (genbank_id, ec_id,)
            stmt = (
                delete(genbanks_ecs).\
                where(genbanks_ecs.c.genbank_id == record[0]).\
                where(genbanks_ecs.c.ec_id == record[1])
            )
            connection.execute(stmt)

    return


def add_pdb_accessions(uniprot_dict, gbk_dict, connection, args):
    """Add PDB accessions to the local CAZyme database

    :param uniprot_dict: dict containing data retrieved from UniProt
    :param gbk_dict: dict representing data from the Genbanks table
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing.
    """
    # load in PDB objects in the local CAZyme db
    pdb_table_dict = get_table_dicts.get_pdb_table_dict(connection)

    pdb_insert_values = set()  # new records to add to db
    update_pdbs = set()  # records to update
    existing_pdbs = set(pdb_table_dict.values())  # records already in the db
    if args.delete_old_pdbs:
        delete_pdbs = set()  # records to delete, PDB acc no longer listed in UniProt for the given protein

    for uniprot_acc in tqdm(uniprot_dict, desc="Identifying new and updated PDB records"):
        genbank_acc = uniprot_dict[uniprot_acc]["genbank_accession"]
        gbk_id = gbk_dict[genbank_acc]

        retrieved_pdbs = uniprot_dict[uniprot_acc]["pdb"]
        if len(retrieved_pdbs) == 0:
            continue

        for pdb in retrieved_pdbs:
            if pdb_acc in existing_pdbs:
                # check if the gbk_id is the same
                current_pdbs = pdb_table_dict[gbk_id]
                if pdb_acc not in current_pdbs:
                    update_pdbs.add( (pdb_acc, gbk_id) )
            
            else:
                pdb_insert_values.add( (pdb_acc, gbk_id) )
        
        if args.delete_old_pdbs:
            existing_pdb_records = pdb_table_dict[gbk_id]
            for pdb_acc in existing_pdb_records:
                if pdb_acc not in retrieved_pdbs:
                    delete_pdbs.add( (pdb_acc,) )

    if len(pdb_insert_values) != 0:
        insert_data(connection, "Pdbs", ["pdb_accession", "genbank_id"], list(pdb_insert_values))

    if len(update_pdbs) != 0:
        for record in tqdm(update_pdbs, desc="Updating PDB records"):
            connection.execute(
                text(
                    "UPDATE Pdbs "
                    f"SET genbank_id = {record[1]} "
                    f"WHERE pdb_accession = '{record[0]}'"
                )
            )

    if args.delete_old_pdbs and len(delete_pdbs) != 0:
        for record in tqdm(delete_pdbs, desc="Deleteing old PDB accessions"):
            # record = (pdb_acc,)
            connection.execute(
                text(
                    "DELETE Pdbs "
                    f"WHERE pdb_accession = '{record[0]}'"
                )
            )
            
    return
