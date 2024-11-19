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
"""Add data retrieved from UniProt to a local SQLite database"""


import logging
from sqlalchemy import delete, text
from tqdm import tqdm

from cazy_webscraper.sql.sql_interface import insert_data
from cazy_webscraper.sql.sql_interface.get_data import get_table_dicts
from cazy_webscraper.sql.sql_orm import genbanks_ecs, Ec, Pdb


def add_uniprot_accessions(uniprot_dict, connection, args):
    """Add UniProt data to the local CAZyme database
    
    :param uniprot_dict: dict containing data retrieved from UniProt
        uniprot_dict[ncbi_acc] = {
            'uniprot_acc': uniprot_acc, - str
            'protein_name': protein_name, -str
            'ec_numbers': ec_numbers, - set
            'sequence': sequence, - str
            'sequence_date': date seq was last updated yyyy-mm-dd
            'pdbs': all_pdbs, - set
        }
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # load the current Uniprot records in the local CAZyme db
    # {acc: {name: str, gbk_id: int, seq: str, seq_date:str } }
    uniprot_table_dict = get_table_dicts.get_uniprot_table_dict(connection)

    uniprot_insert_values = set()   # new rows to add to the local CAZyme db Uniprots table
    record_names_to_update = set()  # ((uniprot_accession, retrieved_name), )
    record_seqs_to_update = set()   # ((uniprot_accession, retrieved_sequence), )

    for ncbi_acc in tqdm(uniprot_dict, desc="Separating new and existing Uniprots Table records"):
        # check if the uniprot accession is in the local CAZyme datbase
        uniprot_acc = uniprot_dict[ncbi_acc]['uniprot_acc']

        try:
            uniprot_table_dict[uniprot_acc]
        except KeyError:
            # uniprot acc not in Uniprots table, so will need to add
            if args.sequence:  # add seq to local db
                uniprot_insert_values.add(
                    (
                        uniprot_acc,
                        uniprot_dict[ncbi_acc]['protein_name'],
                        uniprot_dict[ncbi_acc]['sequence'],
                        uniprot_dict[ncbi_acc]['sequence_date'],
                    )
                )
            else:  # do not add sequence
                uniprot_insert_values.add(
                    (
                        uniprot_acc,
                        uniprot_dict[ncbi_acc]['protein_name'],
                    )
                )
            continue

        # uniprot acc already in the db
        # check if need to update name or sequence

        if args.update_name:
            if uniprot_table_dict[uniprot_acc]['name'] != uniprot_dict[ncbi_acc]['protein_name']:
                logger.warning(
                    f"Updating protein name for UniProt acc {uniprot_acc} from:\n"
                    f"{uniprot_table_dict[uniprot_acc]['name']}\n"
                    "to:\n"
                    f"{uniprot_dict[ncbi_acc]['protein_name']}"
                )
                record_names_to_update.add(
                    (
                        uniprot_acc,
                        uniprot_dict[ncbi_acc]['protein_name'],
                    )
                )

        if args.update_seq:
            if uniprot_table_dict[uniprot_acc]['seq'] != uniprot_dict[ncbi_acc]['sequence']:
                logger.warning(
                    f"Updating protein sequence for UniProt acc {uniprot_acc} from:\n"
                    f"{uniprot_table_dict[uniprot_acc]['seq']}\n"
                    "to:\n"
                    f"{uniprot_dict[ncbi_acc]['sequence']}"
                )
                record_seqs_to_update.add(
                    (
                        uniprot_acc,
                        uniprot_dict[ncbi_acc]['sequence'],
                        uniprot_dict[ncbi_acc]['sequence_date'],
                    )
                )
    
    if len(uniprot_insert_values) != 0:
        logger.warning(f"Inserting {len(uniprot_insert_values)} ew records into the Uniprots table")
        if args.sequence:
            columns = ['uniprot_accession', 'uniprot_name', 'sequence', 'seq_update_date']
        else:
            columns = ['uniprot_accession', 'uniprot_name']
    
        insert_data(connection, "Uniprots", columns, list(uniprot_insert_values))

        
    if len(record_names_to_update) != 0:
        logger.warning(f"Updating {len(record_names_to_update)} UniProt protein names in the local CAZyme db")
        with connection.begin():
            for record in tqdm(record_names_to_update, desc="Updating UniProt protein names"):
                connection.execute(
                    text(
                        "UPDATE Uniprots "
                        f"SET uniprot_name = {record[1]} "
                        f"WHERE uniprot_accession = '{record[0]}'"
                    )
                )

    if len(record_seqs_to_update) != 0:
        logger.warning(f"Updating {len(record_seqs_to_update)} UniProt protein sequences in the local CAZyme db")
        with connection.begin():
            for record in tqdm(record_seqs_to_update, desc="Updating UniProt protein seqs"):
                connection.execute(
                    text(
                        "UPDATE Uniprots "
                        f"SET sequence = {record[1]}, seq_update_date = {record[2]} "
                        f"WHERE uniprot_accession = '{record[0]}'"
                    )
                )

    return


def add_uniprot_genbank_relationships(uniprot_dict, connection):
    """Add Uniprots local db IDs to Genbanks table to link Genbanks and Uniprots records.

    :param uniprot_dict: dict of data from uniprot
        uniprot_dict[ncbi_acc] = {
            'uniprot_acc': uniprot_acc,
            'uniprot_entry_id': uniprot_entry_id,
            'protein_name': protein_name,
            'ec_numbers': ec_numbers,
            'sequence': sequence,
            'pdbs': all_pdbs,
        }
    :param connection: open sqlite3 engine connection

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    uniprot_table_dict = get_table_dicts.get_uniprot_table_dict(connection)
    # uniprot_table_dict = {acc: {db_id: int, name: str, seq: str, seq_date:str } }

    logger.warning(
        f"Linking {len(list(uniprot_dict.keys()))} Genbank records in the local CAZyme db"
        "to records in the Uniprots table"
    )
    with connection.begin():
        for ncbi_acc in tqdm(uniprot_dict, desc="Updating Genbanks records"):
            uniprot_acc = uniprot_dict[ncbi_acc]['uniprot_acc']
            try:
                uniprot_db_id = uniprot_table_dict[uniprot_acc]['db_id']
            except KeyError:
                logger.error(
                    f"Could not find a local db record for UniProt accession {uniprot_acc}\n"
                    "Not linking the accession to a record in the Genbanks table"
                )
                continue

            connection.execute(
                text(
                    "UPDATE Genbanks "
                    f"SET uniprot_id = {uniprot_db_id} "
                    f"WHERE genbank_accession = '{ncbi_acc}'"
                )
            )
    

def add_ec_numbers(uniprot_dict, connection, args):
    """Add EC numbers to the local CAZyme database

    :param uniprot_dict: dict containing data retrieved from UniProt
        uniprot_dict[ncbi_acc] = {
            'uniprot_acc': uniprot_acc, - str
            'protein_name': protein_name, -str
            'ec_numbers': ec_numbers, - set
            'sequence': sequence, - str
            'sequence_date': date seq was last updated yyyy-mm-dd
            'pdbs': all_pdbs, - set
        }
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    
    # load EC records in the local CAZyme db
    ec_table_dict = get_table_dicts.get_ec_table_dict(connection)

    # identify new EC numbers to add to the EC table
    existing_ecs = list(ec_table_dict.keys())
    ec_insert_values = set()

    for ncbi_acc in tqdm(uniprot_dict, desc="Identifying EC numbers to add to the local CAZyme db"):
        for ec_num in uniprot_dict[ncbi_acc]['ec_numbers']:
            if ec_num not in existing_ecs:
                ec_insert_values.add((ec_num,))

    if len(ec_insert_values) != 0:
        insert_data(connection, "Ecs", ["ec_number"], list(ec_insert_values))


def add_genbank_ec_relationships(uniprot_dict, gbk_dict, connection, args):
    """Add and update relationships between the Genbanks and Ecs table.

    :param uniprot_dict: dict containing data retrieved from UniProt
        uniprot_dict[ncbi_acc] = {
            'uniprot_acc': uniprot_acc, - str
            'protein_name': protein_name, -str
            'ec_numbers': ec_numbers, - set
            'sequence': sequence, - str
            'sequence_date': date seq was last updated yyyy-mm-dd
            'pdbs': all_pdbs, - set
        }
    :param gbk_dict: {gbk_acc: gbk_id}
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing
    """
    # load in EC records in the local CAZyme db
    # {ec_number: ec_id}
    ec_table_dict = get_table_dicts.get_ec_table_dict(connection)
    
    # load in gbk_ec table, contains the gbk-ec number relationships
    # {ec_id: {gbk ids}}
    ec_gbk_table_dict = get_table_dicts.get_ec_gbk_table_dict(connection)  

    # compile list of tuples to insert or delete into the Genbanks_Ecs table
    gbk_ec_insert_values = set()

    for ncbi_acc in tqdm(uniprot_dict, desc="Adding and updating EC-Genbanks relationships"):
        for ec_number in uniprot_dict[ncbi_acc]['ec_numbers']:
            gbk_db_id = gbk_dict[ncbi_acc]
            ec_id = ec_table_dict[ec_number]

            try:
                existing_relationships = ec_gbk_table_dict[ec_id]
                if gbk_db_id not in existing_relationships:
                    gbk_ec_insert_values.add( (gbk_db_id, ec_id) )

            except KeyError: # EC number is not linked to any gbk records so add new records
                gbk_ec_insert_values.add( (gbk_db_id, ec_id) )

    if len(gbk_ec_insert_values) != 0:
        insert_data(connection, "Genbanks_Ecs", ["genbank_id", "ec_id"], list(gbk_ec_insert_values))


def add_pdb_accessions(uniprot_dict, gbk_dict, connection, args):
    """Add PDB accessions to the local CAZyme database

    :param uniprot_dict: dict containing data retrieved from UniProt
        uniprot_dict[ncbi_acc] = {
            'uniprot_acc': uniprot_acc, - str
            'protein_name': protein_name, -str
            'ec_numbers': ec_numbers, - set
            'sequence': sequence, - str
            'sequence_date': date seq was last updated yyyy-mm-dd
            'pdbs': all_pdbs, - set
        }
    :param gbk_dict: dict representing data from the Genbanks table
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    pdbs_to_delete = set()  # only used if args.delete_old_pdbs is true
    relationships_to_delete = set()  # only used if args.delete_old_pdbs is true

    # load in PDB objects in the local CAZyme db
    pdb_table_dict = get_table_dicts.get_pdb_table_dict(connection)
    # {pdb_accession: pdb_db_id}

    # First, identify new PDB accessions to add to the database
    pdb_insert_values = set()
    for ncbi_acc in tqdm(uniprot_dict, desc="Identifying new PDBs to add to db"):
        for pdb in uniprot_dict[ncbi_acc]["pdbs"]:
            try:
                pdb_table_dict[pdb]
            except KeyError:
                pdb_insert_values.add( (pdb,) )
    
    if len(pdb_insert_values) != 0:
        logger.warning(f"Adding {len(pdb_insert_values)} PDB accessions to the database")
        insert_data(connection, "Pdbs", ["pdb_accession"], list(pdb_insert_values))


def add_pdb_gbk_relationships(uniprot_dict, gbk_dict, connection, args):
    """Add relationships between PDB accessions and proteins in the Genbanks table

    :param uniprot_dict: dict containing data retrieved from UniProt
        uniprot_dict[ncbi_acc] = {
            'uniprot_acc': uniprot_acc, - str
            'protein_name': protein_name, -str
            'ec_numbers': ec_numbers, - set
            'sequence': sequence, - str
            'sequence_date': date seq was last updated yyyy-mm-dd
            'pdbs': all_pdbs, - set
        }
    :param gbk_dict: dict representing data from the Genbanks table
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing.
    """
    # load the updated pdb table
    # {pdb_acc: pdb_db_id}
    pdb_table_dict = get_table_dicts.get_pdb_table_dict(connection)
    
    # load in Genbanks_Pdbs relationship table
    gbk_pdb_rel_table_dict = get_table_dicts.get_gbk_pdb_table_dict(connection)
    # {gbk_db_id: set(pdb_db_ids) }

    # convert the data from UniProt into a dict of {gbk_db_id: pdb_db_id} 
    # to identify pdb-protein relationships retrieved from UniProt
    gbk_pdb_insert_values = set()

    for ncbi_acc in tqdm(uniprot_dict, desc="Identifying new protein-PDB relationships to add to db"):

        if len(uniprot_dict[ncbi_acc]['pdbs']) == 0:  # not relationships to add
            continue

        uniprot_acc = uniprot_dict[ncbi_acc]['uniprot_acc']
        try:
            gbk_db_id = gbk_dict[ncbi_acc]
        except KeyError:
            logger.error(
                f"Mapped the GenBank accession '{ncbi_acc}' to the UniProt accession\n"
                f"'{uniprot_acc}' but the GenBank accession is not in the local CAZyme database\n"
                f"therefore, not adding protein data for GBK:{ncbi_acc}/UniProt:{uniprot_acc}"
                "to the local CAZyme database."
            )
            continue

        # check if there any existing relationships for the gbk record
        try:
            existing_pdb_relationships = gbk_pdb_rel_table_dict[gbk_db_id]  
        except KeyError:
            existing_pdb_relationships = set()

        for pdb_acc in uniprot_dict[ncbi_acc]["pdbs"]:
            try:
                pdb_db_id = pdb_table_dict[pdb_acc]
            except KeyError:
                logger.error(
                    f"Retrieved PDB:{pdb_acc} from UniProt.\n"
                    "Cannot link to Genbanks table records because PDB not listed in the Pdbs table"
                )
                continue
            
            if pdb_db_id not in gbk_pdb_rel_table_dict:
                gbk_pdb_insert_values.add( (gbk_db_id, pdb_db_id) )

    if len(gbk_pdb_insert_values) != 0:
        insert_data(connection, "Genbanks_Pdbs", ["genbank_id", "pdb_id"], list(gbk_pdb_insert_values))


def add_uniprot_taxs(uniprot_dict, connection, args):
    """Add taxonomic classifications (genus, species) to the local CAZyme database UniprotTaxs table.

    :param uniprot_dict: dict containing data retrieved from UniProt
        uniprot_dict[ncbi_acc] = {
            'uniprot_acc': uniprot_acc, - str
            'protein_name': protein_name, -str
            'ec_numbers': ec_numbers, - set
            'sequence': sequence, - str
            'sequence_date': date seq was last updated yyyy-mm-dd
            'pdbs': all_pdbs, - set
        }
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing.
    """
    # add tax data to the local db
    # load the table
    ut_table_dict, ut_tax_dict = get_table_dicts.get_uniprottax_table_dict(connection)
    # {db id: {'genus': str, 'species': str}}
    # {'genus species': db id}

    taxs_to_add = set()

    for ncbi_acc in tqdm(uniprot_dict, desc="Identifying taxs to add to db"):
        genus = uniprot_dict[ncbi_acc]['genus']
        species = uniprot_dict[ncbi_acc]['species']

        try:
            ut_tax_dict[f"{genus} {species}"]
        except KeyError:
            taxs_to_add.add( (genus, species) )

    if len(taxs_to_add) > 0:
        insert_data(connection, "UniprotTaxs", ["genus", "species"], list(taxs_to_add))

    add_uniprot_tax_relationships(uniprot_dict, connection, args)

    
def add_uniprot_tax_relationships(uniprot_dict, connection, args):
    """Link Uniprot records to UniprotTaxs table.

    :param uniprot_dict: dict containing data retrieved from UniProt
        uniprot_dict[ncbi_acc] = {
            'uniprot_acc': uniprot_acc, - str
            'protein_name': protein_name, -str
            'ec_numbers': ec_numbers, - set
            'sequence': sequence, - str
            'sequence_date': date seq was last updated yyyy-mm-dd
            'pdbs': all_pdbs, - set
        }
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing.
    """
    ut_table_dict, ut_tax_dict = get_table_dicts.get_uniprottax_table_dict(connection)
    # {db id: {'genus': str, 'species': str}}
    # {'genus species': db id}
    
    uniprot_table_dict = get_table_dicts.get_uniprot_table_dict(connection)
    # {acc: {name: str, gbk_id: int, seq: str, seq_date:str } }

    print('ut_tax_dict:', ut_tax_dict)

    relationships = {}  # {uniprot db id: uniprot tax db id}

    for ncbi_acc in tqdm(uniprot_dict, desc="Identifying Uniprot-Tax relationships"):
        genus = uniprot_dict[ncbi_acc]['genus']
        species = uniprot_dict[ncbi_acc]['species']
        uniprot_acc = uniprot_dict[ncbi_acc]['uniprot_acc']

        uni_db_id = uniprot_table_dict[uniprot_acc]['db_id']
        tax_db_id = ut_tax_dict[f"{genus} {species}"]

        relationships[uni_db_id] = tax_db_id

    with connection.begin():
        for uni_db_id in tqdm(relationships, desc="Adding Uniprot-Tax relationships"):
            connection.execute(
                text(
                    "UPDATE Uniprots "
                    f"SET uniprot_tax_id = {relationships[uni_db_id]} "
                    f"WHERE uniprot_id = '{uni_db_id}'"
                )
            )
