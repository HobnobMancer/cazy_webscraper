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
"""Add CAZyme data to a local SQLite database"""


import logging

from sqlalchemy import delete, text
from tqdm import tqdm

from cazy_webscraper.sql.sql_interface import insert_data, get_table_dicts
from cazy_webscraper.sql.sql_orm import genbanks_families


def add_kingdoms(taxa_dict, connection):
    """Add new Kingdoms objects to database.
    
    Check existing kingdom objects in the db against kingdoms retrieved from the 
    CAZy txt file, so as to only add new kingdoms.

    :param taxa_dict: dict of kingdoms and organisms from the cazy_data dict
        {kingdom: {organism}}
    :param connection: open sqlalchemy connection to a local SQLite db engine
    
    Return nothing
    """
    kingdom_table_dict = get_table_dicts.get_kingdom_table_dict(connection)
    # dict {kingdom: {organisms}}

    # retrieve the Kingdoms retrieved from the CAZy txt file
    cazy_kingdoms = list(taxa_dict.keys())
    existing_kingdom_records = list(kingdom_table_dict.keys())

    # create list of tuples for db insert
    kingdoms_db_insert_values = [
        (kngdm,) for kngdm in cazy_kingdoms if kngdm not in existing_kingdom_records
    ]

    if len(kingdoms_db_insert_values) != 0:
        insert_data(connection, 'Kingdoms', ['kingdom'], kingdoms_db_insert_values)
    
    return


def add_source_organisms(taxa_dict, connection):
    """Add taxonomy (source organism) data to the local CAZyme database
    
    :param taxa_dict: dict of taxa data {kingdom: set(organism)} to be added to the db
    :param connection: open sqlalchemy connection to SQLite engine
    
    Return nothing
    """
    logger = logging.getLogger(__name__)

    # retrieve db kingdom objects for retrieving the kingdom_id for the Taxs table
    kingdom_table_dict = get_table_dicts.get_kingdom_table_dict(connection)  # {kingdom: kingdom_id}
    tax_table_dict = get_table_dicts.get_taxs_table_dict(connection)
    # {genus species: {'tax_id': db_tax_id, 'kingdom_id': kingdom_id}

    # compare taxa already in the db against taxa retrieved from the CAZy txt file
    # to identify new taxa objects to be added to the db

    taxonomy_db_insert_values = []
    records_to_update = set()  # used incase kingdom has changed for a species

    for kingdom in tqdm(
        taxa_dict,
        total=len(list(taxa_dict.keys())),
        desc='Create tax objects per Kingdom',
    ):
        kingdom_id = kingdom_table_dict[kingdom]
        
        existing_taxa_records = list(tax_table_dict.keys())  # list of scientific names

        organisms = taxa_dict[kingdom]
        
        for organism in organisms:  # organisms from the CAZy txt file
            
            if organism not in existing_taxa_records:  # new record to add
                genus = organism.split(" ")[0]
                species = ' '.join(organism.split(" ")[1:])
                taxonomy_db_insert_values.append( (genus, species, kingdom_id,) )
            
            else:  # check kingdom is correct
                existing_record_kngdm_id = tax_table_dict[organism]['kingdom_id']
                if existing_record_kngdm_id != kingdom_id:
                    records_to_update.add( (genus, species, kingdom_id,) )
                
    if len(taxonomy_db_insert_values) != 0:
        logger.info(
            f"Adding {len(taxonomy_db_insert_values)} new tax records to the db"
        )
        insert_data(
            connection,
            'Taxs',
            ['genus', 'species', 'kingdom_id'],
            taxonomy_db_insert_values,
        )
    else:
        logger.info(
            "No new tax records to add to the db"
        )
    if len(records_to_update) != 0:
        logger.info(
            f"Updating the parent Kingdom for {len(records_to_update)} tax records in the db"
        )
        for record in records_to_update:
            connection.execute(
                text(
                    "UPDATE Taxs "
                    f"SET kingdom_id = {record[2]} "
                    f"WHERE genus = '{record[0]}' AND species = '{record[1]}'"
                )
            )

    return


def add_cazy_families(cazy_data, connection):
    """Add CAZy families and subfamilies to local CAZyme database
    
    :param cazy_data: dict of data extracted from the txt file
        {gbk_accession: {kingdom:str, organism:str, families{fam:subfam}}}
    :param connection: open sqlalchemy connection to an SQLite db engine
    
    Return nothing"""
    logger = logging.getLogger(__name__)

    # get list of CAZy (sub)families already present in the db
    fam_table_dict = get_table_dicts.get_fams_table_dict(connection)  # {family subfamily: db_family_id}

    existing_fam_records = list(fam_table_dict.keys()) 

    families_db_insert_values = set()  # new fam records to add to db

    for genbank_accession in tqdm(cazy_data, desc='Extracting CAZy fams from CAZy data'):
        for cazy_fam in cazy_data[genbank_accession]["families"]:
            subfamilies = cazy_data[genbank_accession]["families"][cazy_fam]
            
            for subfam in subfamilies:
                if subfam is None:
                    cazy_subfam = "_"
                else:
                    cazy_subfam = subfam
                fam_key = f"{cazy_fam} {cazy_subfam}"
                
                if fam_key in existing_fam_records:
                    families_db_insert_values.add( (cazy_fam, cazy_subfam) )
            
    if len(families_db_insert_values) != 0:
        logger.info(
            f"Inserting {len(families_db_insert_values)} "
            "new family records into the CazyFamilies table"
        )
        insert_data(connection, 'CazyFamilies', ['family', 'subfamily'], list(families_db_insert_values))
    else:
        logger.info("Found no new CAZy family records to add the CazyFamilies table")

    return


def add_genbanks(cazy_data, connection):
    """Add GenBank accessions with tax data to the db
    
    :param cazy_data: dict of CAZy data
    :param connection: open sqlalchemy connection to an SQLite db
    
    Return Nothing
    """
    logger = logging.getLogger(__name__)

    # retrieve existing records from the db
    gbk_table_dict = get_table_dicts.get_gbk_table_dict(connection)
    existing_gbk_records = list(gbk_table_dict.keys())

    taxa_table_dict = get_table_dicts.get_taxs_table_dict(connection)
    # {genus species: {'tax_id': db_tax_id, 'kingdom_id': kingdom_id}

    gbk_record_updates = set()  # {gbk_accession: 'taxa_id': (new taxa_id) int, 'gbk_id': int}
    gbk_db_insert_values = set()

    for gbk_accession in tqdm(cazy_data, desc="Creating Genbank records to insert"):
        if gbk_accession not in existing_gbk_records:
            organism = cazy_data[gbk_accession]['organism']
            taxa_id = taxa_table_dict[organism]['tax_id']
            gbk_db_insert_values.add( (gbk_accession, taxa_id,) )
        
        else:  # check if need to update taxa_id
            # get the taxa_id for the existing record
            existing_record_id = gbk_table_dict[gbk_accession]['taxa_id']
            
            # get the taxa_id for the organism listed in the CAZy txt file
            organism = cazy_data[gbk_accession]['organism']
            cazy_data_taxa_id = taxa_table_dict[organism]
            
            if cazy_data_taxa_id != existing_record_id:
                # need to update the record
                gbk_record_updates.add( (cazy_data_taxa_id, gbk_table_dict[gbk_accession]['gbk_id']) )

    if len(gbk_db_insert_values) != 0:
        logger.info(
            f"Inserting {len(gbk_db_insert_values)} GenBank accessions into the db"
        )
        insert_data(connection, 'Genbanks', ['genbank_accession', 'taxonomy_id'], list(gbk_db_insert_values))
    
    if len(gbk_record_updates) != 0:
        logger.info(
            f"Updating {len(gbk_record_updates)} Genbank table records with new taxonomy IDs"
        )
        for record in gbk_record_updates:
            connection.execute(
                text(
                    "UPDATE Taxs "
                    f"SET taxonomy_id = {record[1]} "
                    f"WHERE genbank_id = '{record[0]}'"
                )
            )
    
    return


def add_genbank_fam_relationships(cazy_data, connection, args):
    """Add GenBank accession - CAZy family relationships to db
    
    :param cazy_data: dict of data extracted from the txt file
        {gbk_accession: {kingdom:str, organism:str, families{fam:subfam}}}
    :param connection: open sqlalchemy connection to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing
    """
    logger = logging.getLogger(__name__)

    gbk_fam_db_insert_values = set()  # new records to add
    gbk_fam_records_to_del = set()  # records/relationships to delete

    # get dict of GenBank and CazyFamilies tables, used for getting gbk_ids  and fam_ids of accessions and
    # families without entries in the CazyFamilies_Genbanks table
    gbk_table_dict = get_table_dicts.get_gbk_table_dict(connection) 
    # {genbank_accession: 'taxa_id': int, 'gbk_id': int}
    fam_table_dict = get_table_dicts.get_fams_table_dict(connection) 
    # {fam: fam_id}

    # load current relationships in the db
    gbk_fam_table_dict = get_table_dicts.get_gbk_fam_table_dict(connection)
    # {genbank_accession: families: {fam: fam_id}, id: int (db_id)}
    
    gbks_with_existing_relationships = list(gbk_fam_table_dict.keys())

    for genbank_accession in tqdm(cazy_data, desc="Extract Genbank-Family relationships from CAZy data"):
        
        if genbank_accession in gbks_with_existing_relationships:
            # identify which genbank-fam relationships in the CAZy data to delete or add
            existing_families_rel_dict = gbk_fam_table_dict[genbank_accession]['families']  # {fam: fam_id}
            
            families = cazy_data[genbank_accession]['families']
            
            cazy_families = set()  # used for checking which existing relationships are no longer in CAZy
            
            for fam in families:
                if families[fam] is None:
                    family = f"{fam} _"
                else:
                    family = f"{fam} {families[fam]}"
                cazy_families.add(family)
                
                if family not in (list(existing_families_rel_dict.keys())):
                    # add new relationship
                    genbank_id = gbk_table_dict[genbank_accession]['id']
                    fam_id = fam_table_dict[family]
                    gbk_fam_db_insert_values.add( (genbank_id, fam_id,) )
                    
                # else: gbk-fam already in db
                
            # check for relationships to delete: a relationship that is in the db but no longer in CAZy
            for family in list(existing_families_rel_dict.keys()):
                if family not in cazy_families:  # relationship no longer in CAZy
                    genbank_id = gbk_table_dict[genbank_accession]['id']
                    fam_id = fam_table_dict[family]
                    gbk_fam_records_to_del.add( (genbank_id, fam_id) )
                
                # else: fam relationship in the db is still in CAZy
        
        else:
            # add all CAZy (sub)family associations to db for this GenBank accession
            families = cazy_data[genbank_accession]['families']
            for fam in families:
                if families[fam] is None:
                    family = f"{fam} _"
                else:
                    family = f"{fam} {families[fam]}"
                
                fam_id = fam_table_dict[family]
                
                genbank_id = gbk_table_dict[genbank_accession]['id']
                
                gbk_fam_db_insert_values.add( (genbank_id, fam_id,) )
        
        
    if len(gbk_fam_db_insert_values) != 0:
        logger.info(
            f"Adding {len(gbk_fam_db_insert_values)} new GenBank accession - "
            "CAZy (sub)family relationships to the db"
        )
        insert_data(
            connection,
            'Genbanks_CazyFamilies',
            ['genbank_id', 'family_id'],
            list(gbk_fam_db_insert_values),
        )
    else:
        logger.info(
            "No new Genbank accession-CAZy (sub)family relationships to add to the db"
        )
    
    if (len(gbk_fam_records_to_del) != 0) and args.delete_old_relationships:
        logger.info(
            "Deleting {(len(gbk_fam_records_to_del)} GenBank accession - "
            "CAZy (sub)family relationships\n"
            "that are the db but are no longer in CAZy"
        )
        for record in gbk_fam_records_to_del:
            # record = (genbank_id, fam_id,)
            stmt = (
                delete(genbanks_families).\
                where(genbanks_families.c.genbank_id == record[0]).\
                where(genbanks_families.c.family_id == record[1])
            )
            connection.execute(stmt)

    return
