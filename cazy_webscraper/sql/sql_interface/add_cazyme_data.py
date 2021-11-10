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

from tqdm import tqdm

from cazy_webscraper.sql.sql_interface import insert_data, get_gbk_table_dict, get_table_dicts
from cazy_webscraper.sql.sql_orm import Session, CazyFamily, Kingdom, Taxonomy


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
    
    :param taxa_dict: dict of taxa data {kingdom: set(organism)}
    :param connection: open sqlalchemy connection to SQLite engine
    
    Return nothing
    """
    # retrieve db kingdom objects for retrieving the kingdom_id for the Taxs table
    kingdom_table_dict = get_table_dicts.get_kingdom_table_dict(connection)  # {kingdom: kingdom_id}

    # compare taxa already in the db against taxa retrieved from the CAZy txt file
    # to identify new taxa objects to be added to the db

    taxonomy_db_insert_values = []

    for kingdom in tqdm(
        taxa_dict,
        total=len(list(taxa_dict.keys())),
        desc='Create tax objects per Kingdom',
    ):
        kingdom_id = kingdom_table_dict[kingdom]
        
        existing_taxa_records = taxa_dict[kingdom]
        
        for organism in taxa_dict[kingdom]:  # organisms from the CAZy txt file
            if organism not in existing_taxa_records:  # new record to add
                genus = organism.split(" ")[0]
                species = ' '.join(organism.split(" ")[1:])
                taxonomy_db_insert_values.append( (genus, species, kingdom_id,) )

    if len(taxonomy_db_insert_values) != 0:
        insert_data(connection, 'Taxs', ['genus', 'species', 'kingdom_id'], taxonomy_db_insert_values)


def add_cazy_families(cazy_data, connection):
    """Add CAZy families and subfamilies to local CAZyme database
    
    :param cazy_data: dict of data retrieved from CAZy
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


def load_taxa_fam_data(connection):
    """Load taxonomy and CAZy fam tables into memory.
    
    :param connection: open sqlalchemy connection to an SQLite db engine
    
    Return dict representing the Taxs table and a dict of the Families table
    """
    with Session(bind=connection) as session:
        all_families = session.query(CazyFamily).all()

    db_fam_dict = {}

    for fam in all_families:
        if fam.subfamily is None:
            subfam = '_'
        else:
            subfam = fam.subfamily
            
        db_fam_dict[f"{fam.family} {subfam}"] = fam.family_id

    # create dict of tax data from db, valued by the tax obj id number
    with Session(bind=connection) as session:
        all_taxa = session.query(Taxonomy).all()
        
    db_tax_dict = {}
    for taxa in all_taxa:
        if len(taxa.species) == 0:
            db_tax_dict[f"{taxa.genus}"] = taxa.taxonomy_id
        else:
            db_tax_dict[f"{taxa.genus} {taxa.species}"] = taxa.taxonomy_id
    
    return db_tax_dict, db_fam_dict


def add_genbanks(cazy_data, db_tax_dict, db_fam_dict, connection):
    """Add GenBank accessions with tax data to the db
    
    :param cazy_data: dict of CAZy data
    :param connection: open sqlalchemy connection to an SQLite db
    
    Return set of tuples, of (genbank_acc, fam_id)
    """
    # create set of tuples for inserting into the Genbank db table
    # associates a genbank accession with its source organisms db tax id num
    gbk_db_insert_values = set()

    # ...and create set of tuples for inserting into the CazyFamilies_Genbank table
    # associates genbank acc db id nums with CAZy fam db ids
    gbk_fam_values = set()  # { (gbk_accession, fam_db_id) } 

    for genbank_accession in tqdm(cazy_data, desc='Adding GenBank to db'):
        organism = list(cazy_data[genbank_accession]['organism'])[0]
        tax_id = db_tax_dict[organism]
        
        gbk_db_insert_values.add( (genbank_accession, tax_id) )
        
        for cazy_fam in  cazy_data[genbank_accession]['families']:
            subfamilies = cazy_data[genbank_accession]['families'][cazy_fam]
            
            for cazy_subfam in subfamilies:
                if cazy_subfam is None:
                    fam_name = f"{cazy_fam} _"
                
                else:
                    fam_name = f"{cazy_fam} {cazy_subfam}"
            
            fam_id = db_fam_dict[fam_name]
            gbk_fam_values.add( (genbank_accession, fam_id) )

    insert_data(connection, 'Genbanks', ['genbank_accession', 'taxonomy_id'], list(gbk_db_insert_values))

    return gbk_fam_values


def add_genbank_fam_relationships(gbk_fam_values, connection):
    """Add GenBank accession - CAZy family relationships to db
    
    :param gbk_fam_values: set of tuples (gbk_acc, fam_id)
    :param connection: open sqlalchemy connection to an SQLite db engine

    Return nothing
    """
    # retrieve genbank id numbers from the local database
    db_gbk_dict = get_gbk_table_dict(connection)

    gbk_fam_db_insert_values = set()  # { (gbk_id(int), fam_id(int), ), }

    for gbk_tuple in tqdm(gbk_fam_values, desc='Adding Protein-Fam relationships'):
        gbk_id = db_gbk_dict[gbk_tuple[0]]
        fam_id = gbk_tuple[1]
        gbk_fam_db_insert_values.add( (gbk_id, fam_id) )

    insert_data(connection, 'Genbanks_CazyFamilies', ['genbank_id', 'family_id'], list(gbk_fam_db_insert_values))
