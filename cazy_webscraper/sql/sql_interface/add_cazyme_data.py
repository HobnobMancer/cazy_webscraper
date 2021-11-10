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

from cazy_webscraper.sql.sql_interface import insert_data, get_gbk_table_dict
from cazy_webscraper.sql.sql_orm import Session, CazyFamily, Kingdom, Taxonomy


def add_kingdoms(cazy_data, connection):
    """Create a dict of taxa information and add kingdoms to db.
    
    :param cazy_data: dict of CAZy data
    :param connection: open sqlalchemy conenction to an SQLite db engine

    Return taxa_dict, dict of taxonomy data {kingdom: set(organism)}
    """
    taxa_dict = {}  # {kingdom: {organism,}}
    
    for genbank_accession in tqdm(cazy_data, "Compiling taxa data"):
        kingdom = list(cazy_data[genbank_accession]['kingdom'])[0]
        organism = list(cazy_data[genbank_accession]['organism'])[0]
        
        try:
            taxa_dict[kingdom].add(organism)
        except KeyError:
            taxa_dict[kingdom] = {organism}

    kingdoms_db_insert_values = [(kingdom,) for kingdom in taxa_dict.keys()]

    insert_data(connection, 'Kingdoms', ['kingdom'], kingdoms_db_insert_values)

    return taxa_dict


def add_source_organisms(taxa_data, connection):
    """Add taxonomy (source organism) data to the local CAZyme database
    
    :param taxa_data: dict of taxa data {kingdom: set(organism)}
    :param connection: open sqlalchemy connection to SQLite engine
    
    Return nothing
    """
    taxonomy_db_insert_values = []
    with Session(bind=connection) as session:
        for kingdom in tqdm(
            taxa_data,
            total=len(list(taxa_data.keys())),
            desc='Adding Tax data to db per Kingdom',
        ):
            # query db to get the kingdom db object, retrieve only the kingdom ID number
            found_kingdom = session.query(Kingdom.kingdom_id).\
                filter(Kingdom.kingdom==kingdom).\
                first()[0]
            
            for organism in taxa_data[kingdom]:
                genus = organism.split(" ")[0]
                species = ' '.join(organism.split(" ")[1:])
                taxonomy_db_insert_values.append((genus, species, found_kingdom))

    insert_data(connection, 'Taxs', ['genus', 'species', 'kingdom_id'], taxonomy_db_insert_values)


def add_cazy_families(cazy_data, connection):
    """Add CAZy families and subfamilies to local CAZyme database
    
    :param cazy_data: dict of data retrieved from CAZy
    :param connection: open sqlalchemy connection to an SQLite db engine
    
    Return nothing"""
    families_db_insert_values = set()

    for genbank_accession in tqdm(cazy_data, desc='Adding CAZy fams to db'):
        for cazy_fam in cazy_data[genbank_accession]["families"]:
            subfamilies = cazy_data[genbank_accession]["families"][cazy_fam]
            if None not in subfamilies:
                cazy_data[genbank_accession]["families"][cazy_fam].add(None)
            
            for cazy_subfam in subfamilies:
                families_db_insert_values.add( (cazy_fam, cazy_subfam) )

    insert_data(connection, 'CazyFamilies', ['family', 'subfamily'], list(families_db_insert_values))


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
