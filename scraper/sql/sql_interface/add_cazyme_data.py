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

from scraper.sql.sql_interface import insert_data
from scraper.sql.sql_orm import Session, Kingdom


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


def add_source_organisms(taxa_data, connection):
    """Add taxonomy (source organism) data to the local CAZyme database
    
    :param taxa_data: dict of taxa data {kingdom: set(organism)}
    :param connection: open sqlalchemy connection to SQLite engine
    
    Return nothing
    """
    taxonomy_db_insert_values = []
    with Session(bind=connection) as session:
        for kingdom in tqdm(taxa_data, desc='Adding Tax data to db per Kingdom'):
            print(kingdom)
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


