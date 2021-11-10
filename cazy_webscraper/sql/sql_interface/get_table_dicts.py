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
"""Retrieve all objects from a db table and parse the data to build a dict, repr the current table state."""


from cazy_webscraper.sql.sql_orm import CazyFamily, Genbank, Kingdom, Session


def get_fams_table_dict(connection):
    """Create dict of objects present in the CazyFamilies table.
    
    :param connection: open sqlalchemy db engine connection
    
    Return dict {family subfamily: db_family_id}
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
    
    return db_fam_dict


def get_gbk_table_dict(connection):
    """Compile a dict of the data in the Genbanks table
    
    :param connection: open connection to an SQLite3 database
    
    Return dict {genbank_accession: 'taxa_id': int, 'gbk_id': int}
    """
    with Session(bind=connection) as session:
        all_genbank = session.query(Genbank).all()

    db_gbk_dict = {}  # {genbank_accession: 'taxa_id': str, 'id': int}
    
    for gbk in all_genbank:
        db_gbk_dict[f"{gbk.genbank_accession}"] = {
            'taxa_id': gbk.taxonomy_id,
            'id': gbk.genbank_id
        }
    
    return db_gbk_dict


def get_kingdom_table_dict(connection):
    """Load and parse the Kingdoms table from the db and compile a dict {kgnd: id}
    
    :param connection:
    
    Return dict {kingdom: kindom_db_id}
    """
    with Session(bind=connection) as session:
        kingdom_table = session.query(Kingdom).all()
    
    kingdom_dict = {}  # {kingdom: kindom_db_id}
    
    for kingdom_obj in kingdom_table:
        kingdom_dict[kingdom_obj.kingdom] = kingdom_obj.kingdom_id
        
    return kingdom_dict


