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
"""Test script cazy_webscraper.sql.sql_interface.add_data.add_gtdb_tax.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace

from cazy_webscraper.sql.sql_interface.add_data import add_gtdb_tax


from sqlalchemy import (
    MetaData,
    create_engine,
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

# Use the declarative system
# Database structured in NF1
metadata_obj = MetaData()
Base = declarative_base()
Session = sessionmaker()


def test_add_gtdb_tax(monkeypatch):
    """Test add_gtdb_taxs()"""
    def mock_exisiting_taxs(*args, **kwards):
        return {}

    def mock_return_none(*args, **kwards):
        return

    monkeypatch.setattr(add_gtdb_tax, "get_gtdb_table_dict", mock_exisiting_taxs)
    monkeypatch.setattr(add_gtdb_tax, "insert_data", mock_return_none)

    gtdb_lineages = {
        'GCA_0000000000.1': {'lineage': 'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Pectobacteriaceae; g__Pectobacterium; s__Pectobacterium brasiliense',
        'release': '20221101'},
        'GCA_0000000001.1': {'lineage': 'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Pectobacterium; s__Pectobacterium aquaticum',
         'release': '20221101'},
    }

    connection = None

    add_gtdb_tax.add_gtdb_taxs(gtdb_lineages, connection)


def test_relate_gtdb(monkeypatch, db_path):
    """Test add_genome_gtdb_relations()"""
    def mock_genome_table(*args, **kwards):
        return {
            'GCA_0000000000.1': {'db_id': 1, 'gtdb_id': 1},
            'GCA_0000000001.1': {'db_id': 2, 'gtdb_id': 2},
            'GCA_0000000002.1': {'db_id': 3, 'gtdb_id': 3},
        }

    def mock_gtdb_table_dict(*args, **kwards):
        lin_1 = 'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Pectobacteriaceae; g__Pectobacterium; s__Pectobacterium brasiliense'
        lin_1 = tuple(lin_1.split(";"))
        lin_2 = 'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Pectobacterium; s__Pectobacterium aquaticum'
        lin_2 = tuple(lin_2.split(";"))
        return {
            1: lin_1,
            2: lin_2,
        }

    def mock_return_none(*args, **kwards):
        return

    monkeypatch.setattr(add_gtdb_tax, "get_genome_table", mock_genome_table)
    monkeypatch.setattr(add_gtdb_tax, "get_gtdb_table_dict", mock_gtdb_table_dict)

    gtdb_lineages = {
        'GCA_0000000000.1': {'lineage': 'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Pectobacteriaceae; g__Pectobacterium; s__Pectobacterium brasiliense',
        'release': '20221101'},
        'GCA_0000000001.1': {'lineage': 'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Pectobacterium; s__Pectobacterium aquaticum',
         'release': '20221101'},
        'GCA_0000000006.1': {'lineage': 'd__Bacteriaaaaa; p__Proteobacteriaaaaa; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Pectobacterium; s__Pectobacterium aquaticum',
         'release': '20221101'},
    }

    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
            update_genome_lineage=True,
        )
    }

    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=False, future=True)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)  # allows for calls to session later on when required

    with engine.connect() as conn:
        savepoint = conn.begin_nested()
        add_gtdb_tax.add_genome_gtdb_relations(gtdb_lineages, argsdict['args'], conn)
        savepoint.rollback()
