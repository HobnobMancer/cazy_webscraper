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
"""Test script cazy_webscraper.sql.sql_interface.add_data.add_ncbi_tax_data.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace

from cazy_webscraper.sql.sql_interface.add_data import add_ncbi_tax_data

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


def test_add_ncbi_taxs(monkeypatch):
    """Test add_ncbi_taxonomies()"""
    tax_dict = {
        1 : {
            'kingdom': 'kingdom',
            'phylum': 'phylum',
            'class': 'class',
            'order': 'order',
            'family': 'family',
            'genus': 'genus',
            'species': 'species',
            'strain': 'strain',
        },
        2 : {
            'kingdom': 'kingdom',
            'phylum': 'phylum',
            'class': 'class',
            'order': 'order',
            'family': 'family',
            'genus': 'genus',
            'species': 'species',
            'strain': 'strain',
        }
    }
    connection = None
    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
            seq_update=True,
            update_taxs=True,
        )
    }

    def mock_get_tax_table(*args, **kwards):
        return {
            1 : {
                'kingdom': 'kingdom',
                'phylum': 'phylum',
                'class': 'class',
                'order': 'order',
                'family': 'family',
                'genus': 'genus',
                'species': 'species',
                'strain': 'strain',
            },
        }

    def mock_return_none(*args, **kwards):
        return

    monkeypatch.setattr(add_ncbi_tax_data, "get_ncbi_tax_table", mock_get_tax_table)
    monkeypatch.setattr(add_ncbi_tax_data, "update_ncbi_tax_records", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "insert_data", mock_return_none)

    add_ncbi_tax_data.add_ncbi_taxonomies(
        tax_dict,
        connection,
        argsdict['args'],
    )


def test_update_taxs(db_path):
    """Test update_ncbi_tax_records()"""
    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=False, future=True)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)  # allows for calls to session later on when required
    connection = engine.connect()

    add_ncbi_tax_data.update_ncbi_tax_records(
        [('1', 'kingdom1', 'phylum1', 'tax1', 'order1', 'family1', 'genus1', 'species1', 'strain1')],
        connection,
        unit_test=True,
    )


def test_update_gbk_tax(db_path, monkeypatch):
    """Test update_genbank_ncbi_tax()"""
    def mock_ncbi_tax_table(*args, **kwards):
        return {
            1: 1,
        }

    def mock_tax_gbk_table(*args, **kwards):
        return {
            1: 2,
            2: 2,
        }

    def mock_get_tax_table(*args, **kwards):
        return {
            1 : {
                'kingdom': 'kingdom',
                'phylum': 'phylum',
                'class': 'class',
                'order': 'order',
                'family': 'family',
                'genus': 'genus',
                'species': 'species',
                'strain': 'strain',
            },
        }

    monkeypatch.setattr(add_ncbi_tax_data, "get_ncbi_tax_table", mock_get_tax_table)
    monkeypatch.setattr(add_ncbi_tax_data, "get_ncbi_tax_table", mock_ncbi_tax_table)
    monkeypatch.setattr(add_ncbi_tax_data, "get_no_tax_gbk_table_dict", mock_tax_gbk_table)

    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=False, future=True)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)  # allows for calls to session later on when required
    connection = engine.connect()

    tax_prot_dict = {
        1: {'proteins': [1]},
    }

    argsdict = {
        "args": Namespace(
            update_gbk=True,
        )
    }

    add_ncbi_tax_data.update_genbank_ncbi_tax(
        tax_prot_dict,
        connection,
        argsdict['args'],
        unit_test=True,
    )
