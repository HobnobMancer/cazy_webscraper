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


from multiprocessing import connection
from cattr import gen
from matplotlib.pyplot import connect
import pytest

from argparse import Namespace

from cazy_webscraper.sql.sql_interface.add_data import add_genome_data


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


@pytest.fixture
def ncbi_genome_dict():
    ncbi_genome_dict = {
        'Assemb1': {
            'gbk_acc': 'GCA_000000000.1',
            'gbk_uid': 'GCA_000000000.1',
            'refseq_acc': 'GCA_000000000.1',
            'refseq_uid': 'GCA_000000000.1',
        },
        'Assemb2': {
            'gbk_acc': 'GCA_000000000.1',
            'gbk_uid': 'GCA_000000000.1',
            'refseq_acc': 'GCA_000000000.1',
            'refseq_uid': 'GCA_000000000.1',
        },
    }
    return ncbi_genome_dict


@pytest.fixture
def genomes_of_interest():
    genomes_of_interest = ['Assemb1']
    return genomes_of_interest


def test_add_assembly_data(monkeypatch):
    """Test add_assembly_data()"""

    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
            update=True,
        )
    }

    def mock_genome_table(*args, **kwards):
        return {'Assemb1': 1}

    def mock_return_none(*args, **kwards):
        return

    monkeypatch.setattr(add_genome_data, "get_assembly_table", mock_genome_table)
    monkeypatch.setattr(add_genome_data, "update_genomic_data", mock_return_none)
    monkeypatch.setattr(add_genome_data, "add_genomic_data", mock_return_none)
    monkeypatch.setattr(add_genome_data, "add_prot_genome_relationships", mock_return_none)

    ncbi_genome_dict = {
        'Assemb1': '',
        'Assemb2': '',
    }
    connection = None
    assembly_prot_dict = {}
    gbk_dict = {}

    add_genome_data.add_assembly_data(
        assembly_prot_dict,
        ncbi_genome_dict,
        gbk_dict,
        connection,
        argsdict['args'],
    )


def test_add_genomic_data(monkeypatch, ncbi_genome_dict, genomes_of_interest):
    """Test add_genomic_data()"""
    connection = None

    def mock_return_none(*args, **kwards):
        return

    monkeypatch.setattr(add_genome_data, "insert_data", mock_return_none)

    add_genome_data.add_genomic_data(
        ncbi_genome_dict,
        genomes_of_interest,
        connection,
    )


def test_update_genomic_data(genomes_of_interest, db_path):
    """Test update_genomic_data()"""
    genome_table_dict = {'ASM75530v1': 1}

    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=False, future=True)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)  # allows for calls to session later on when required

    connection = engine.connect()

    ncbi_genome_dict = {'ASM75530v1': {
        'gbk_acc': 'GCA_000755305.1',
        'gbk_uid': '1203388',
        'refseq_acc': 'GCF_000755305.1',
        'refseq_uid': '1523468',
    }}

    add_genome_data.update_genomic_data(
        genomes_of_interest,
        genome_table_dict,
        ncbi_genome_dict,
        connection,
    )


def test_add_prot_g_relations(monkeypatch):
    """Test add_prot_genome_relationships()"""
    def mock_gbk_genome_table(*args, **kwards):
        return []
    
    def mock_return_none(*args, **kwards):
        return

    monkeypatch.setattr(add_genome_data, "get_gbk_genome_table_data", mock_gbk_genome_table)
    monkeypatch.setattr(add_genome_data, "insert_data", mock_return_none)


    assembly_prot_dict = {
        'Assemb1' : ['prot1'],
        'Assemb2' : ['prot2'],
    }
    gbk_dict = {'prot1': 1, 'prot2': 2}

    connection = None

    db_genome_table_dict = {'Assemb1': 1, 'Assemb2': 2}

    add_genome_data.add_prot_genome_relationships(
        assembly_prot_dict,
        gbk_dict,
        db_genome_table_dict,
        connection,
    )
