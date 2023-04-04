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


from pathlib import Path
import pytest
import shutil
import sys

from argparse import Namespace

from cazy_webscraper.sql.sql_interface.add_data import add_genbank_data

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


def test_add_gbk_seqs(monkeypatch, db_path):
    """Test add_gbk_seqs_to_db()"""
    seq_dict = {
        'gbk1': 'BBBB',
        'gbk2': 'CCCC',
    }
    retrieval_data = "2022-11-08"
    gbk_dict = {
        'gbk1': 1,
        'gbk2': 2,
    }
    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
            seq_update=True,
        )
    }

    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=False, future=True)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)  # allows for calls to session later on when required
    connection = engine.connect()

    def mock_gbk_table_seq_dict(*args, **kwards):
        return {
            'gbk1': {'sequence': None, 'seq_date': None},
            'gbk2': {'sequence': 'AAAAAA', 'seq_date': '2022-11-01'},
        }

    monkeypatch.setattr(add_genbank_data, "get_gbk_table_seq_dict", mock_gbk_table_seq_dict)
    Path("tests/test_outputs/test_outputs_sql").mkdir(exist_ok=True)
    cache_dir = Path("tests/test_outputs/test_outputs_sql/temp_dir")
    cache_dir.mkdir(exist_ok=True, parents=True)

    add_genbank_data.add_gbk_seqs_to_db(
        seq_dict,
        retrieval_data,
        gbk_dict,
        connection,
        cache_dir,
        argsdict['args'],
        unit_test=True,
    )

    shutil.rmtree(cache_dir, ignore_errors=True)
