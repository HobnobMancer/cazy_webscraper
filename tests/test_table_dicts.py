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
"""Tests the sql.sql_interface.get_data.get_table_dicts.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace

from cazy_webscraper.sql import sql_orm
from cazy_webscraper.sql.sql_interface.get_data import get_table_dicts


def test_get_ec_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_ec_table_dict(db_connection)


def test_get_ec_gbk_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_ec_gbk_table_dict(db_connection)
    

def test_get_gbk_ec_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_gbk_ec_table_dict(db_connection)


def test_get_fam_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_fams_table_dict(db_connection)


def test_get_gbk_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_gbk_table_dict(db_connection)


def test_get_gbk_no_tax_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_no_tax_gbk_table_dict(db_connection)


def test_get_gbk_seq_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_gbk_table_seq_dict(db_connection)


def test_get_gbk_fam_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_gbk_fam_table_dict(db_connection)


def test_get_kingdom_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_kingdom_table_dict(db_connection)


def test_get_pdb_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_pdb_table_dict(db_connection)


def test_get_gbk_pdb_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_gbk_pdb_table_dict(db_connection)


def test_get_taxs_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_taxs_table_dict(db_connection)


def test_get_uniprot_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_uniprot_table_dict(db_connection)


def test_get_gbk_kingdom_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_gbk_kingdom_dict(db_connection)


def test_get_ncbi_tax_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_ncbi_tax_table(db_connection)


def test_get_gtdb_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_gtdb_table_dict(db_connection)


def test_get_gtdb_table(db_path):
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_table_dicts.get_genome_table(db_connection)
