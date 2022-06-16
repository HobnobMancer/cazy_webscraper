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
"""Tests sql.sql_interface.get_data.get_taxonomies.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


from argparse import Namespace, ArgumentParser
from datetime import datetime
from pathlib import Path

import pytest

from cazy_webscraper.sql import sql_orm
from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks, get_taxonomies
from cazy_webscraper.sql.sql_orm import (
    Ec,
    Genbank,
    Session,
    Taxonomy,
    Uniprot,
)


def test_get_tax_user_acc(monkeypatch):
    argsdict = {"args": Namespace(
        genbank_accessions="tests/test_inputs/test_inputs_sql_interface/test_accs.txt",
        uniprot_accessions="tests/test_inputs/test_inputs_sql_interface/test_accs.txt",
    )}

    def mock_get_tax(*args, **kwards):
        return [1, 2, 3]

    def mock_get_table_dicts(*args, **kwards):
        return {}, {}

    monkeypatch.setattr(get_taxonomies, "get_taxs_for_user_gbks", mock_get_tax)
    monkeypatch.setattr(get_taxonomies, "get_taxs_for_uniprots", mock_get_tax)
    monkeypatch.setattr(get_taxonomies, "get_uni_gbk_tax_dict", mock_get_table_dicts)

    get_taxonomies.get_taxonomies(set(), set(), {}, set(), set(), 'connection', argsdict['args'])


def test_get_tax_db_tax(monkeypatch):
    argsdict = {"args": Namespace(
        genbank_accessions=None,
        uniprot_accessions=None,
    )}

    def mock_get_tax(*args, **kwards):
        return [1, 2, 3]

    monkeypatch.setattr(get_taxonomies, "get_taxs_for_user_gbks", mock_get_tax)
    monkeypatch.setattr(get_taxonomies, "get_taxs_for_uniprots", mock_get_tax)
    monkeypatch.setattr(get_taxonomies, "get_filtered_taxs", mock_get_tax)

    get_taxonomies.get_taxonomies(set(), set(), {}, set(), set(), 'connection', argsdict['args'])


def test_get_uni_gbk_dict(db_path):
    argsdict = {"args": Namespace(
        sql_echo=True,
        uniprot_accessions=None,
    )}
    db_connection = sql_orm.get_db_connection(db_path, argsdict["args"], False)

    assert ({}, {}) == get_taxonomies.get_uni_gbk_tax_dict(db_connection)


def test_get_user_gbks_fail():
    argsdict = {"args": Namespace(
        sql_echo=True,
        uniprot_accessions=None,
        genbank_accessions="tests/test_inputs/test_inputs_sql_interface/test_accs_FAIL.txt",
    )}

    gbk_table_dict = {'test_gbk': 1, 'gbk_acc': 2}

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_taxonomies.get_taxs_for_user_gbks(gbk_table_dict, argsdict['args'])
    assert pytest_wrapped_e.type == SystemExit


def test_get_user_gbks():
    argsdict = {"args": Namespace(
        sql_echo=True,
        uniprot_accessions=None,
        genbank_accessions="tests/test_inputs/test_inputs_sql_interface/test_accs.txt",
    )}

    gbk_table_dict = {'test_gbk': 1, 'gbk_acc': 2}

    assert [2] == get_taxonomies.get_taxs_for_user_gbks(gbk_table_dict, argsdict['args'])
