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
"""Tests sql.sql_interface.get_data.get_records.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


from argparse import Namespace, ArgumentParser
from datetime import datetime
from pathlib import Path

import pytest

from cazy_webscraper.sql.sql_interface.get_data import get_records


def test_get_user_gbk():
    argsdict = {"args": Namespace(
        genbank_accessions="tests/test_inputs/test_inputs_sql_interface/test_accs.txt",
    )}
    # test_gbk

    gbk_table_dict = {'test_gbk': 1, 'gbk_acc': 2}

    assert {'gbk_acc': 2} == get_records.get_user_genbank_sequences(gbk_table_dict, argsdict['args'])


def test_get_user_gbk_fail():
    argsdict = {"args": Namespace(
        genbank_accessions="tests/test_inputs/test_inputs_sql_interface/test_FAIL_accs.txt",
    )}
    # test_gbk

    gbk_table_dict = {'test_gbk': 1, 'gbk_acc': 2}

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_records.get_user_genbank_sequences(gbk_table_dict, argsdict['args'])
    assert pytest_wrapped_e.type == SystemExit


def test_get_user_uniprot_fail():
    argsdict = {"args": Namespace(
        uniprot_accessions="tests/test_inputs/test_inputs_sql_interface/test_FAIL_accs.txt",
    )}
    # test_gbk

    gbk_table_dict = {'test_gbk': 1, 'gbk_acc': 2}

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_records.get_user_uniprot_sequences(gbk_table_dict, {}, argsdict['args'])
    assert pytest_wrapped_e.type == SystemExit


def test_get_user_uniprot():
    argsdict = {"args": Namespace(
        uniprot_accessions="tests/test_inputs/test_inputs_sql_interface/test_accs.txt",
    )}
    # test_gbk

    gbk_table_dict = {'test_gbk_': '1', 'gbk_acc_': '2'}
    uni_table_dict = {'test_acc': {'genbank_id': '1'}, 'gbk_acc': {'genbank_id': '2'}}

    assert {'gbk_acc_': '2', 'test_gbk_': '1'} == get_records.get_user_uniprot_sequences(gbk_table_dict, uni_table_dict, argsdict['args'])
