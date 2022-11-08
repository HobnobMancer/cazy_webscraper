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
"""Tests expand.extract_seqs.extract_db_seqs.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


from argparse import Namespace, ArgumentParser
from datetime import datetime
from pathlib import Path

import pytest

from saintBioutils.utilities import logger as saint_logger
from sqlalchemy.exc import IntegrityError

from cazy_webscraper.expand.extract_seqs import extract_db_seqs


def test_validate_opts_no_output():
    argsdict = {
        "args": Namespace(
            source=['genbank', 'uniprot'],
            blastdb=None,
            fasta_dir=None,
            fasta_file=None,
        )
    }

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        extract_db_seqs.validate_user_options(argsdict['args'])
    assert pytest_wrapped_e.type == SystemExit


def test_validate_opts_all_outputs():
    argsdict = {
        "args": Namespace(
            source=['genbank', 'uniprot'],
            blastdb='blast_db',
            fasta_dir='fasta_dir',
            fasta_file='fasta_file',
        )
    }

    extract_db_seqs.validate_user_options(argsdict['args'])


def test_get_seqs():
    gbk_seq_dict = {
        'gbk_acc_1': {'sequence': 'abc'},
        'gbk_acc_2': {'sequence': 'abc'},
    }
    gbk_dict = {
        'gbk_acc_1': {'db': 'genbank', 'seq': 'abc'},
        'gbk_acc_2': {'db': 'genbank', 'seq': 'abc'},
    }
    assert {
        'gbk_acc_1': {'db': 'GenBank', 'seq': 'abc'},
        'gbk_acc_2': {'db': 'GenBank', 'seq': 'abc'},
    } == extract_db_seqs.get_genbank_sequences(gbk_seq_dict, gbk_dict)


def test_get_uniprot_seqs():
    uniprot_table_dict = {
        'uni_acc_1': {'name': 'name', 'genbank_id': 1, 'seq': 'abc', 'seq_date': 'date'},
        'uni_acc_2': {'name': 'name', 'genbank_id': 2, 'seq': 'abc', 'seq_date': 'date'},
    }
    gbk_dict = {
        'gbk_acc_1': 1,
        'gbk_acc_2': 2,
    }
    assert {
        'uni_acc_1': {'db': 'UniProt', 'seq': 'abc'},
        'uni_acc_2': {'db': 'UniProt', 'seq': 'abc'},
    } == extract_db_seqs.get_uniprot_sequences(uniprot_table_dict, gbk_dict)


