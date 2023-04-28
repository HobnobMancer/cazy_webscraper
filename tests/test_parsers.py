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
"""Tests the module utilities which builds the args parser.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace

from cazy_webscraper.utilities.parsers import (
    api_parser,
    cazy_webscraper_parser,
    extract_seq_parser,
    gbk_seq_parser,
    pdb_strctre_parser,
    uniprot_parser,
    tax_ncbi_parser,
    get_genomes_parser,

)


@pytest.fixture
def args_v_false():
    args_dict = {
        "args": Namespace(
            verbose=False,
            log=None,
        )
    }
    return args_dict


@pytest.fixture
def args_v_true():
    args_dict = {
        "args": Namespace(
            verbose=True,
            log="tests/test_outputs/test_log",
        )
    }
    return args_dict


# test building the parser


# test cazy_webscraper parser

def test_parser_cw():
    """Test building the parser when argsv is None"""
    cazy_webscraper_parser.build_parser()


def test_parser_cw_arsv():
    """Test building the parser when argsv is not None"""
    cazy_webscraper_parser.build_parser(["dummy_email"])


def test_parser_extract():
    """Test building the parser when argsv is None"""
    extract_seq_parser.build_parser()


def test_parser_arsv_extract():
    """Test building the parser when argsv is not None"""
    extract_seq_parser.build_parser(["dummy_email", "genbank", "uniprot"])


def test_parser_api():
    """Test building the parser when argsv is None"""
    api_parser.build_parser()


def test_parser_arsv_api():
    """Test building the parser when argsv is not None"""
    api_parser.build_parser(["db", "csv"])


def test_parser_gbk():
    """Test building the parser when argsv is None"""
    gbk_seq_parser.build_parser()


def test_parser_arsv_gbk():
    """Test building the parser when argsv is not None"""
    gbk_seq_parser.build_parser(["db", "dummy_email"])


def test_parser_pdb():
    """Test building the parser when argsv is None"""
    pdb_strctre_parser.build_parser()


def test_parser_arsv_pdb():
    """Test building the parser when argsv is not None"""
    pdb_strctre_parser.build_parser(["dummy_email", "pdb"])


def test_parser_uniprot():
    """Test building the parser when argsv is None"""
    uniprot_parser.build_parser()


def test_parser_arsv_uniprot():
    """Test building the parser when argsv is not None"""
    uniprot_parser.build_parser(["dbpath"])


def test_parser_ncbi_tax():
    """Test building the parser when argsv is None"""
    tax_ncbi_parser.build_parser()


def test_parser_arsv_ncbi_tax():
    """Test building the parser when argsv is not None"""
    tax_ncbi_parser.build_parser(["db", "dummy_email"])


def test_parser_genomes():
    """Test building the parser when argsv is None"""
    get_genomes_parser.build_parser()


def test_parser_arsv_genomes():
    """Test building the parser when argsv is not None"""
    get_genomes_parser.build_parser(["db", "dummy_email"])