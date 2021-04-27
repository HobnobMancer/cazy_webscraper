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
"""Tests the module utilities which builds the logger and args parser.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace

from scraper.utilities import config_logger, parsers


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


def test_parser():
    """Test building the parser when argsv is None"""
    parsers.build_parser()


def test_parser_arsv():
    """Test building the parser when argsv is not None"""
    parsers.build_parser(["-f"])


# test building the logger


def test_verbose_false(args_v_false):
    """Test build_logger when args.verbose and args.log are false"""
    config_logger(args_v_false["args"])


def test_verbose_true(args_v_true):
    """Test build_logger when args.verbose and args.log are true"""
    config_logger(args_v_true["args"])


# test building genbank_sequencases parser


def test_genbank_seq_parser():
    """Test building the parser when argsv is None"""
    parsers.build_genbank_sequences_parser()


def test_genbank_seq_parser_arsv():
    """Test building the parser when argsv is not None"""
    parsers.build_genbank_sequences_parser(["database", "email"])


# test building genbank_sequencases parser


def test_pdb_struc_parser():
    """Test building the parser when argsv is None"""
    parsers.build_pdb_structures_parser()


def test_pdb_struc_parser_arsv():
    """Test building the parser when argsv is not None"""
    parsers.build_pdb_structures_parser(["database", "pdb"])
