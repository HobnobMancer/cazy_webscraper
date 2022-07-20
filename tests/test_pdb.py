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
"""Tests the expand.pdb for downloading structure files from RSCB PDB

These test are intened to be run from the root of the repository using:
pytest -v
"""


from pathlib import Path

import logging
import pytest

from argparse import Namespace, ArgumentParser

from cazy_webscraper import connect_existing_db
from cazy_webscraper.expand.pdb import get_pdb_structures
from cazy_webscraper.utilities import parse_configuration
from cazy_webscraper.utilities.parsers import pdb_strctre_parser


def test_main_no_pdb(monkeypatch):
    """Test main() when an the database file cannot be found."""

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_pdb_structures.py",
            usage=None,
            description="Download structure files from RSCB PDB",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            database=Path("--"),
            verbose=False,
            log=None,
            force=False,
            nodelete=False,
            outdir=None,
            genbank_accessions=None,
            uniprot_accessions=None,
        )
        return parser

    def mock_no_return(*args, **kwargs):
        return

    def mock_config(*args, **kwargs):
        return None, set(), set(), set(), set(), set()

    def mock_db_connection(*args, **kwards):
        return "connection", "logger", "cache_dir"
    
    def mock_selected_pdb(*args, **kwards):
        return []

    def mock_gbk_table(*args, **kwards):
        return {}

    monkeypatch.setattr(pdb_strctre_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(parse_configuration, "get_expansion_configuration", mock_config)
    monkeypatch.setattr(get_pdb_structures, "make_output_directory", mock_no_return)
    monkeypatch.setattr(get_pdb_structures, "connect_existing_db", mock_db_connection)
    monkeypatch.setattr(get_pdb_structures, "get_pdb_accessions", mock_selected_pdb)
    monkeypatch.setattr(get_pdb_structures, "get_gbk_table_dict", mock_gbk_table)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_pdb_structures.main()
    assert pytest_wrapped_e.type == SystemExit


def test_main(monkeypatch):
    """Test main() when an the database file cannot be found."""

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_pdb_structures.py",
            usage=None,
            description="Download structure files from RSCB PDB",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            database=Path("--"),
            verbose=False,
            log=None,
            force=False,
            nodelete=False,
            outdir=None,
            genbank_accessions="",
            uniprot_accessions="",
            cache_dir=Path("tests/test_outputs/test_outputs_pdb"),
            nodelete_cache=False,
        )
        return parser

    def mock_no_return(*args, **kwargs):
        return

    def mock_config(*args, **kwargs):
        return None, set(), set(), set(), set(), set()

    def mock_db_connection(*args, **kwards):
        return "connection", "logger", "cache_dir"

    def mock_selected_pdb(*args, **kwards):
        return [1, 2, 3]

    def mock_gbk_table(*args, **kwards):
        return {}

    monkeypatch.setattr(pdb_strctre_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(parse_configuration, "get_expansion_configuration", mock_config)
    monkeypatch.setattr(get_pdb_structures, "make_output_directory", mock_no_return)
    monkeypatch.setattr(get_pdb_structures, "connect_existing_db", mock_db_connection)
    monkeypatch.setattr(get_pdb_structures, "get_pdb_accessions", mock_selected_pdb)
    monkeypatch.setattr(get_pdb_structures, "download_pdb_structures", mock_no_return)
    monkeypatch.setattr(get_pdb_structures, "get_gbk_table_dict", mock_gbk_table)
    monkeypatch.setattr(get_pdb_structures, "get_user_genbank_sequences", mock_gbk_table)
    monkeypatch.setattr(get_pdb_structures, "get_uniprot_table_dict", mock_gbk_table)
    monkeypatch.setattr(get_pdb_structures, "get_user_uniprot_sequences", mock_gbk_table)

    get_pdb_structures.main()
