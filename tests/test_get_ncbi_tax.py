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
"""Tests expand.taxonomy.get_ncbi_taxs.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


from argparse import Namespace, ArgumentParser
from datetime import datetime
from pathlib import Path

import pytest

from saintBioutils.utilities import logger as saint_logger
from sqlalchemy.exc import IntegrityError

from cazy_webscraper.expand.genbank.taxonomy import get_ncbi_taxs
from cazy_webscraper.sql import sql_interface
from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks
from cazy_webscraper.sql.sql_interface.add_data import add_ncbi_tax_data
from cazy_webscraper.utilities.parsers import tax_ncbi_parser


def test_main(
    mock_building_parser,
    mock_return_logger,
    config_dict,
    db_connection,
    monkeypatch,
    test_dir,
):
    """Test main()"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            cache_dir=(test_dir / "test_outputs" / "test_outputs_uniprot"),
            nodelete_cache=False,
            batch_size=150,
            config=None,
            classes=None,
            database="fake_database_path",
            ec_filter=None,
            force=False,
            families=None,
            genbank_accessions=None,
            genera=None,
            get_pages=True,
            kingdoms=None,
            log=None,
            nodelete=False,
            output=None,
            retries=10,
            sequence=True,
            seq_update=True,
            subfamilies=True,
            species=None,
            strains=None,
            uniprot_accessions=None,
            use_lineage_cache=None,
            verbose=False,
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return

    def mock_connect_existing_db(*args, **kwards):
        return db_connection, None, "cache_dir"

    def mock_get_expansion_configuration(*args, **kwards):
        return config_dict, set(), set(), set(), dict(), set()

    def mock_get_genbank_accessions(*args, **kwards):
        return {1: 1, 2: 2, 3: 3}

    def mock_get_ncbi_data(*args, **kwards):
        return (
            {
                1: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                2: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                3: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}}
            },
            {1, 2, 3},
        )

    monkeypatch.setattr(tax_ncbi_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger, "config_logger", mock_return_logger)
    monkeypatch.setattr(get_ncbi_taxs, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr("cazy_webscraper.expand.uniprot.get_uniprot_data.make_output_directory", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "get_expansion_configuration", mock_get_expansion_configuration)
    monkeypatch.setattr(sql_interface, "log_scrape_in_db", mock_return_none)
    # not using cached lineages
    monkeypatch.setattr(get_ncbi_taxs, "get_db_proteins", mock_get_genbank_accessions)
    monkeypatch.setattr(get_ncbi_taxs, "get_ncbi_ids", mock_get_ncbi_data)
    monkeypatch.setattr(get_ncbi_taxs, "get_lineage_protein_data", mock_get_genbank_accessions)
    # mock adding data to the local CAZyme database
    monkeypatch.setattr(get_ncbi_taxs, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "closing_message", mock_return_none)

    get_ncbi_taxs.main()


def test_main_using_lineage_cache(
    mock_building_parser,
    mock_return_logger,
    config_dict,
    db_connection,
    monkeypatch,
    test_dir,
):
    """Test main()"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            cache_dir=(test_dir / "test_outputs" / "test_outputs_uniprot"),
            nodelete_cache=False,
            batch_size=150,
            config=None,
            classes=None,
            database="fake_database_path",
            ec_filter=None,
            force=False,
            families=None,
            genbank_accessions=None,
            genera=None,
            get_pages=True,
            kingdoms=None,
            log=None,
            nodelete=False,
            output=None,
            retries=10,
            sequence=True,
            seq_update=True,
            subfamilies=True,
            species=None,
            strains=None,
            uniprot_accessions=None,
            use_lineage_cache=Path("tests/test_inputs/test_inputs_ncbi_tax/test_lineage_cache.json"),
            verbose=False,
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return

    def mock_connect_existing_db(*args, **kwards):
        return db_connection, None, "cache_dir"

    def mock_get_expansion_configuration(*args, **kwards):
        return config_dict, set(), set(), set(), dict(), set()

    def mock_get_genbank_accessions(*args, **kwards):
        return {1: 1, 2: 2, 3: 3}

    def mock_get_ncbi_data(*args, **kwards):
        return (
            {
                1: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                2: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                3: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}}
            },
            {1, 2, 3},
        )

    monkeypatch.setattr(tax_ncbi_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger, "config_logger", mock_return_logger)
    monkeypatch.setattr(get_ncbi_taxs, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr("cazy_webscraper.expand.uniprot.get_uniprot_data.make_output_directory", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "get_expansion_configuration", mock_get_expansion_configuration)
    monkeypatch.setattr(sql_interface, "log_scrape_in_db", mock_return_none)
    # not using cached lineages
    monkeypatch.setattr(get_ncbi_taxs, "get_db_proteins", mock_get_genbank_accessions)
    monkeypatch.setattr(get_ncbi_taxs, "get_ncbi_ids", mock_get_ncbi_data)
    monkeypatch.setattr(get_ncbi_taxs, "get_lineage_protein_data", mock_get_genbank_accessions)
    # mock adding data to the local CAZyme database
    monkeypatch.setattr(get_ncbi_taxs, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "closing_message", mock_return_none)

    get_ncbi_taxs.main()


def test_main_using_lineage_cache_fails(
    mock_building_parser,
    mock_return_logger,
    config_dict,
    db_connection,
    monkeypatch,
    test_dir,
):
    """Test main() when cannot find linaege cache file"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            cache_dir=(test_dir / "test_outputs" / "test_outputs_uniprot"),
            nodelete_cache=False,
            batch_size=150,
            config=None,
            classes=None,
            database="fake_database_path",
            ec_filter=None,
            force=False,
            families=None,
            genbank_accessions=None,
            genera=None,
            get_pages=True,
            kingdoms=None,
            log=None,
            nodelete=False,
            output=None,
            retries=10,
            sequence=True,
            seq_update=True,
            subfamilies=True,
            species=None,
            strains=None,
            uniprot_accessions=None,
            use_lineage_cache=Path("tests/test_inputs/test_inputs_ncbi_tax/test_lineage_cache_NOT_EXIST.json"),
            verbose=False,
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return

    def mock_connect_existing_db(*args, **kwards):
        return db_connection, None, "cache_dir"

    def mock_get_expansion_configuration(*args, **kwards):
        return config_dict, set(), set(), set(), dict(), set()

    def mock_get_genbank_accessions(*args, **kwards):
        return {1: 1, 2: 2, 3: 3}

    def mock_get_ncbi_data(*args, **kwards):
        return (
            {
                1: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                2: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                3: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}}
            },
            {1, 2, 3},
        )

    monkeypatch.setattr(tax_ncbi_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger, "config_logger", mock_return_logger)
    monkeypatch.setattr(get_ncbi_taxs, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr("cazy_webscraper.expand.uniprot.get_uniprot_data.make_output_directory", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "get_expansion_configuration", mock_get_expansion_configuration)
    monkeypatch.setattr(sql_interface, "log_scrape_in_db", mock_return_none)
    # not using cached lineages
    monkeypatch.setattr(get_ncbi_taxs, "get_db_proteins", mock_get_genbank_accessions)
    monkeypatch.setattr(get_ncbi_taxs, "get_ncbi_ids", mock_get_ncbi_data)
    monkeypatch.setattr(get_ncbi_taxs, "get_lineage_protein_data", mock_get_genbank_accessions)
    # mock adding data to the local CAZyme database
    monkeypatch.setattr(get_ncbi_taxs, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "closing_message", mock_return_none)

    get_ncbi_taxs.main()
