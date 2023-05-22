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
"""Tests the script expand/uniprot/get_uniprot_data.py.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import pytest
import shutil

from argparse import Namespace, ArgumentParser
from pathlib import Path

from saintBioutils import uniprot
from saintBioutils.utilities import file_io
from saintBioutils.utilities import logger as saint_logger

import cazy_webscraper

from cazy_webscraper import utilities
from cazy_webscraper.cazy_scraper import connect_existing_db
from cazy_webscraper.expand.uniprot import get_uniprot_data
from cazy_webscraper.ncbi.gene_names import get_linked_ncbi_accessions
from cazy_webscraper.sql import sql_interface
from cazy_webscraper.sql.sql_interface.add_data import add_uniprot_data
from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks, get_table_dicts
from cazy_webscraper.utilities.parsers import uniprot_parser



@pytest.fixture
def mock_building_parser(*args, **kwargs):
    parser_args = ArgumentParser(
        prog="cazy_webscraper.py",
        usage=None,
        description="Scrape the CAZy database",
        conflict_handler="error",
        add_help=True,
    )
    return parser_args


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
            config=None,
            classes=None,
            database="fake_database_path",
            ec=True,
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
            update_seq=True,
            subfamilies=True,
            species=None,
            strains=None,
            streamline=None,
            timeout=45,
            uniprot_accessions=None,
            uniprot_batch_size=150,
            uniprot_data=None,
            verbose=False,
            pdb=True,
            skip_uniprot_accessions=None,
            use_uniprot_cache=None,
            email="dummy_email",
            taxonomy=True,
            skip_download=False,
            delete_old_ec_relationships=True,
            delete_old_ecs=True,
            delete_old_pdb_relationships=True,
            delete_old_pdbs=True,
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return

    def mock_connect_existing_db(*args, **kwards):
        return db_connection, None, "cache_dir"

    def mock_get_expansion_configuration(*args, **kwards):
        return config_dict, set(), set(), set(), dict(), set()

    def mock_get_genbank_accessions(*args, **kwards):
        return {1: 1, 2:2, 3:3}
    
    def mock_get_uniprot_data(*args, **kwards):
        return  {'ncbi_acc': {
            'uniprot_acc': 'uniprot_acc',
            'uniprot_entry_id': 'uniprot_entry_id',
            'protein_name': 'protein_name',
            'ec_numbers': {'ec_numbers'},
            'sequence': 'sequence',
            'pdbs': {'all_pdbs'},
        }}

    def mock_uniprot_data(*args, **kwards):
        return {
                1: {'genbank_accession': 1, 'gene_name': 1, 'ec': {1,2,3}, 'pdb': {1,2,3}},
                2: {'gene_name': 'nan', 'ec': {1,2,3}, 'pdb': {1,2,3}}, 3: {'ec': {1,2,3}, 'pdb': {1,2,3}}
            }
    
    monkeypatch.setattr(uniprot_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger, "config_logger", mock_return_logger)
    monkeypatch.setattr(get_uniprot_data, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr("cazy_webscraper.expand.uniprot.get_uniprot_data.make_output_directory", mock_return_none)
    monkeypatch.setattr(get_uniprot_data, "get_expansion_configuration", mock_get_expansion_configuration)
    monkeypatch.setattr(sql_interface, "log_scrape_in_db", mock_return_none)
    monkeypatch.setattr(get_selected_gbks, "get_genbank_accessions", mock_get_genbank_accessions)
    monkeypatch.setattr(get_uniprot_data, "get_uniprot_accessions", mock_get_genbank_accessions)
    monkeypatch.setattr(get_uniprot_data, "add_uniprot_genbank_relationships", mock_return_none)
    monkeypatch.setattr(get_uniprot_data, "get_uniprot_data", mock_get_uniprot_data)
    monkeypatch.setattr(get_uniprot_data, "add_uniprot_accessions", mock_return_none)
    monkeypatch.setattr(get_uniprot_data, "add_ec_numbers", mock_return_none)
    monkeypatch.setattr(get_uniprot_data, "add_genbank_ec_relationships", mock_return_none)
    monkeypatch.setattr(get_uniprot_data, "delete_old_relationships", mock_return_none)
    monkeypatch.setattr(get_uniprot_data, "delete_old_annotations", mock_return_none)
    monkeypatch.setattr(get_uniprot_data, "add_pdb_accessions", mock_return_none)
    monkeypatch.setattr(get_uniprot_data, "add_pdb_gbk_relationships", mock_return_none)
    monkeypatch.setattr(get_uniprot_data, "delete_old_annotations", mock_return_none)
    monkeypatch.setattr(get_uniprot_data, "delete_old_relationships", mock_return_none)
    monkeypatch.setattr(get_uniprot_data, "add_uniprot_taxs", mock_return_none)
    monkeypatch.setattr(cazy_webscraper, "closing_message", mock_return_none)

    monkeypatch.setattr(get_table_dicts, "get_ec_table_dict", mock_get_genbank_accessions)
    monkeypatch.setattr(get_table_dicts, "get_ec_gbk_table_dict", mock_get_genbank_accessions)
    monkeypatch.setattr(get_table_dicts, "get_pdb_table_dict", mock_get_genbank_accessions)
    monkeypatch.setattr(get_table_dicts, "get_gbk_pdb_table_dict", mock_get_genbank_accessions)

    output = test_dir / "test_outputs" / "test_outputs_uniprot"
    output.mkdir(parents=True, exist_ok=True)
    get_uniprot_data.main()
    shutil.rmtree((test_dir / "test_outputs" / "test_outputs_uniprot"))
    output.mkdir(parents=True, exist_ok=True)
