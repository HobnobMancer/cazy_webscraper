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
"""Tests the get_gtdb_taxs.py from the expand module.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import ArgumentParser, Namespace
from pathlib import Path

from saintBioutils.utilities import logger as saint_logger

from cazy_webscraper.expand.gtdb import get_gtdb_tax
from cazy_webscraper.expand import gtdb
from cazy_webscraper.utilities.parsers import get_gtdb_parser


@pytest.fixture
def archaea_data_file_path():
    _path = Path("tests/test_inputs/test_inputs_gtdb/ar53_taxonomy.tsv.gz")
    return _path


@pytest.fixture
def test_genomes():
    _genomes = [
        'GCF_000979375.1',
        'GCF_11111111.1',
        'GCA_002506415.1',
    ]
    return _genomes


@pytest.fixture
def mock_building_parser(*args, **kwargs):
    parser_args = ArgumentParser(
        prog="get_gtdb_taxs.py",
        usage=None,
        description="Retrieve GTDB taxonomic classifications",
        conflict_handler="error",
        add_help=True,
    )
    return parser_args


def test_parse_gtdb_datafile(archaea_data_file_path, test_genomes):
    result = get_gtdb_tax.get_lineage_data(archaea_data_file_path, test_genomes)
    assert result == {
        'GCF_000979375.1': {'lineage': 'd__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei', 'release': 'ar53_taxonomy.tsv'},
        'GCA_002506415.1': {'lineage': 'd__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei', 'release': 'ar53_taxonomy.tsv'}
    }


def test_gtdb_main_no_genomes(monkeypatch, mock_config_logger, db_path, config_dict):
    
    def mock_parser(*args, **kwargs):
        parser = Namespace(
            database=db_path,
            taxs=['archaea', 'bacteria'],
            archaea_file=None,
            bacterial_file=None,
            classes=None,
            config=None,
            force=False,
            families=None,
            genera=None,
            kingdoms=None,
            log=None,
            nodelete=False,
            nodelete_cache=False,
            nodelete_log=False,
            retries=10,
            sql_echo=False,
            subfamilies=False,
            species=None,
            strains=None,
            timeout=45,
            validate=False,
            verbose=False,
            version=False,
            uniprot_accessions=None,
            update_genom_lineage=False,
            cache_dir=Path("tests/test_output/mock_gtdb_cache"),
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return

    def mock_get_expansion_configuration(*args, **kwards):
        return config_dict, set(), set(), set(), dict(), set()
    
    def mock_empty_dict(*args, **kwards):
        return {}

    def mock_get_links(*args, **kwards):
        return 'url', 'url'
    
    def mock_get_genomes(*args, **kwards):
        return {}, []

    monkeypatch.setattr(get_gtdb_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger, "config_logger", mock_config_logger)
    monkeypatch.setattr(get_gtdb_tax, "make_output_directory", mock_return_none)
    monkeypatch.setattr(get_gtdb_tax, "get_expansion_configuration", mock_get_expansion_configuration)
    monkeypatch.setattr(get_gtdb_tax, "get_gbks_of_interest", mock_return_none)
    monkeypatch.setattr(get_gtdb_tax, "get_genomes", mock_get_genomes)
    # monkeypatch.setattr(get_gtdb_tax, "get_gtdb_data", mock_get_links)
    # monkeypatch.setattr(get_gtdb_tax, "get_lineage_data", mock_empty_dict)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_gtdb_tax.main()
    assert pytest_wrapped_e.type == SystemExit
    


def test_gtdb_main_no_proteins(monkeypatch, mock_config_logger, db_path, config_dict):
    
    def mock_parser(*args, **kwargs):
        parser = Namespace(
            database=db_path,
            taxs=['archaea', 'bacteria'],
            archaea_file=None,
            bacterial_file=None,
            classes=None,
            config=None,
            force=False,
            families=None,
            genera=None,
            kingdoms=None,
            log=None,
            nodelete=False,
            nodelete_cache=False,
            nodelete_log=False,
            retries=10,
            sql_echo=False,
            subfamilies=False,
            species=None,
            strains=None,
            timeout=45,
            validate=False,
            verbose=False,
            version=False,
            uniprot_accessions=None,
            update_genom_lineage=False,
            cache_dir=Path("tests/test_output/mock_gtdb_cache"),
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return

    def mock_get_expansion_configuration(*args, **kwards):
        return config_dict, set(), set(), set(), dict(), set()
    
    def mock_empty_dict(*args, **kwards):
        return {}

    def mock_get_links(*args, **kwards):
        return 'url', 'url'
    
    def mock_get_genomes(*args, **kwards):
        return {}, ['g1', 'g2']

    monkeypatch.setattr(get_gtdb_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger, "config_logger", mock_config_logger)
    monkeypatch.setattr(get_gtdb_tax, "make_output_directory", mock_return_none)
    monkeypatch.setattr(get_gtdb_tax, "get_expansion_configuration", mock_get_expansion_configuration)
    monkeypatch.setattr(get_gtdb_tax, "get_gbks_of_interest", mock_return_none)
    monkeypatch.setattr(get_gtdb_tax, "get_genomes", mock_get_genomes)
    monkeypatch.setattr(get_gtdb_tax, "get_gtdb_data", mock_get_links)
    monkeypatch.setattr(get_gtdb_tax, "get_lineage_data", mock_empty_dict)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_gtdb_tax.main()
    assert pytest_wrapped_e.type == SystemExit


def test_gtdb_main_successful(monkeypatch, mock_config_logger, db_path, config_dict):
    
    def mock_parser(*args, **kwargs):
        parser = Namespace(
            database=db_path,
            taxs=['archaea', 'bacteria'],
            archaea_file=None,
            bacterial_file=None,
            classes=None,
            config=None,
            force=False,
            families=None,
            genera=None,
            kingdoms=None,
            log=None,
            nodelete=False,
            nodelete_cache=False,
            nodelete_log=False,
            retries=10,
            sql_echo=False,
            subfamilies=False,
            species=None,
            strains=None,
            timeout=45,
            validate=False,
            verbose=False,
            version=False,
            uniprot_accessions=None,
            update_genom_lineage=False,
            cache_dir=Path("tests/test_output/mock_gtdb_cache"),
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return

    def mock_get_expansion_configuration(*args, **kwards):
        return config_dict, set(), set(), set(), dict(), set()
    
    def mock_dict(*args, **kwards):
        return {1:1, 2:2}

    def mock_get_links(*args, **kwards):
        return 'url', 'url'
    
    def mock_get_genomes(*args, **kwards):
        return {}, ['g1', 'g2']

    monkeypatch.setattr(get_gtdb_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger, "config_logger", mock_config_logger)
    monkeypatch.setattr(get_gtdb_tax, "make_output_directory", mock_return_none)
    monkeypatch.setattr(get_gtdb_tax, "get_expansion_configuration", mock_get_expansion_configuration)
    monkeypatch.setattr(get_gtdb_tax, "get_gbks_of_interest", mock_return_none)
    monkeypatch.setattr(get_gtdb_tax, "get_genomes", mock_get_genomes)
    monkeypatch.setattr(get_gtdb_tax, "get_gtdb_data", mock_get_links)
    monkeypatch.setattr(get_gtdb_tax, "get_lineage_data", mock_dict)
    monkeypatch.setattr(get_gtdb_tax, "add_gtdb_taxs", mock_return_none)
    monkeypatch.setattr(get_gtdb_tax, "add_genome_gtdb_relations", mock_return_none)
    monkeypatch.setattr(get_gtdb_tax, "closing_message", mock_return_none)

    get_gtdb_tax.main()


def test_getting_download_links_fail(monkeypatch):
    _local_page = Path("tests/test_inputs/test_inputs_gtdb/gtdb_page.html")

    def mock_get_page(*args, **kwards):
        return None, 'error'

    _args = {"args": Namespace(
        retries=10,
    )}

    _output_dir = Path("tests/test_outputs/test_gtdb_output")

    monkeypatch.setattr(gtdb, "get_page", mock_get_page)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        gtdb.get_gtdb_data(_args['args'], _output_dir, True, True)
    assert pytest_wrapped_e.type == SystemExit
