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
"""Tests the script cazy_webscraper.py which coordinates the scraping of CAZy.

These test are intened to be run from the root of the repository using:
pytest -v
"""

from logging import getLogger
import logging
import os

from cazy_webscraper.utilities.parse_configuration.cazy_class_synonym_dict import cazy_synonym_dict
from saintBioutils.utilities import logger as saint_logger
import pytest
import sys

import pandas as pd

from argparse import Namespace, ArgumentParser
from datetime import datetime
from pathlib import Path

from cazy_webscraper import cazy_scraper, sql, closing_message
from cazy_webscraper.sql.sql_interface.add_data import add_cazyme_data
from cazy_webscraper.sql import sql_interface, sql_orm
from cazy_webscraper.utilities import parse_configuration, parsers
from cazy_webscraper.utilities.parsers import cazy_webscraper_parser
from cazy_webscraper.utilities.parse_configuration import cazy_synonym_dict


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


@pytest.fixture
def mock_parser(*args, **kwargs):
    parser = {"args": Namespace(
        config=None,
        classes=None,
        database="fake_database_path",
        ec=None,
        force=False,
        families=None,
        genera=None,
        get_pages=True,
        kingdoms=None,
        log=None,
        nodelete=False,
        output=output_dir,
        retries=10,
        scrape_files=None,
        subfamilies=True,
        species=None,
        strains=None,
        streamline=None,
        timeout=45,
        verbose=False,
    )}
    return parser


@pytest.fixture
def taxonomic_filter_dict():
    """Dict returned from parse_configuration when no tax filters given."""
    taxonomy_filter = {"genera": [], "species": [], "strains": []}
    return taxonomy_filter


@pytest.fixture
def input_dir(test_input_dir):
    dir_path = test_input_dir / "test_inputs_webscraper"
    return dir_path


@pytest.fixture
def output_dir(test_dir):
    path = test_dir / "test_outputs"
    return path


@pytest.fixture
def db_path():
    db_path = "tests/test_inputs/unit_test_database/unit_test_2021-04-27--11-54-58.db"
    return db_path


@pytest.fixture
def logs_dir(output_dir):
    path_ = output_dir / "test_webscraper" / "test_logs"
    return path_


@pytest.fixture
def args_get_cazy_data(logs_dir):
    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
            output=logs_dir,
        )
    }
    return argsdict


@pytest.fixture
def args_get_cazy_data_stdout(logs_dir):
    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
            output=sys.stdout,
        )
    }
    return argsdict


@pytest.fixture
def config_dict():
    configuration_dict = {
        'classes': [],
        "Glycoside Hydrolases (GHs)": ["GH3"],
        'GlycosylTransferases (GTs)': [],
        "Polysaccharide Lyases (PLs)": None,
        'Carbohydrate Esterases (CEs)': [],
        'Auxiliary Activities (AAs)': [],
        'Carbohydrate-Binding Modules (CBMs)': [],
    }
    return configuration_dict


@pytest.fixture
def empty_config_dict():
    configuration_dict = {
        'classes': [],
        "Glycoside Hydrolases (GHs)": [],
        'GlycosylTransferases (GTs)': [],
        "Polysaccharide Lyases (PLs)": [],
        'Carbohydrate Esterases (CEs)': [],
        'Auxiliary Activities (AAs)': [],
        'Carbohydrate-Binding Modules (CBMs)': [],
    }
    return configuration_dict


# test main()


def test_main_citation(mock_building_parser, mock_config_logger, monkeypatch):
    """Test main() with args.citation"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
            skip_ncbi_tax=False,
            cache_dir=None,
            cazy_data=None,
            cazy_synonyms=None,
            classes=None,
            config=None,
            citation=True,
            db_output=None,
            database=None,
            delete_old_relationships=False,
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
        )
        return parser

    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger, "config_logger", mock_config_logger)

    cazy_scraper.main()


def test_main_version(mock_building_parser, mock_config_logger, monkeypatch):
    """Test main() with args.version"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
            skip_ncbi_tax=False,
            cazy_synonyms=None,
            classes=None,
            config=None,
            citation=False,
            db_output=None,
            database=None,
            delete_old_relationships=False,
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
            version=True,
        )
        return parser

    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger,  "config_logger", mock_config_logger)

    cazy_scraper.main()


def test_main_version_argv(mock_building_parser, mock_config_logger, monkeypatch):
    """Test main() with args.version"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
            cazy_synonyms=None,
            skip_ncbi_tax=False,
            classes=None,
            config=None,
            citation=False,
            db_output=None,
            database=None,
            delete_old_relationships=False,
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
            version=True,
        )
        return parser

    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger,  "config_logger", mock_config_logger)

    cazy_scraper.main()


def test_main_double_db(mock_building_parser, mock_config_logger, monkeypatch):
    """Test main() with args.database and args.db_output"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
            cazy_synonyms=None,
            skip_ncbi_tax=False,
            classes=None,
            config=None,
            citation=False,
            db_output=Path('database.db'),
            database=Path('database.db'),
            delete_old_relationships=False,
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
        )
        return parser

    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger,  "config_logger", mock_config_logger)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy_scraper.main()
    assert pytest_wrapped_e.type == SystemExit


def test_main_new_db_exists_force(db_path, mock_building_parser, mock_config_logger, monkeypatch):
    """Test main() with args.db_output and the database already exists"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
            cazy_synonyms=None,
            classes=None,
            config=None,
            citation=False,
            skip_ncbi_tax=False,
            db_output=Path(db_path),  ###
            database=None,
            delete_old_relationships=False,
            force=True,   ###
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
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return
    
    def mock_connect_db(*args, **kwargs):
        return (
            'connection',
            Path('logger_name'),
            Path('cache_dir'),
        )
    
    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger,  "config_logger", mock_config_logger)
    monkeypatch.setattr(cazy_scraper, "connect_to_new_db", mock_connect_db)
    monkeypatch.setattr(cazy_scraper, "connect_existing_db", mock_connect_db)
    monkeypatch.setattr(cazy_scraper, "make_output_directory", mock_return_none)
    monkeypatch.setattr(sql_interface, "log_scrape_in_db", mock_return_none)
    monkeypatch.setattr("cazy_webscraper.cazy_scraper.make_output_directory", mock_return_none)
    monkeypatch.setattr(cazy_scraper, "get_cazy_data", mock_return_none)

    cazy_scraper.main()


def test_main_new_db_exists(db_path, mock_building_parser, mock_config_logger, monkeypatch):
    """Test main() with args.db_output and the database already exists"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=Path('cazy_data.json'),
            cazy_synonyms=None,
            classes=None,
            skip_ncbi_tax=False,
            config=None,
            citation=False,
            db_output=Path(db_path),  ###
            database=None,
            delete_old_relationships=False,
            force=False,   ###
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
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return

    def mock_connect_db(*args, **kwargs):
        return (
            'connection',
            Path('logger_name'),
            Path('cache_dir'),
        )

    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger,  "config_logger", mock_config_logger)
    monkeypatch.setattr(cazy_scraper, "connect_to_new_db", mock_connect_db)
    monkeypatch.setattr(cazy_scraper, "connect_existing_db", mock_connect_db)
    monkeypatch.setattr(cazy_scraper, "make_output_directory", mock_return_none)
    monkeypatch.setattr(sql_interface, "log_scrape_in_db", mock_return_none)
    monkeypatch.setattr("cazy_webscraper.cazy_scraper.make_output_directory", mock_return_none)
    monkeypatch.setattr(cazy_scraper, "get_cazy_data", mock_return_none)

    cazy_scraper.main()


def test_main_database(
    db_path,
    mock_building_parser,
    mock_config_logger,
    config_dict,
    monkeypatch,
):
    """Test main() with args.db_output and the database already exists"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
            cache_dir='cache',
            cazy_data=None,
            cazy_synonyms=None,
            classes=None,
            skip_ncbi_tax=False,
            config=None,
            citation=False,
            db_output=None,
            database=Path(db_path),  ###
            delete_old_relationships=False,
            force=False,
            families=None,
            genera=None,
            kingdoms=None,
            log=Path('log_dir'),
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
        )
        return parser

    def mock_parse_config(*args, **kwards):
        return (
            None,
            config_dict,
            {},
            {'GH'},
            {'CE1'},
            {'Bacteria'},
            {'genus': 'a genus'},
            {'a genus', 'a species'}
        )

    def mock_connect_existing_db(*args, **kwards):
        return (
            'connection',
            'logger_name',
            'cache_dir',
        )

    def mock_return_none(*args, **kwargs):
        return

    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger,  "config_logger", mock_config_logger)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_parse_config)
    monkeypatch.setattr(cazy_scraper, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr(sql.sql_interface, "log_scrape_in_db", mock_return_none)
    monkeypatch.setattr("cazy_webscraper.cazy_scraper.make_output_directory", mock_return_none)
    monkeypatch.setattr(cazy_scraper, "get_cazy_data", mock_return_none)

    cazy_scraper.main()


def test_main_db_output(
    db_path,
    mock_building_parser,
    mock_config_logger,
    config_dict,
    monkeypatch,
):
    """Test main() with args.db_output and the database already exists"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
            cazy_synonyms=None,
            classes=None,
            config=None,
            citation=False,
            db_output=Path(db_path),  ###
            database=None,
            delete_old_relationships=False,
            force=False,
            families=None,
            skip_ncbi_tax=False,
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
        )
        return parser

    def mock_parse_config(*args, **kwards):
        return (
            None,
            config_dict,
            {},
            {'GH'},
            {'CE1'},
            {'Bacteria'},
            {'genus': 'a genus'},
            {'a genus', 'a species'}
        )

    def mock_return_none(*args, **kwargs):
        return

    def mock_connect_db(*args, **kwargs):
        return (
            'connection',
            Path('logger_name'),
            Path('cache_dir'),
        )

    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger,  "config_logger", mock_config_logger)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_parse_config)
    monkeypatch.setattr(cazy_scraper, "connect_to_new_db", mock_connect_db)
    monkeypatch.setattr(cazy_scraper, "connect_existing_db", mock_connect_db)
    monkeypatch.setattr(cazy_scraper, "make_output_directory", mock_return_none)
    monkeypatch.setattr(sql_interface, "log_scrape_in_db", mock_return_none)
    monkeypatch.setattr("cazy_webscraper.cazy_scraper.make_output_directory", mock_return_none)
    monkeypatch.setattr(cazy_scraper, "get_cazy_data", mock_return_none)


    cazy_scraper.main()


def test_closing_message():
    """Test closing_message() from CW __init__"""
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    argsdict = {
        "args": Namespace(
            verbose=False,
        )
    }

    closing_message('cazy_webscraper', start_time, argsdict['args'])


def test_closing_message_verbose():
    """Test closing_message() from CW __init__"""
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    argsdict = {
        "args": Namespace(
            verbose=True,
        )
    }

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        closing_message('cazy_webscraper', start_time, argsdict['args'], early_term=True)
    assert pytest_wrapped_e.type == SystemExit


# # # test get_cazy_data()


def test_get_cazy_data_no_filters(db_path, empty_config_dict, monkeypatch):
    """Test get_cazy_data() when no user defined features are provided"""

    argsdict = {
        "args": Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
            cazy_synonyms=None,
            classes=None,
            config=None,
            citation=False,
            db_output=Path(db_path),  ###
            database=None,
            delete_old_relationships=False,
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
            validate=True,
            verbose=False,
            version=False,
        )
    }

    def mock_return_logger(*args, **kwards):
        return logging.getLogger('mock_logger')

    def mock_return_none(*args, **kwards):
        return

    def mock_tax_dict(*args, **kwargs):
        return {}, {}

    def mock_multiple_taxa(*args, **kwards):
        return [1,2,3,4]

    def mock_replace_taxa(*args, **kwards):
        return [], True

    def mock_cazy_data(*args, **kwargs):
        return {'protein_1': 1, 'protein_2': 2}

    monkeypatch.setattr(cazy_scraper, "build_logger", mock_return_logger)
    monkeypatch.setattr(cazy_scraper, "get_validation_data", mock_return_none)
    monkeypatch.setattr(cazy_scraper, "get_cazy_txt_file_data", mock_multiple_taxa)
    monkeypatch.setattr(cazy_scraper, "parse_all_cazy_data", mock_cazy_data)
    monkeypatch.setattr(cazy_scraper, "identify_multiple_taxa", mock_multiple_taxa)
    monkeypatch.setattr(cazy_scraper, "replace_multiple_tax", mock_replace_taxa)
    monkeypatch.setattr(cazy_scraper, "build_taxa_dict", mock_tax_dict)
    monkeypatch.setattr(add_cazyme_data, "add_kingdoms", mock_return_none)
    monkeypatch.setattr(add_cazyme_data, "add_source_organisms", mock_return_none)
    monkeypatch.setattr(add_cazyme_data, "add_cazy_families", mock_return_none)
    monkeypatch.setattr(add_cazyme_data, "add_genbanks", mock_return_none)
    monkeypatch.setattr(add_cazyme_data, "add_genbank_fam_relationships", mock_return_none)

    cazy_scraper.get_cazy_data(
        "cazy_home_page",
        None,
        cazy_synonym_dict,
        empty_config_dict,
        set(),
        set(),
        set(),
        set(),
        'connection',
        Path('tests/test_outputs/test_webscraper/test_logs'),
        'logger_name',
        'time_stamp',
        argsdict['args'],
    )


def test_term_message():
    """Test closing_message when terminating early"""
    argsdict = {
        "args": Namespace(
            verbose="False",
        )
    }

    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    start_time = pd.to_datetime(start_time)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        closing_message(
            job="Test",
            start_time=start_time,
            args=argsdict["args"],
            early_term=True
        )
    assert pytest_wrapped_e.type == SystemExit
