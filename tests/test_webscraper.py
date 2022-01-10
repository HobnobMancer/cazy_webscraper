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
"""Tests the script cazy_webscraper.py which coordinates the scraping of CAZy.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import os
from cazy_webscraper.utilities.parse_configuration.cazy_class_synonym_dict import cazy_synonym_dict
import pytest
import sys

from argparse import Namespace, ArgumentParser


from cazy_webscraper import cazy_scraper, crawler, sql, utilities
from cazy_webscraper.cazy import parse_all_cazy_data, parse_cazy_data_with_filters
from cazy_webscraper.sql import sql_interface, sql_orm
from cazy_webscraper.utilities import file_io, parse_configuration, parsers
from cazy_webscraper.utilities.parsers import cazy_webscraper_parser


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
    parser = Namespace(
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
    )
    return parser


@pytest.fixture
def mock_config_logger(*args, **kwargs):
    return


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
        "Glycoside Hydrolases (GHs)": ["GH3"],
        "Polysaccharide Lyases (PLs)": None,
    }
    return configuration_dict


# test main()


def test_main_citation(mock_building_parser, mock_config_logger, monkeypatch):
    """Test main() with args.citation"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
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
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)

    cazy_scraper.main()


def test_main_version(mock_building_parser, mock_config_logger, monkeypatch):
    """Test main() with args.version"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
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
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)

    cazy_scraper.main()


def test_main_version_argv(mock_building_parser, mock_config_logger, monkeypatch):
    """Test main() with args.version"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
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
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)

    cazy_scraper.main()


def test_main_double_db(mock_building_parser, mock_config_logger, monkeypatch):
    """Test main() with args.database and args.db_output"""

    argv = [
        "dummy@domain.com",
        None,
        None,
        None,
        None,
        None,
        False,
        'database',
        'database',
        False,
        False,
        None,
        None,
        None,
        None,
        False,
        False,
        False,
        10,
        False,
        False,
        None,
        None,
        45,
        False,
        False,
        True,
    ]
    def mock_parser(*args, **kwargs):
        parser = Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
            cazy_synonyms=None,
            classes=None,
            config=None,
            citation=False,
            db_output='database',
            database='database',
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
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)

    cazy_scraper.main(argv=argv)


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
            db_output=db_path,  ###
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
            version=True,
        )
        return parser

    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)

    cazy_scraper.main()


def test_main_new_db_exists(db_path, mock_building_parser, mock_config_logger, monkeypatch):
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
            db_output=db_path,  ###
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
            version=True,
        )
        return parser

    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)

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
            config=None,
            citation=False,
            db_output=None,
            database=db_path,  ###
            delete_old_relationships=False,
            force=False,
            families=None,
            genera=None,
            kingdoms=None,
            log='log_dir',
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

    def mock_log_scrape_in_db(*args, **kwargs):
        return

    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_parse_config)
    monkeypatch.setattr(cazy_scraper, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr(sql.sql_interface, "log_scrape_in_db", mock_log_scrape_in_db)

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
            db_output=db_path,  ###
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

    def mock_log_scrape_in_db(*args, **kwargs):
        return

    monkeypatch.setattr(cazy_webscraper_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_parse_config)
    monkeypatch.setattr(cazy_scraper, "connect_to_new_db", mock_connect_existing_db)
    monkeypatch.setattr(sql.sql_interface, "log_scrape_in_db", mock_log_scrape_in_db)

    cazy_scraper.main()






# # # test get_cazy_data()


# # def test_get_cazy_data_no_fam_urls(
# #     cazy_home_url,
# #     cazy_dictionary,
# #     config_dict,
# #     time_stamp,
# #     args_get_cazy_data,
# #     logs_dir,
# #     monkeypatch,
# #     null_logger
# # ):
# #     """Test get_cazy_data() when no family URLS are retrieved."""
# #     os.makedirs(logs_dir, exist_ok=True)

# #     def mock_get_classes(*args, **kwargs):
# #         class1 = crawler.CazyClass("test_class", "test_class_url.html", 0)
# #         return [class1]

# #     def mock_get_families(*args, **kwargs):
# #         return None, "test error message", ["test_url1", "test_url2"]

# #     def mock_logger(*args, **kwargs):
# #         return null_logger

# #     monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
# #     monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_families)
# #     monkeypatch.setattr(utilities, "build_logger", mock_logger)

# #     cazy_scraper.get_cazy_data(
# #         cazy_home=cazy_home_url,
# #         excluded_classes=None,
# #         config_dict=config_dict,
# #         cazy_dict=cazy_dictionary,
# #         taxonomy_filters=set(),
# #         kingdoms="all",
# #         ec_filters=[],
# #         time_stamp="timestamp",
# #         session="session_representative",
# #         args=args_get_cazy_data["args"],
# #     )
# #     file_io.make_output_directory(logs_dir, True, False)


# # def test_get_cazy_data_no_all(
# #     time_stamp,
# #     cazy_home_url,
# #     cazy_dictionary,
# #     args_get_cazy_data,
# #     logs_dir,
# #     monkeypatch,
# # ):
# #     """Test get_cazy_data() when no kingdoms are specified and config_dict is None."""
# #     # prepare dir for log files
# #     os.makedirs(logs_dir, exist_ok=True)

# #     fam1 = crawler.Family("test_fam", "test_class", "test_url")

# #     def mock_get_classes(*args, **kwargs):
# #         class1 = crawler.CazyClass(
# #             name="test_class",
# #             url="test_class_url.html",
# #             tries=0,
# #         )
# #         return [class1]

# #     def mock_get_fam_urls(*args, **kwargs):
# #         return [fam1], "error message", ["in", "cor", "rect", "urls"]

# #     def mock_parse_family(*args, **kwargs):
# #         return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"], ["format error"], "session"

# #     monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
# #     monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_fam_urls)
# #     monkeypatch.setattr(scrape_all, "parse_family_via_all_pages", mock_parse_family)

# #     cazy_scraper.get_cazy_data(
# #         cazy_home=cazy_home_url,
# #         excluded_classes=None,
# #         config_dict=None,
# #         cazy_dict=None,
# #         taxonomy_filters=set(),
# #         kingdoms="all",
# #         ec_filters=[],
# #         time_stamp="timestamp",
# #         session="session_representative",
# #         args=args_get_cazy_data["args"],
# #     )
# #     file_io.make_output_directory(logs_dir, True, False)


# # def test_get_cazy_data_no_config_dict_kingdom(
# #     time_stamp,
# #     cazy_home_url,
# #     cazy_dictionary,
# #     args_get_cazy_data,
# #     logs_dir,
# #     monkeypatch,
# # ):
# #     """Test get_cazy_data() when kingdoms are specified and config_dict is None."""
# #     # prepare dir for log files
# #     os.makedirs(logs_dir, exist_ok=True)

# #     fam1 = crawler.Family("test_fam", "test_class", "test_url")

# #     def mock_get_classes(*args, **kwargs):
# #         class1 = crawler.CazyClass(
# #             name="test_class",
# #             url="test_class_url.html",
# #             tries=0,
# #             failed_families={fam1: 0},
# #         )
# #         return [class1]

# #     def mock_parse_family(*args, **kwargs):
# #         return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"], ["format error"], "session"

# #     monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
# #     monkeypatch.setattr(scrape_by_kingdom, "parse_family_by_kingdom", mock_parse_family)

# #     cazy_scraper.get_cazy_data(
# #         cazy_home=cazy_home_url,
# #         excluded_classes=None,
# #         config_dict=None,
# #         cazy_dict=None,
# #         taxonomy_filters=set(),
# #         kingdoms=["Bacteria", "Viruses"],
# #         ec_filters=[],
# #         time_stamp="timestamp",
# #         session="session_representative",
# #         args=args_get_cazy_data["args"],
# #     )
# #     file_io.make_output_directory(logs_dir, True, False)


# # def test_get_cazy_data_config_data_all(
# #     time_stamp,
# #     cazy_home_url,
# #     cazy_dictionary,
# #     args_get_cazy_data,
# #     logs_dir,
# #     monkeypatch,
# # ):
# #     """Test get_cazy_data() when no kingdoms are specified and configuration given."""
# #     # prepare dir for log files
# #     os.makedirs(logs_dir, exist_ok=True)

# #     fam1 = crawler.Family("GH3_1", "test_class", "test_url")

# #     config_dict = {"Glycoside Hydrolases": ["GH3"]}

# #     def mock_get_classes(*args, **kwargs):
# #         class1 = crawler.CazyClass(
# #             name="Glycoside Hydrolases",
# #             url="test_class_url.html",
# #             tries=0,
# #             failed_families={fam1: 0},
# #         )
# #         return [class1]

# #     def mock_parse_family(*args, **kwargs):
# #         return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"], ["format error"], "session"

# #     monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
# #     monkeypatch.setattr(scrape_all, "parse_family_via_all_pages", mock_parse_family)

# #     cazy_scraper.get_cazy_data(
# #         cazy_home=cazy_home_url,
# #         excluded_classes=None,
# #         config_dict=config_dict,
# #         cazy_dict=cazy_dictionary,
# #         taxonomy_filters=set(),
# #         kingdoms="all",
# #         ec_filters=[],
# #         time_stamp="timestamp",
# #         session="session_representative",
# #         args=args_get_cazy_data["args"],
# #     )
# #     file_io.make_output_directory(logs_dir, True, False)


# # def test_get_cazy_data_config_data_kingdom(
# #     time_stamp,
# #     cazy_home_url,
# #     cazy_dictionary,
# #     args_get_cazy_data,
# #     logs_dir,
# #     monkeypatch,
# # ):
# #     """Test get_cazy_data() when kingdoms are specified and configuration given."""
# #     # prepare dir for log files
# #     os.makedirs(logs_dir, exist_ok=True)

# #     fam1 = crawler.Family("GH1", "test_class", "test_url")

# #     config_dict = {"Glycoside Hydrolases": ["GH1"]}

# #     def mock_get_classes(*args, **kwargs):
# #         class1 = crawler.CazyClass(
# #             name="Glycoside Hydrolases",
# #             url="test_class_url.html",
# #             tries=0,
# #             failed_families={fam1: 0},
# #         )
# #         return [class1]

# #     def mock_parse_family(*args, **kwargs):
# #         return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"], ["format error"], {}

# #     monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
# #     monkeypatch.setattr(scrape_by_kingdom, "parse_family_by_kingdom", mock_parse_family)

# #     cazy_scraper.get_cazy_data(
# #         cazy_home=cazy_home_url,
# #         excluded_classes=None,
# #         config_dict=config_dict,
# #         cazy_dict=cazy_dictionary,
# #         taxonomy_filters=set(),
# #         kingdoms=["Bacteria", "Viruses"],
# #         ec_filters=[],
# #         time_stamp="timestamp",
# #         session={},
# #         args=args_get_cazy_data["args"],
# #     )
# #     file_io.make_output_directory(logs_dir, True, False)


# # def test_get_cazy_data_config_data_kingdom_stdout(
# #     time_stamp,
# #     cazy_home_url,
# #     cazy_dictionary,
# #     args_get_cazy_data_stdout,
# #     logs_dir,
# #     monkeypatch,
# # ):
# #     """Test get_cazy_data() when kingdoms are specified and configuration given."""
# #     # prepare dir for log files
# #     os.makedirs(logs_dir, exist_ok=True)

# #     fam1 = crawler.Family("GH1", "test_class", "test_url")

# #     config_dict = {"Glycoside Hydrolases": ["GH1"]}

# #     def mock_get_classes(*args, **kwargs):
# #         class1 = crawler.CazyClass(
# #             name="Glycoside Hydrolases",
# #             url="test_class_url.html",
# #             tries=0,
# #             failed_families={fam1: 0},
# #         )
# #         return [class1]

# #     def mock_parse_family(*args, **kwargs):
# #         return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"], ["format error"], {}

# #     monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
# #     monkeypatch.setattr(scrape_by_kingdom, "parse_family_by_kingdom", mock_parse_family)

# #     cazy_scraper.get_cazy_data(
# #         cazy_home=cazy_home_url,
# #         excluded_classes=None,
# #         config_dict=config_dict,
# #         cazy_dict=cazy_dictionary,
# #         taxonomy_filters=set(),
# #         kingdoms=["Bacteria", "Viruses"],
# #         ec_filters=[],
# #         time_stamp="timestamp",
# #         session={},
# #         args=args_get_cazy_data_stdout["args"],
# #     )
# #     file_io.make_output_directory(logs_dir, True, False)
