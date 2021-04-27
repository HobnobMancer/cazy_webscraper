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
import pytest
import sys

from argparse import Namespace, ArgumentParser

from scraper import cazy_webscraper, crawler, sql, utilities
from scraper.crawler.cazy_html_pages import get_cazy_pages, parse_local_pages
from scraper.crawler.parse_cazy_families import scrape_all, scrape_by_kingdom
from scraper.utilities import file_io, parse_configuration, parsers


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
    db_path = "tests/test_inputs/test_inputs_sql/unit_test_db_2021-03-01--15-06-59.db"
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


def test_main_get_pages(output_dir, cazy_dictionary, taxonomic_filter_dict, monkeypatch):
    """Test function main() when retrieval of CAZy HTML pages is enabled.

    Argv is None, logger is None, args.output is not sys.stdout, args.subfamilies is True.
    """

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="cazy_webscraper.py",
            usage=None,
            description="Scrape the CAZy database",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

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

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        # excluded_classes, config_dict, cazy_dict, taxonomy_filter, kingdoms, ec_filter
        return None, None, cazy_dictionary, taxonomic_filter_dict, [], []

    def mock_filter_set(*args, **kwargs):
        return set()

    def mock_get_pages(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "get_filter_set", mock_filter_set)
    monkeypatch.setattr(get_cazy_pages, "get_cazy_pages", mock_get_pages)

    cazy_webscraper.main()


def test_main_invalid_db_path(output_dir, cazy_dictionary, taxonomic_filter_dict, monkeypatch):
    """Test function main() when an invalid db path is given.

    Argv is None, logger is None, args.output is not sys.stdout, args.subfamilies is True.
    """

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="cazy_webscraper.py",
            usage=None,
            description="Scrape the CAZy database",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            config=None,
            classes=None,
            database="fake_database_path",
            ec=None,
            force=False,
            families=None,
            genera=None,
            get_pages=False,
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

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        # excluded_classes, config_dict, cazy_dict, taxonomy_filter, kingdoms, ec_filter
        return None, None, cazy_dictionary, taxonomic_filter_dict, [], []

    def mock_filter_set(*args, **kwargs):
        return set()

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "get_filter_set", mock_filter_set)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy_webscraper.main()
    assert pytest_wrapped_e.type == SystemExit


def test_main_db_raises_error(output_dir, cazy_dictionary, taxonomic_filter_dict, monkeypatch):
    """Test function main() when aan error is raised when building the local db."""

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="cazy_webscraper.py",
            usage=None,
            description="Scrape the CAZy database",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            config=None,
            classes=None,
            database=None,
            ec=None,
            force=False,
            families=None,
            genera=None,
            get_pages=False,
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

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        # excluded_classes, config_dict, cazy_dict, taxonomy_filter, kingdoms, ec_filter
        return None, None, cazy_dictionary, taxonomic_filter_dict, [], []

    def mock_filter_set(*args, **kwargs):
        return set()

    def mock_db_build(*args, **kwargs):
        raise TypeError

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "get_filter_set", mock_filter_set)
    monkeypatch.setattr(sql.sql_orm, "build_db", mock_db_build)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy_webscraper.main()
    assert pytest_wrapped_e.type == SystemExit


def test_main_existing_db_scrape_local_pages(
    output_dir,
    null_logger,
    cazy_dictionary,
    db_path,
    input_dir,
    monkeypatch,
):
    """Test function main() when passed an existing database scraping data from local HTML pages.

    Argv is not None, logger is not None, args.output is output_dir, args.subfamilies is True,
    and valid db path is given by db_path.
    """

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="cazy_webscraper.py",
            usage=None,
            description="Scrape the CAZy database",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            config=None,
            classes=None,
            database=db_path,
            ec=None,
            force=False,
            families=None,
            genera=None,
            get_pages=False,
            kingdoms=None,
            log=None,
            nodelete=False,
            output=output_dir,
            retries=10,
            scrape_files=input_dir,
            subfamilies=True,
            species=None,
            strains=None,
            streamline="streamline_args",
            timeout=45,
            verbose=False,
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        # excluded_classes, config_dict, cazy_dict, taxonomy_filter, kingdoms, ec_filter
        return None, None, cazy_dictionary, taxonomic_filter_dict, [], []

    def mock_filter_set(*args, **kwargs):
        return set()

    def mock_none(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "get_filter_set", mock_filter_set)
    monkeypatch.setattr(parse_configuration, "create_streamline_scraping_warning", mock_none)
    monkeypatch.setattr(sql.sql_interface, "log_scrape_in_db", mock_none)
    monkeypatch.setattr(parse_local_pages, "parse_local_pages", mock_none)

    cazy_webscraper.main(["argv"])


def test_main_existing_db_scrape_direct(
    output_dir,
    null_logger,
    cazy_dictionary,
    db_path,
    input_dir,
    monkeypatch,
):
    """Test function main() when passed an existing database scraping CAZy directly.

    Argv is not None, logger is not None, args.output is output_dir, args.subfamilies is True,
    and valid db path is given by db_path.
    """

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="cazy_webscraper.py",
            usage=None,
            description="Scrape the CAZy database",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            config=None,
            classes=None,
            database=db_path,
            ec=None,
            force=False,
            families=None,
            genera=None,
            get_pages=False,
            kingdoms=None,
            log=None,
            nodelete=False,
            output=output_dir,
            retries=10,
            scrape_files=None,
            subfamilies=True,
            species=None,
            strains=None,
            streamline="streamline_args",
            timeout=45,
            verbose=False,
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        # excluded_classes, config_dict, cazy_dict, taxonomy_filter, kingdoms, ec_filter
        return None, None, cazy_dictionary, taxonomic_filter_dict, [], []

    def mock_filter_set(*args, **kwargs):
        return set()

    def mock_none(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "get_filter_set", mock_filter_set)
    monkeypatch.setattr(parse_configuration, "create_streamline_scraping_warning", mock_none)
    monkeypatch.setattr(sql.sql_interface, "log_scrape_in_db", mock_none)
    monkeypatch.setattr(cazy_webscraper, "get_cazy_data", mock_none)

    cazy_webscraper.main(["argv"])


def test_main_new_database(output_dir, cazy_dictionary, monkeypatch):
    """Test main() when a new database_file is created"""

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="cazy_webscraper.py",
            usage=None,
            description="Scrape the CAZy database",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    output_path = output_dir / "test_webscraper" / "temp_dir_for_db"
    os.makedirs(output_path, exist_ok=True)

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            config=None,
            classes=None,
            database=None,
            ec=None,
            force=True,
            families=None,
            genera=None,
            get_pages=False,
            kingdoms=None,
            log=None,
            nodelete=False,
            output=output_path,
            retries=10,
            scrape_files=None,
            subfamilies=True,
            species=None,
            strains=None,
            streamline="streamline_args",
            timeout=45,
            verbose=False,
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        # excluded_classes, config_dict, cazy_dict, taxonomy_filter, kingdoms, ec_filter
        return None, None, cazy_dictionary, taxonomic_filter_dict, [], []

    def mock_filter_set(*args, **kwargs):
        return set()

    def mock_none(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "get_filter_set", mock_filter_set)
    monkeypatch.setattr(parse_configuration, "create_streamline_scraping_warning", mock_none)
    monkeypatch.setattr(sql.sql_interface, "log_scrape_in_db", mock_none)
    monkeypatch.setattr(cazy_webscraper, "get_cazy_data", mock_none)

    cazy_webscraper.main(["argv"])

    # delete newly build db
    file_io.make_output_directory(output_dir, True, False)


# test get_filter_set


def test_get_filter_set():
    """Test get_filter_set."""
    taxonomy_filter = {"genera": ["Aspergillus", "Trichoderma"], "strains": []}
    cazy_webscraper.get_filter_set(taxonomy_filter)


def test_get_filter_set_none():
    """Test get_filter_set when no filters are provided."""
    taxonomy_filter = {}
    assert cazy_webscraper.get_filter_set(taxonomy_filter) is None


# test get_cazy_data()


def test_get_cazy_data_no_fam_urls(
    cazy_home_url,
    cazy_dictionary,
    config_dict,
    time_stamp,
    args_get_cazy_data,
    logs_dir,
    monkeypatch,
    null_logger
):
    """Test get_cazy_data() when no family URLS are retrieved."""

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass("test_class", "test_class_url.html", 0)
        return [class1]

    def mock_get_families(*args, **kwargs):
        return None, "test error message", ["test_url1", "test_url2"]

    def mock_logger(*args, **kwargs):
        return null_logger

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_families)
    monkeypatch.setattr(utilities, "build_logger", mock_logger)

    cazy_webscraper.get_cazy_data(
        cazy_home=cazy_home_url,
        excluded_classes=None,
        config_dict=config_dict,
        cazy_dict=cazy_dictionary,
        taxonomy_filters=set(),
        kingdoms="all",
        ec_filters=[],
        time_stamp="timestamp",
        session="session_representative",
        args=args_get_cazy_data["args"],
    )


def test_get_cazy_data_no_all(
    time_stamp,
    cazy_home_url,
    cazy_dictionary,
    args_get_cazy_data,
    logs_dir,
    monkeypatch,
):
    """Test get_cazy_data() when no kingdoms are specified and config_dict is None."""
    # prepare dir for log files
    os.makedirs(logs_dir, exist_ok=True)

    fam1 = crawler.Family("test_fam", "test_class", "test_url")

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass(
            name="test_class",
            url="test_class_url.html",
            tries=0,
        )
        return [class1]

    def mock_get_fam_urls(*args, **kwargs):
        return [fam1], "error message", ["in", "cor", "rect", "urls"]

    def mock_parse_family(*args, **kwargs):
        return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"], ["format error"], "session"

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_fam_urls)
    monkeypatch.setattr(scrape_all, "parse_family_via_all_pages", mock_parse_family)

    cazy_webscraper.get_cazy_data(
        cazy_home=cazy_home_url,
        excluded_classes=None,
        config_dict=None,
        cazy_dict=None,
        taxonomy_filters=set(),
        kingdoms="all",
        ec_filters=[],
        time_stamp="timestamp",
        session="session_representative",
        args=args_get_cazy_data["args"],
    )
    file_io.make_output_directory(logs_dir, True, False)


def test_get_cazy_data_no_config_dict_kingdom(
    time_stamp,
    cazy_home_url,
    cazy_dictionary,
    args_get_cazy_data,
    logs_dir,
    monkeypatch,
):
    """Test get_cazy_data() when kingdoms are specified and config_dict is None."""
    # prepare dir for log files
    os.makedirs(logs_dir, exist_ok=True)

    fam1 = crawler.Family("test_fam", "test_class", "test_url")

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass(
            name="test_class",
            url="test_class_url.html",
            tries=0,
            failed_families={fam1: 0},
        )
        return [class1]

    def mock_parse_family(*args, **kwargs):
        return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"], ["format error"], "session"

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(scrape_by_kingdom, "parse_family_by_kingdom", mock_parse_family)

    cazy_webscraper.get_cazy_data(
        cazy_home=cazy_home_url,
        excluded_classes=None,
        config_dict=None,
        cazy_dict=None,
        taxonomy_filters=set(),
        kingdoms=["Bacteria", "Viruses"],
        ec_filters=[],
        time_stamp="timestamp",
        session="session_representative",
        args=args_get_cazy_data["args"],
    )
    file_io.make_output_directory(logs_dir, True, False)


def test_get_cazy_data_config_data_all(
    time_stamp,
    cazy_home_url,
    cazy_dictionary,
    args_get_cazy_data,
    logs_dir,
    monkeypatch,
):
    """Test get_cazy_data() when no kingdoms are specified and configuration given."""
    # prepare dir for log files
    os.makedirs(logs_dir, exist_ok=True)

    fam1 = crawler.Family("GH3_1", "test_class", "test_url")

    config_dict = {"Glycoside Hydrolases": ["GH3"]}

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass(
            name="Glycoside Hydrolases",
            url="test_class_url.html",
            tries=0,
            failed_families={fam1: 0},
        )
        return [class1]

    def mock_parse_family(*args, **kwargs):
        return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"], ["format error"], "session"

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(scrape_all, "parse_family_via_all_pages", mock_parse_family)

    cazy_webscraper.get_cazy_data(
        cazy_home=cazy_home_url,
        excluded_classes=None,
        config_dict=config_dict,
        cazy_dict=cazy_dictionary,
        taxonomy_filters=set(),
        kingdoms="all",
        ec_filters=[],
        time_stamp="timestamp",
        session="session_representative",
        args=args_get_cazy_data["args"],
    )
    file_io.make_output_directory(logs_dir, True, False)


def test_get_cazy_data_config_data_kingdom(
    time_stamp,
    cazy_home_url,
    cazy_dictionary,
    args_get_cazy_data,
    logs_dir,
    monkeypatch,
):
    """Test get_cazy_data() when kingdoms are specified and configuration given."""
    # prepare dir for log files
    os.makedirs(logs_dir, exist_ok=True)

    fam1 = crawler.Family("GH1", "test_class", "test_url")

    config_dict = {"Glycoside Hydrolases": ["GH1"]}

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass(
            name="Glycoside Hydrolases",
            url="test_class_url.html",
            tries=0,
            failed_families={fam1: 0},
        )
        return [class1]

    def mock_parse_family(*args, **kwargs):
        return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"], ["format error"], {}

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(scrape_by_kingdom, "parse_family_by_kingdom", mock_parse_family)

    cazy_webscraper.get_cazy_data(
        cazy_home=cazy_home_url,
        excluded_classes=None,
        config_dict=config_dict,
        cazy_dict=cazy_dictionary,
        taxonomy_filters=set(),
        kingdoms=["Bacteria", "Viruses"],
        ec_filters=[],
        time_stamp="timestamp",
        session={},
        args=args_get_cazy_data["args"],
    )
    file_io.make_output_directory(logs_dir, True, False)


def test_get_cazy_data_config_data_kingdom_stdout(
    time_stamp,
    cazy_home_url,
    cazy_dictionary,
    args_get_cazy_data_stdout,
    logs_dir,
    monkeypatch,
):
    """Test get_cazy_data() when kingdoms are specified and configuration given."""
    # prepare dir for log files
    os.makedirs(logs_dir, exist_ok=True)

    fam1 = crawler.Family("GH1", "test_class", "test_url")

    config_dict = {"Glycoside Hydrolases": ["GH1"]}

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass(
            name="Glycoside Hydrolases",
            url="test_class_url.html",
            tries=0,
            failed_families={fam1: 0},
        )
        return [class1]

    def mock_parse_family(*args, **kwargs):
        return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"], ["format error"], {}

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(scrape_by_kingdom, "parse_family_by_kingdom", mock_parse_family)

    cazy_webscraper.get_cazy_data(
        cazy_home=cazy_home_url,
        excluded_classes=None,
        config_dict=config_dict,
        cazy_dict=cazy_dictionary,
        taxonomy_filters=set(),
        kingdoms=["Bacteria", "Viruses"],
        ec_filters=[],
        time_stamp="timestamp",
        session={},
        args=args_get_cazy_data_stdout["args"],
    )
    file_io.make_output_directory(logs_dir, True, False)
