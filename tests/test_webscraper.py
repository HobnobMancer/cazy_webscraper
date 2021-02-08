#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License

"""Tests the script cazy_webscraper.py which coordinates the scraping of CAZy.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import json
import pytest
import sys

from argparse import Namespace, ArgumentParser
from requests.exceptions import MissingSchema

from scraper import cazy_webscraper, crawler, file_io, utilities


@pytest.fixture
def input_dir(test_input_dir):
    dir_path = test_input_dir / "test_inputs_webscraper"
    return dir_path


@pytest.fixture
def output_dir(test_dir):
    path = test_dir / "test_outputs"
    return path


@pytest.fixture
def db_path(test_input_dir):
    db_path = test_input_dir / "test_db" / "test_cazy.db"
    return db_path


@pytest.fixture
def args_get_cazy_data():
    argsdict = {
        "args": Namespace(
            subfamilies=False,
        )
    }
    return argsdict


@pytest.fixture
def get_cazy_data_args_family():
    argsdict = {
        "args": Namespace(
            data_split="family",
            subfamilies=True,
        )
    }
    return argsdict


@pytest.fixture
def get_cazy_data_args_class():
    argsdict = {
        "args": Namespace(
            data_split="class",
            subfamilies=True,
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


@pytest.fixture
def mock_family():
    family = crawler.Family("GH1", "Glycoside_Hydrolases(GH)")
    return family


# test main()


def test_main_invalid_db_path(output_dir, null_logger, cazy_dictionary, monkeypatch):
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
            output=output_dir,
            subfamilies=True,
            force=False,
            nodelete=False,
            retries=1,
            database="fake_database_path",
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        return None, None, cazy_dictionary

    monkeypatch.setattr(utilities, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(file_io, "parse_configuration", mock_retrieving_configuration)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy_webscraper.main()
    assert pytest_wrapped_e.type == SystemExit


def test_main_existing_database(output_dir, null_logger, cazy_dictionary, db_path, monkeypatch):
    """Test function main() when passed an existing database.

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
            output=output_dir,
            subfamilies=True,
            force=False,
            nodelete=False,
            retries=1,
            database=db_path,
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        return None, None, cazy_dictionary

    def mock_retrieving_cazy_data(*args, **kwargs):
        return

    monkeypatch.setattr(utilities, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(file_io, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "get_cazy_data", mock_retrieving_cazy_data)

    cazy_webscraper.main(["argv"])


def test_main_new_database(output_dir, null_logger, cazy_dictionary, db_path, monkeypatch):
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

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            output=output_dir,
            subfamilies=True,
            force=False,
            nodelete=False,
            retries=1,
            database=None,
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        return None, None, cazy_dictionary

    def mock_retrieving_cazy_data(*args, **kwargs):
        return

    monkeypatch.setattr(utilities, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(file_io, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "get_cazy_data", mock_retrieving_cazy_data)

    cazy_webscraper.main()


# test get_cazy_data()


def test_get_cazy_data_no_fam_urls(
    cazy_home_url,
    cazy_dictionary,
    config_dict,
    time_stamp,
    args_get_cazy_data,
    monkeypatch
):
    """Test get_cazy_data() when no family URLS are retrieved, and fails to scrape families from
    a class.
    """

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass("test_class", "test_class_url.html", 0)
        class2 = crawler.CazyClass("test_class2", "test_class_url2.html", 0)
        return [class1, class2]

    def mock_get_families(*args, **kwargs):
        return None, "test error message", ["test_url1", "test_url2"]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_families)
    monkeypatch.setattr(file_io, "write_out_failed_scrapes", mock_no_return)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        config_dict,
        cazy_dictionary,
        2,
        time_stamp,
        "session",
        args_get_cazy_data["args"],
    )


def test_get_cazy_data_no_config_dict(
    time_stamp,
    cazy_home_url,
    cazy_dictionary,
    args_get_cazy_data,
    monkeypatch
):
    """Test get_cazy_data() when some families aren't scraped, and config_dict is None."""

    fam1 = crawler.Family("test_fam", "test_class", "test_url")

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass("test_class", "test_class_url.html", 0, {fam1: 0})
        return [class1]

    def mock_parse_family(*args, **kwargs):
        return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(crawler, "parse_family", mock_parse_family)
    monkeypatch.setattr(file_io, "write_out_failed_proteins", mock_no_return)
    monkeypatch.setattr(file_io, "write_out_failed_scrapes", mock_no_return)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        None,
        cazy_dictionary,
        2,
        time_stamp,
        "session",
        args_get_cazy_data["args"],
    )


def test_get_cazy_data_with_config_dict(
    time_stamp,
    cazy_home_url,
    cazy_dictionary,
    config_dict,
    args_get_cazy_data,
    monkeypatch
):
    """Test get_cazy_data() when some families aren't scraped, and a config_dict is given."""

    fam1 = crawler.Family("GH3_1", "Glycoside Hydrolases (GHs)", "test_url")

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass("Glycoside Hydrolases (GHs)", "test_class_url.html", 0, {fam1: 0})
        return [class1]

    def mock_parse_family(*args, **kwargs):
        return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(crawler, "parse_family", mock_parse_family)
    monkeypatch.setattr(file_io, "write_out_failed_proteins", mock_no_return)
    monkeypatch.setattr(file_io, "write_out_failed_scrapes", mock_no_return)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        None,
        cazy_dictionary,
        2,
        time_stamp,
        "session",
        args_get_cazy_data["args"],
    )
