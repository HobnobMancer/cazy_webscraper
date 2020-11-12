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

"""Tests the script cazy_webscraper.py which crawls the CAZy database and retrieves protein data.

These test are intened to be run from the root of the repository using:
pytest -v

"""

import json
import pytest
import sys

from argparse import Namespace, ArgumentParser

from scraper import cazy_webscraper, file_io, utilities


@pytest.fixture
def cazy_dictionary(test_input_dir):
    dict_path = test_input_dir / "test_inputs_crawler" / "cazy_dictionary.json"
    with open(dict_path, "r") as fh:
        cazy_dict = json.load(fh)
    return cazy_dict


@pytest.fixture
def args_datasplit_none():
    argsdict = {
        "args": Namespace(
            data_split=None
        )
    }


# test main()


def test_main_one(test_dir, null_logger, cazy_dictionary, monkeypatch):
    """Test function main().

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
            output=test_dir,
            subfamilies=True,
            force=False,
            nodelete=False,
        )
        return parser

    def mock_build_logger(*args, **kwargs):
        return null_logger

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        return None, None, cazy_dictionary

    def mock_retrieving_cazy_data(*args, **kwargs):
        return

    monkeypatch.setattr(utilities, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "build_logger", mock_build_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(file_io, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "get_cazy_data", mock_retrieving_cazy_data)

    cazy_webscraper.main()


def test_main_two(null_logger, cazy_dictionary, monkeypatch):
    """Test function main() using opposite scenariors for if statements than in test_main_one().

    Argv is not None, logger is not None, args.output is sys.stdout, args.subfamilies is False.
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
            output=sys.stdout,
            subfamilies=False,
            force=False,
            nodelete=False,
        )
        return parser

    def mock_build_logger(*args, **kwargs):
        return null_logger

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        return None, None, cazy_dictionary

    def mock_retrieving_cazy_data(*args, **kwargs):
        return

    monkeypatch.setattr(utilities, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "build_logger", mock_build_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(file_io, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "get_cazy_data", mock_retrieving_cazy_data)

    cazy_webscraper.main(["argv"], null_logger)


# test get_cazy_data


def test_get_cazy_data_class_urls_none(cazy_home_url, null_logger, monkeypatch):
    """Test get_cazy_data, checking that when no class page urls are returned the program exits."""

    def mock_failed_class_pages(*args, **kwargs):
        return

    monkeypatch.setattr(cazy_webscraper, "get_cazy_class_urls", mock_failed_class_pages)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy_webscraper.get_cazy_data(cazy_home_url, None, None, None, null_logger, None)
    assert pytest_wrapped_e.type == SystemExit


# def test_get_cazy_family_urls_none(cazy_home_url, args_datasplit_none, null_logger, monkeypatch):
#     """Test get_cazy_data when nothing returned from get_cazy_family_urls."""

#     def mock_class_pages(*args, **kwargs):
#         return ["http://www.cazy.org/Glycoside-Hydrolases.html"]
    
#     def mock_family_urls_none(*args, **kwargs):
#         return
    
#     monkeypatch.setattr(cazy_webscraper, "get")
