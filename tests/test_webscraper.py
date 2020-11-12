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

from scraper import cazy_webscraper, file_io, parse, utilities


@pytest.fixture
def cazy_dictionary(test_input_dir):
    dict_path = test_input_dir / "test_inputs_webscraper" / "cazy_dictionary.json"
    return dict_path


@pytest.fixture
def args_datasplit_none():
    argsdict = {
        "args": Namespace(
            data_split=None
        )
    }
    return argsdict


@pytest.fixture
def args_datasplit_family():
    argsdict = {
        "args": Namespace(
            data_split="family"
        )
    }
    return argsdict


# test main()


def test_main_one(test_dir, null_logger, cazy_dictionary, monkeypatch):
    """Test function main().

    Argv is None, logger is None, args.output is not sys.stdout, args.subfamilies is True.
    """
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

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
        return None, None, cazy_dict

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
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

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
        return None, None, cazy_dict

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


def test_get_cazy_data_class_urls_empty_list(cazy_home_url, null_logger, monkeypatch):
    """Test get_cazy_data, checking that when no class page urls are returned the program exits."""

    def mock_failed_class_pages(*args, **kwargs):
        return []

    monkeypatch.setattr(cazy_webscraper, "get_cazy_class_urls", mock_failed_class_pages)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy_webscraper.get_cazy_data(cazy_home_url, None, None, None, null_logger, None)
    assert pytest_wrapped_e.type == SystemExit


def test_get_cazy_data_family_urls_none(
    cazy_home_url,
    cazy_dictionary,
    args_datasplit_none,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_data when nothing returned from get_cazy_family_urls."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    def mock_class_pages(*args, **kwargs):
        return ["http://www.cazy.org/Glycoside-Hydrolases.html"]

    def mock_family_urls_none(*args, **kwargs):
        return

    monkeypatch.setattr(cazy_webscraper, "get_cazy_class_urls", mock_class_pages)
    monkeypatch.setattr(cazy_webscraper, "get_cazy_family_urls", mock_family_urls_none)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        None,
        cazy_dict,
        null_logger,
        args_datasplit_none["args"]
    )


def test_get_cazy_data_no_config_ds_family(
    cazy_home_url,
    cazy_dictionary,
    args_datasplit_family,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_data wit no config dict and data split is 'family'"""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    def mock_class_pages(*args, **kwargs):
        return ["http://www.cazy.org/Glycoside-Hydrolases.html"]

    def mock_family_urls(*args, **kwargs):
        return ["http://www.cazy.org/GH1.html"]

    def mock_parsing_family(*args, **kwargs):
        family = Namespace(
            name="cazy fam",
            cazy_class="cazy_class",
            members=set([1, 2, 3])
        )
        return family

    def mock_building_dataframe(*args, **kwargs):
        return

    monkeypatch.setattr(cazy_webscraper, "get_cazy_class_urls", mock_class_pages)
    monkeypatch.setattr(cazy_webscraper, "get_cazy_family_urls", mock_family_urls)
    monkeypatch.setattr(cazy_webscraper, "parse_family", mock_parsing_family)
    monkeypatch.setattr(parse, "proteins_to_dataframe", mock_building_dataframe)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        None,
        cazy_dict,
        null_logger,
        args_datasplit_family["args"]
    )
