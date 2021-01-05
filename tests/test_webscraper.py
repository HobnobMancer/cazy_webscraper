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

from scraper import cazy_webscraper, crawler, file_io, parse, utilities


@pytest.fixture
def input_dir(test_input_dir):
    dir_path = test_input_dir / "test_inputs_webscraper"
    return dir_path


@pytest.fixture
def output_dir(test_dir):
    path = test_dir / "test_outputs"
    return path


@pytest.fixture
def get_cazy_data_args_none():
    argsdict = {
        "args": Namespace(
            data_split=None,
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
        "Glycoside Hydrolases (GHs)": ["GH1"],
        "Polysaccharide Lyases (PLs)": None,
    }
    return configuration_dict


@pytest.fixture
def mock_family():
    family = crawler.Family("GH1", "Glycoside_Hydrolases(GH)")
    return family


# test classes


def test_protein_get_protein_dict():
    """Test the Protein class __str__ and __repr__."""
    protein = crawler.Protein(
        "protein_name",
        "GH1",
        ["1.2.3.4"],
        "organism",
        {"GenBank": ["link1", "linka", "linkb"], "UniProt": ["link2"], "PDB/3D": ["link3"]},
    )

    protein_dict = protein.get_protein_dict()

    expected_protein_dict = {
        'Protein_name': ['protein_name'],
        'CAZy_family': ['GH1'],
        'EC#': ['1.2.3.4'],
        'Source_organism': ['organism'],
        'GenBank': ['link1,\nlinka,\nlinkb'],
        'UniProt': ['link2'],
        'PDB/3D': ['link3'],
    }

    assert expected_protein_dict == protein_dict


def test_family_get_name():
    """Tests get family name for Family."""
    family = crawler.Family("GH1", "Glycoside_Hydrolases(GH)")

    family_name = family.get_family_name()

    exepected_name = "GH1"

    assert exepected_name == family_name


# test main()


def test_main_one(test_dir, output_dir, null_logger, cazy_dictionary, monkeypatch):
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
            genbank="dummy_email",
            genbank_output=output_dir,
            pdb="pdb",
            pdb_output=output_dir,
            retries=0,
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


def test_main_two(output_dir, null_logger, cazy_dictionary, monkeypatch):
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
            genbank="dummy_email",
            genbank_output=output_dir,
            pdb="pdb",
            pdb_output=output_dir,
            retries=0,
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


def test_get_cazy_data_no_fam_urls(
    cazy_home_url,
    cazy_dictionary,
    config_dict,
    time_stamp,
    null_logger,
    get_cazy_data_args_none,
    monkeypatch
):
    """Test get_cazy_data."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    def mock_get_class_urls(*args, **kwargs):
        return [
            ["http://www.cazy.org/Glycoside-Hydrolases.html", 0],
            ["http://www.cazy.org/Polysaccharide-Lyases.html", 0],
        ]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(cazy_webscraper, "get_class_urls", mock_get_class_urls)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_no_return)
    monkeypatch.setattr(parse, "proteins_to_dataframe", mock_no_return)
    monkeypatch.setattr(file_io, "write_out_failed_scrapes", mock_no_return)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        config_dict,
        cazy_dict,
        1,
        time_stamp,
        null_logger,
        get_cazy_data_args_none["args"],
    )


def test_get_cazy_data_no_fam_returned(
    time_stamp,
    cazy_home_url,
    cazy_dictionary,
    config_dict,
    null_logger,
    get_cazy_data_args_none,
    monkeypatch
):
    """Test get_cazy_data() when data_split is None, subfamilies is False."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    def mock_get_class_urls(*args, **kwargs):
        return [
            ["http://www.cazy.org/Glycoside-Hydrolases.html", 0],
            ["http://www.cazy.org/Polysaccharide-Lyases.html", 0],
        ]

    def mock_get_family_urls(*args, **kwargs):
        return [
            "dddwww.cazddddy.org/GH1.html",
            "http://www.cazy.org/GH1.html",
            "http://www.cazy.org/PL1.html",
        ]

    def mock_parse_family(*args, **kwargs):
        return [None, "mock error message"]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(cazy_webscraper, "get_class_urls", mock_get_class_urls)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_family_urls)
    monkeypatch.setattr(crawler, "parse_family", mock_parse_family)
    monkeypatch.setattr(parse, "proteins_to_dataframe", mock_no_return)
    monkeypatch.setattr(file_io, "write_out_failed_scrapes", mock_no_return)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        config_dict,
        cazy_dict,
        1,
        time_stamp,
        null_logger,
        get_cazy_data_args_none["args"],
    )


def test_get_cazy_data_subfams_split_fam(
    time_stamp,
    cazy_home_url,
    cazy_dictionary,
    config_dict,
    mock_family,
    null_logger,
    get_cazy_data_args_family,
    monkeypatch
):
    """Test get_cazy_data() when data split is family and subfamily is True."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    def mock_get_class_urls(*args, **kwargs):
        return [
            ["http://www.cazy.org/Glycoside-Hydrolases.html", 0],
            ["http://www.cazy.org/Polysaccharide-Lyases.html", 0],
        ]

    def mock_get_family_urls(*args, **kwargs):
        return [
            "dddwww.cazddddy.org/GH1.html",
            "http://www.cazy.org/GH1.html",
            "http://www.cazy.org/PL1.html",
        ]

    def mock_parse_family(*args, **kwargs):
        return [mock_family, None]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(cazy_webscraper, "get_class_urls", mock_get_class_urls)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_family_urls)
    monkeypatch.setattr(crawler, "parse_family", mock_parse_family)
    monkeypatch.setattr(parse, "proteins_to_dataframe", mock_no_return)
    monkeypatch.setattr(file_io, "write_out_failed_scrapes", mock_no_return)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        config_dict,
        cazy_dict,
        1,
        time_stamp,
        null_logger,
        get_cazy_data_args_family["args"],
    )


def test_get_cazy_data_subfams_split_class(
    time_stamp,
    cazy_home_url,
    cazy_dictionary,
    config_dict,
    mock_family,
    null_logger,
    get_cazy_data_args_class,
    monkeypatch
):
    """Test get_cazy_data() when data split is class and subfamily is True."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    def mock_get_class_urls(*args, **kwargs):
        return [
            ["http://www.cazy.org/Glycoside-Hydrolases.html", 0],
            ["http://www.cazy.org/Polysaccharide-Lyases.html", 0],
        ]

    def mock_get_family_urls(*args, **kwargs):
        return [
            "dddwww.cazddddy.org/GH1.html",
            "http://www.cazy.org/GH1.html",
            "http://www.cazy.org/PL1.html",
        ]

    def mock_parse_family(*args, **kwargs):
        return [mock_family, None]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(cazy_webscraper, "get_class_urls", mock_get_class_urls)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_family_urls)
    monkeypatch.setattr(crawler, "parse_family", mock_parse_family)
    monkeypatch.setattr(parse, "proteins_to_dataframe", mock_no_return)
    monkeypatch.setattr(file_io, "write_out_failed_scrapes", mock_no_return)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        config_dict,
        cazy_dict,
        1,
        time_stamp,
        null_logger,
        get_cazy_data_args_class["args"],
    )


def test_get_cazy_data_split_none(
    time_stamp,
    cazy_home_url,
    cazy_dictionary,
    config_dict,
    mock_family,
    null_logger,
    get_cazy_data_args_none,
    monkeypatch
):
    """Test get_cazy_data() when data split is None and subfamily is False."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    def mock_get_class_urls(*args, **kwargs):
        return [
            ["http://www.cazy.org/Glycoside-Hydrolases.html", 0],
            ["http://www.cazy.org/Polysaccharide-Lyases.html", 0],
        ]

    def mock_get_family_urls(*args, **kwargs):
        return [
            "dddwww.cazddddy.org/GH1.html",
            "http://www.cazy.org/GH1.html",
            "http://www.cazy.org/PL1.html",
        ]

    def mock_parse_family(*args, **kwargs):
        return [mock_family, None]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(cazy_webscraper, "get_class_urls", mock_get_class_urls)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_family_urls)
    monkeypatch.setattr(crawler, "parse_family", mock_parse_family)
    monkeypatch.setattr(parse, "proteins_to_dataframe", mock_no_return)
    monkeypatch.setattr(file_io, "write_out_failed_scrapes", mock_no_return)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        config_dict,
        cazy_dict,
        1,
        time_stamp,
        null_logger,
        get_cazy_data_args_none["args"],
    )


# test get_class_urls


def test_get_class_urls_no_urls(cazy_home_url, null_logger, monkeypatch):
    """Test get_class_urls() when an empty list is returned."""

    def mock_class_urls(*args, **kwargs):
        return []

    monkeypatch.setattr(crawler, "get_cazy_class_urls", mock_class_urls)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy_webscraper.get_class_urls(cazy_home_url, None, 1, null_logger)
    assert pytest_wrapped_e.type == SystemExit


def test_get_class_urls_none_urls(cazy_home_url, null_logger, monkeypatch):
    """Test get_class_urls() when a None type object is returned."""

    def mock_class_urls(*args, **kwargs):
        return

    monkeypatch.setattr(crawler, "get_cazy_class_urls", mock_class_urls)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy_webscraper.get_class_urls(cazy_home_url, None, 1, null_logger)
    assert pytest_wrapped_e.type == SystemExit


def test_get_class_urls_success(cazy_home_url, null_logger, monkeypatch):
    """Test get_class_urls() when URLs were returned."""

    def mock_class_urls(*args, **kwargs):
        return ["url1", "url2", "url3"]

    monkeypatch.setattr(crawler, "get_cazy_class_urls", mock_class_urls)

    assert [['url1', 0], ['url2', 0], ['url3', 0]] == cazy_webscraper.get_class_urls(
        cazy_home_url,
        None,
        1,
        null_logger,
    )
