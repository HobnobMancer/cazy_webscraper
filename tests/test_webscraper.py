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

from bs4 import BeautifulSoup

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
            data_split="family",
            subfamilies=True,
        )
    }
    return argsdict


@pytest.fixture
def args_datasplit_class():
    argsdict = {
        "args": Namespace(
            data_split="class"
        )
    }
    return argsdict


@pytest.fixture
def config_dict():
    configuration_dict = {
        "Glycoside Hydrolases (GHs)": ["GH1"]
    }
    return configuration_dict


@pytest.fixture
def args_datasplit_family_sub_false():
    argsdict = {
        "args": Namespace(
            data_split="family",
            subfamilies=False,
        )
    }
    return argsdict


@pytest.fixture
def cazy_home_page(test_input_dir):
    file_path = test_input_dir / "test_inputs_webscraper" / "cazy_homepage.html"
    return file_path


@pytest.fixture
def cazy_class_page(test_input_dir):
    file_path = test_input_dir / "test_inputs_webscraper" / "cazy_classpage.html"
    return file_path


@pytest.fixture
def family_urls(test_input_dir):
    file_path = test_input_dir / "test_inputs_webscraper" / "family_urls.txt"
    with open(file_path, "r") as fh:
        fam_string = fh.read()
    fam_string = fam_string[1:-1]
    fam_string = fam_string.replace("'", "")
    fam_list = fam_string.split(", ")
    return fam_list


@pytest.fixture
def family_h3_element(cazy_class_page):
    with open(cazy_class_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    return [_ for _ in
            soup.find_all("h3", {"class": "spip"}) if
            str(_.contents[0]) == "Tables for Direct Access"][0]


@pytest.fixture
def subfamily_urls(test_input_dir):
    file_path = test_input_dir / "test_inputs_webscraper" / "subfamily_urls.txt"
    with open(file_path, "r") as fh:
        fam_string = fh.read()
    fam_string = fam_string[1:-1]
    fam_string = fam_string.replace("'", "")
    fam_list = fam_string.split(", ")
    return fam_list


@pytest.fixture
def no_subfam_h3_element(test_input_dir):
    file_path = test_input_dir / "test_inputs_webscraper" / "cazy_classpage_no_subfams.html"
    with open(file_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    return [_ for _ in
            soup.find_all("h3", {"class": "spip"}) if
            str(_.contents[0]) == "Tables for Direct Access"][0]


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


def test_get_cazy_data_no_config_ds_class(
    cazy_home_url,
    cazy_dictionary,
    args_datasplit_class,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_data wit no config dict and data split is 'class'"""
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
        args_datasplit_class["args"]
    )


def test_get_cazy_data_no_config_ds_class_lenfamilies_0(
    cazy_home_url,
    cazy_dictionary,
    args_datasplit_class,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_data wit no config dict and data split is class, and len(families) == 0"""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    def mock_class_pages(*args, **kwargs):
        return ["http://www.cazy.org/Glycoside-Hydrolases.html"]

    def mock_family_urls(*args, **kwargs):
        return []

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
        args_datasplit_class["args"]
    )


def test_get_cazy_data_no_config_ds_none_lenfamilies_0(
    cazy_home_url,
    cazy_dictionary,
    args_datasplit_none,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_data wit no config dict and data split is None, and len(families) == 0"""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    def mock_class_pages(*args, **kwargs):
        return ["http://www.cazy.org/Glycoside-Hydrolases.html"]

    def mock_family_urls(*args, **kwargs):
        return []

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
        args_datasplit_none["args"]
    )


def test_get_cazy_data_no_config_ds_none(
    cazy_home_url,
    cazy_dictionary,
    args_datasplit_none,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_data wit no config dict and data split is 'None'"""
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
        args_datasplit_none["args"]
    )


def test_get_cazy_data_ds_none(
    cazy_home_url,
    cazy_dictionary,
    config_dict,
    args_datasplit_family,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_data with a config dict and data split is 'family'"""
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
        config_dict,
        cazy_dict,
        null_logger,
        args_datasplit_family["args"]
    )


def test_get_cazy_data_ds_none_subfams_false(
    cazy_home_url,
    cazy_dictionary,
    config_dict,
    args_datasplit_family_sub_false,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_data with a config dict and data split is 'None' and subfamilies is False"""
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
        config_dict,
        cazy_dict,
        null_logger,
        args_datasplit_family_sub_false["args"]
    )


# test get_cazy_class_urls()


def test_get_class_urls_exclusions_none(cazy_home_url, cazy_home_page, null_logger, monkeypatch):
    """Test get_cazy_class_urls when excluded_classess is None."""
    with open(cazy_home_page, "r") as fp:
        home_page = BeautifulSoup(fp, features="lxml")

        def mock_get_home_page(*args, **kwargs):
            return [home_page, None]

        monkeypatch.setattr(cazy_webscraper, "get_page", mock_get_home_page)

        expected_result = [
            'http://www.cazy.org/Glycoside-Hydrolases.html',
            'http://www.cazy.org/GlycosylTransferases.html',
            'http://www.cazy.org/Polysaccharide-Lyases.html',
            'http://www.cazy.org/Carbohydrate-Esterases.html',
            'http://www.cazy.org/Auxiliary-Activities.html',
            'http://www.cazy.org/Carbohydrate-Binding-Modules.html',
            ]

        assert expected_result == cazy_webscraper.get_cazy_class_urls(
            cazy_home_url,
            None,
            null_logger,
        )


def test_get_class_urls_fail(cazy_home_url, null_logger, monkeypatch):
    """Test get_cazy_class_urls home_page not returned"""

    def mock_get_home_page(*args, **kwargs):
        return [None, "error"]

    monkeypatch.setattr(cazy_webscraper, "get_page", mock_get_home_page)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy_webscraper.get_cazy_class_urls(cazy_home_url, None, null_logger)
    assert pytest_wrapped_e.type == SystemExit


def test_get_class_urls_exclusions_given(
    cazy_home_url,
    cazy_home_page,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_class_urls when excluded_classess is not None."""
    with open(cazy_home_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    exclusions = ["<strong>Glycoside Hydrolases (GHs)</strong>"]

    def mock_get_home_page(*args, **kwargs):
        return [soup, None]

    monkeypatch.setattr(cazy_webscraper, "get_page", mock_get_home_page)

    expected_result = [
        'http://www.cazy.org/GlycosylTransferases.html',
        'http://www.cazy.org/Polysaccharide-Lyases.html',
        'http://www.cazy.org/Carbohydrate-Esterases.html',
        'http://www.cazy.org/Auxiliary-Activities.html',
        'http://www.cazy.org/Carbohydrate-Binding-Modules.html',
        ]

    assert expected_result == cazy_webscraper.get_cazy_class_urls(
        cazy_home_url,
        exclusions,
        null_logger,
    )


# test get_cazy_family_urls


def test_get_family_urls_fail(args_datasplit_family, null_logger, monkeypatch):
    """Test get_cazy_family_urls when no page is returned."""

    def mock_get_page(*args, **kwargs):
        return [None, "error"]

    monkeypatch.setattr(cazy_webscraper, "get_page", mock_get_page)

    assert None is cazy_webscraper.get_cazy_family_urls(
        "class_url",
        "cazy_home_url",
        "class_name",
        args_datasplit_family["args"],
        null_logger,
    )


def test_get_family_urls_success(
    cazy_class_page,
    args_datasplit_family,
    family_urls,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_family_urls when successful, and subfamilies is True."""
    with open(cazy_class_page) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return [page, None]

    def mock_get_subfams(*args, **kwargs):
        return []

    monkeypatch.setattr(cazy_webscraper, "get_page", mock_get_page)
    monkeypatch.setattr(cazy_webscraper, "get_subfamily_links", mock_get_subfams)

    assert family_urls == cazy_webscraper.get_cazy_family_urls(
        "class_url",
        "http://www.cazy.org",
        "class_name",
        args_datasplit_family["args"],
        null_logger,
    )


# test get_subfamily_links


def test_get_subfam_links_len_0(family_h3_element, subfamily_urls, null_logger):
    """Test get_subfamily_links when no links are retrieved."""

    assert subfamily_urls == cazy_webscraper.get_subfamily_links(
        family_h3_element,
        "http://www.cazy.org",
        null_logger,
    )


def test_get_subfam_links_urls(no_subfam_h3_element, null_logger):
    """Test get_subfamily_links when urls are retrieved."""

    assert None is cazy_webscraper.get_subfamily_links(
        no_subfam_h3_element,
        "http://www.cazy.org",
        null_logger,
    )


# test parse_family()


def test_parse_family(null_logger, monkeypatch):
    """Tests parse_family()"""

    def mock_parsing_family_pages(*args, **kwargs):
        protein_list = [
            cazy_webscraper.Protein("protein_name", "GH1", "1.2.3.4", "organism"),
            cazy_webscraper.Protein(
                "protein",
                "GH1",
                "",
                "organism",
                {"GenBank": ["link1"], "UniProt": ["link2"], "PDB": ["link3"]},
            ),
        ]
        return protein_list

    monkeypatch.setattr(cazy_webscraper, "parse_family_pages", mock_parsing_family_pages)

    cazy_webscraper.parse_family(
        "http://www.cazy.org/GH1.html",
        "GH1",
        "http://www.cazy.org",
        null_logger,
    )
