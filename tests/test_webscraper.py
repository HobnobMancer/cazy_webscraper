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

import pytest

from argparse import Namespace, ArgumentParser

from scraper import cazy_webscraper, crawler, sql, utilities
from scraper.utilities import file_io, parse_configuration, parsers


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
def args_get_cazy_data():
    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
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
            verbose=False,
            log=None,
            streamline=None,
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        return None, None, cazy_dictionary, [], []

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_retrieving_configuration)

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
            verbose=True,
            log=None,
            streamline=None,
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        return None, None, cazy_dictionary, None, []

    def mock_adding_log(*args, **kwargs):
        return

    def mock_get_filter_set(*args, **kwargs):
        return set()

    def mock_retrieving_cazy_data(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "log_scrape_in_db", mock_adding_log)
    monkeypatch.setattr(cazy_webscraper, "get_filter_set", mock_get_filter_set)
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
            verbose=False,
            log=None,
            streamline=None,
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        return None, None, cazy_dictionary, None, []

    def mock_getting_db_session(*args, **kwargs):
        return "session"

    def mock_adding_log(*args, **kwargs):
        return

    def mock_tax_set(*args, **kwargs):
        return set()

    def mock_retrieving_cazy_data(*args, **kwargs):
        return

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(sql.sql_orm, "build_db", mock_getting_db_session)
    monkeypatch.setattr(cazy_webscraper, "log_scrape_in_db", mock_adding_log)
    monkeypatch.setattr(cazy_webscraper, "get_filter_set", mock_tax_set)
    monkeypatch.setattr(cazy_webscraper, "get_cazy_data", mock_retrieving_cazy_data)

    cazy_webscraper.main()


def test_main_build_sql_error(output_dir, null_logger, cazy_dictionary, db_path, monkeypatch):
    """Test main() when an error is raised when buildin a new database."""

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
            output=1,
            subfamilies=True,
            force=False,
            nodelete=False,
            retries=1,
            database=None,
            verbose=False,
            log=None,
            streamline=None,
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        return None, None, cazy_dictionary, [], None

    def mock_adding_log(*args, **kwargs):
        return

    def mock_retrieving_cazy_data(*args, **kwargs):
        return

    def mock_building_db(*args, **kwargs):
        raise TypeError

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(sql.sql_orm, "get_db_session", mock_building_db)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "log_scrape_in_db", mock_adding_log)
    monkeypatch.setattr(cazy_webscraper, "get_cazy_data", mock_retrieving_cazy_data)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy_webscraper.main()
    assert pytest_wrapped_e.type == SystemExit


def test_main_get_sql_error(output_dir, null_logger, cazy_dictionary, db_path, monkeypatch):
    """Test main() when an error is raised when buildin a new database."""

    def mock_building_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="cazy_webscraper.py",
            usage=None,
            description="Scrape the CAZy database",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    path_ = output_dir / "test_outputs_sql" / "non_sql_file.html"

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            output=output_dir,
            subfamilies=True,
            force=False,
            nodelete=False,
            retries=1,
            database=path_,
            verbose=True,
            log=None,
            streamline='genbank,pdb,uniprot,ec'
        )
        return parser

    def mock_config_logger(*args, **kwargs):
        return

    def mock_making_output_dir(*args, **kwargs):
        return

    def mock_retrieving_configuration(*args, **kwargs):
        return None, None, cazy_dictionary, [], []

    def mock_retrieving_cazy_data(*args, **kwargs):
        return

    def mock_getting_db_session(*args, **kwargs):
        raise TypeError

    monkeypatch.setattr(parsers, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_config_logger)
    monkeypatch.setattr(file_io, "make_output_directory", mock_making_output_dir)
    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_retrieving_configuration)
    monkeypatch.setattr(cazy_webscraper, "get_cazy_data", mock_retrieving_cazy_data)
    monkeypatch.setattr(sql.sql_orm, "get_db_session", mock_getting_db_session)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy_webscraper.main()
    assert pytest_wrapped_e.type == SystemExit


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

    fam1 = crawler.Family("test_fam", "test_class", "test_url")

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass("test_class", "test_class_url.html", 0)
        class2 = crawler.CazyClass("test_class2", "test_class_url2.html", 0)
        return [class1, class2]

    def mock_get_families(*args, **kwargs):
        return None, "test error message", ["test_url1", "test_url2"]

    def mock_parse_family(*args, **kwargs):
        return fam1, True, ["fail1", "fail2"], ["sqlFail1", "sqlFail2"]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_families)
    monkeypatch.setattr(crawler, "parse_family", mock_parse_family)
    monkeypatch.setattr(file_io, "write_out_failed_scrapes", mock_no_return)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        config_dict,
        cazy_dictionary,
        None,
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
        None,
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

    fam1 = crawler.Family("GH3", "Glycoside Hydrolases (GHs)", "test_url")

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass(
            "Glycoside Hydrolases (GHs)",
            "test_class_url.html",
            0,
        )
        return [class1]

    def mock_parse_family(*args, **kwargs):
        return fam1, False, [], ["sqlFail1", "sqlFail2"]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(crawler, "parse_family", mock_parse_family)
    monkeypatch.setattr(file_io, "write_out_failed_proteins", mock_no_return)
    monkeypatch.setattr(file_io, "write_out_failed_scrapes", mock_no_return)

    cazy_webscraper.get_cazy_data(
        cazy_home_url,
        None,
        config_dict,
        cazy_dictionary,
        None,
        2,
        time_stamp,
        "session",
        args_get_cazy_data["args"],
    )


def test_get_cazy_data_witha_config_dict_subfam(
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
        class1 = crawler.CazyClass(
            "Glycoside Hydrolases (GHs)",
            "test_class_url.html",
            0,
            {fam1: 0},
        )
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
        config_dict,
        cazy_dictionary,
        None,
        2,
        time_stamp,
        "session",
        args_get_cazy_data["args"],
    )


# tests for adding log to the local database


def test_add_db_log_no_config(db_session):
    config_dict = None
    taxonomy_filters = {"genera": None, "species": None, "strains": None}
    kingdoms = ["Archaea", "Bacteria", "Eukaryota", "Viruses", "Unclassified"]
    args = {
        "args": Namespace(
            classes=None,
            families=None,
            genera=None,
            species=None,
            strains=None,
            kingdoms=None,
            streamline=None,
        )
    }

    cazy_webscraper.log_scrape_in_db(
        "YYYY-MM-DD--HH-MM-SS",
        config_dict,
        taxonomy_filters,
        kingdoms,
        db_session,
        args["args"],
    )


def test_add_db_log_with_config(db_session):
    config_dict = {
        'classes': None,
        "Polysaccharide Lyases (PLs)": ["PL2", "PL3"],
        "Glycoside Hydrolases (GHs)": None,
    }
    taxonomy_filters = {
        "genera": ["Caldivirga", "Cuniculiplasma"],
        "species": ["Pyrococcus furiosus"],
        "strains": ["Saccharolobus solfataricus POZ149", "Saccharolobus solfataricus SULB"]
    }
    kingdoms = ["Archaea"]
    args = {
        "args": Namespace(
            classes="GH,PL",
            families="AA1,AA2",
            genera="Trichoderma",
            species="Aspergillus Niger",
            strains="Acidianus ambivalens LEI 10",
            kingdoms="Archaea,Bacteria",
            streamline="uniprot,pdb",
        )
    }

    cazy_webscraper.log_scrape_in_db(
        "YYYY-MM-DD--HH-MM-SS",
        config_dict,
        taxonomy_filters,
        kingdoms,
        db_session,
        args["args"],
    )


def test_add_db_log_all_kingdoms(db_session):
    """Test adding log to database when all Kingdoms are scraped."""
    config_dict = {
        'classes': None,
        "Polysaccharide Lyases (PLs)": ["PL2", "PL3"],
        "Glycoside Hydrolases (GHs)": None,
    }
    taxonomy_filters = {
        "genera": ["Caldivirga", "Cuniculiplasma"],
        "species": ["Pyrococcus furiosus"],
        "strains": ["Saccharolobus solfataricus POZ149", "Saccharolobus solfataricus SULB"]
    }
    kingdoms = ["Archaea"]
    args = {
        "args": Namespace(
            classes="GH,PL",
            families="AA1,AA2",
            genera="Trichoderma",
            species="Aspergillus Niger",
            strains="Acidianus ambivalens LEI 10",
            kingdoms=None,
            streamline="uniprot,pdb",
        )
    }

    cazy_webscraper.log_scrape_in_db(
        "YYYY-MM-DD--HH-MM-SS",
        config_dict,
        taxonomy_filters,
        kingdoms,
        db_session,
        args["args"],
    )


# tests for parsing taxonomy filters


def test_get_filter_set():
    taxonomy_dict = {"genera": None, "species": ["species1", "species2"]}

    assert set(["species1", "species2"]) == cazy_webscraper.get_filter_set(taxonomy_dict)


def test_get_filter_set_none():
    taxonomy_dict = {"genera": None, "species": None}

    assert None is cazy_webscraper.get_filter_set(taxonomy_dict)
