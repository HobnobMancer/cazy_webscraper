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

"""Tests the module sql which builds and interacts with an SQL database.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace, ArgumentParser
from pathlib import Path

from scraper import expand, utilities, file_io
from scraper.expand import get_pdb_structures, get_genbank_sequences


@pytest.fixture
def db_path():
    path_ = Path("tests")
    path_ = path_ / "test_inputs" / "test_inputs_sql" / "unit_test_db_2021-03-01--15-06-59.db"
    return path_


@pytest.fixture
def output_dir(test_dir):
    path_ = test_dir / "test_outputs"
    return path_


# tests for get_pdb_structures


def test_main_no_db(monkeypatch):
    """Test main() when an the database file cannot be found."""

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
            database=Path("--"),
            verbose=False,
            log=None,
            force=False,
            nodelete=False,
            outdir=None,
        )
        return parser

    def mock_no_return(*args, **kwargs):
        return

    def mock_config(*args, **kwargs):
        return None, set()

    monkeypatch.setattr(utilities, "build_pdb_structures_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(file_io, "get_configuration", mock_config)
    monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_pdb_structures.main()
    assert pytest_wrapped_e.type == SystemExit


def test_main_outdir_is_none(db_path, monkeypatch):
    """Test main() when outdir=None."""

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
            database=db_path,
            outdir=None,
            verbose=False,
            log=None,
            force=False,
            nodelete=False,
        )
        return parser

    def mock_no_return(*args, **kwargs):
        return

    def mock_config(*args, **kwargs):
        return None, set()

    monkeypatch.setattr(utilities, "build_pdb_structures_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_no_return)
    monkeypatch.setattr(file_io, "get_configuration", mock_config)
    monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
    monkeypatch.setattr(get_pdb_structures, "get_every_cazymes_structures", mock_no_return)

    get_pdb_structures.main()


def test_main(db_path, output_dir, monkeypatch):
    """Test main()."""

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
            database=db_path,
            outdir=output_dir,
            verbose=False,
            log=None,
            force=True,
            nodelete=True,
        )
        return parser

    def mock_no_return(*args, **kwargs):
        return

    def mock_config(*args, **kwargs):
        return None, set()

    monkeypatch.setattr(utilities, "build_pdb_structures_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_no_return)
    monkeypatch.setattr(file_io, "get_configuration", mock_config)
    monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
    monkeypatch.setattr(get_pdb_structures, "get_every_cazymes_structures", mock_no_return)

    get_pdb_structures.main()


# test for get_genbank_sequences


def test_main_no_db_genbank(monkeypatch):
    """Test main() when an the database file cannot be found."""

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
            database=Path("--"),
            email="dummy_email",
            verbose=False,
            log=None,
            force=False,
            nodelete=False,
            outdir=None,
        )
        return parser

    def mock_no_return(*args, **kwargs):
        return

    def mock_config(*args, **kwargs):
        return None, set()

    monkeypatch.setattr(utilities, "build_genbank_sequences_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(file_io, "get_configuration", mock_config)
    monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_genbank_sequences.main()
    assert pytest_wrapped_e.type == SystemExit


def test_main_no_config_no_update(db_path, monkeypatch):
    """Test main() when outdir=None."""

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
            database=db_path,
            email="dummy_email",
            outdir=None,
            verbose=False,
            log=None,
            force=False,
            nodelete=False,
            update=False,
        )
        return parser

    def mock_no_return(*args, **kwargs):
        return

    def mock_config(*args, **kwargs):
        return None, set()

    monkeypatch.setattr(utilities, "build_genbank_sequences_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_no_return)
    monkeypatch.setattr(file_io, "get_configuration", mock_config)
    monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
    monkeypatch.setattr(
        get_genbank_sequences,
        "get_missing_sequences_for_everything",
        mock_no_return,
    )

    get_genbank_sequences.main()


def test_main_no_config_yes_update(db_path, monkeypatch):
    """Test main() when outdir=None."""

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
            database=db_path,
            email="dummy_email",
            outdir=None,
            verbose=False,
            log=None,
            force=False,
            nodelete=False,
            update=True,
        )
        return parser

    def mock_no_return(*args, **kwargs):
        return

    def mock_config(*args, **kwargs):
        return None, set()

    monkeypatch.setattr(utilities, "build_genbank_sequences_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_no_return)
    monkeypatch.setattr(file_io, "get_configuration", mock_config)
    monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
    monkeypatch.setattr(
        get_genbank_sequences,
        "add_and_update_all_sequences",
        mock_no_return,
    )

    get_genbank_sequences.main()


def test_main_yes_config_no_update(db_path, monkeypatch):
    """Test main() when outdir=None."""

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
            database=db_path,
            email="dummy_email",
            outdir=None,
            verbose=False,
            log=None,
            force=False,
            nodelete=False,
            update=False,
        )
        return parser

    def mock_no_return(*args, **kwargs):
        return

    def mock_config(*args, **kwargs):
        return {}, set()

    monkeypatch.setattr(utilities, "build_genbank_sequences_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_no_return)
    monkeypatch.setattr(file_io, "get_configuration", mock_config)
    monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
    monkeypatch.setattr(
        get_genbank_sequences,
        "get_missing_sequences_for_specific_records",
        mock_no_return,
    )

    get_genbank_sequences.main()


def test_main_yes_config_yes_update(db_path, monkeypatch):
    """Test main() when outdir=None."""

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
            database=db_path,
            email="dummy_email",
            outdir=None,
            verbose=False,
            log=None,
            force=False,
            nodelete=False,
            update=True,
        )
        return parser

    def mock_no_return(*args, **kwargs):
        return

    def mock_config(*args, **kwargs):
        return {}, set()

    monkeypatch.setattr(utilities, "build_genbank_sequences_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(utilities, "config_logger", mock_no_return)
    monkeypatch.setattr(file_io, "get_configuration", mock_config)
    monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
    monkeypatch.setattr(
        get_genbank_sequences,
        "update_sequences_for_specific_records",
        mock_no_return,
    )

    get_genbank_sequences.main()
