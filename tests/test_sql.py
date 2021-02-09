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

from argparse import Namespace
from datetime import datetime
from pathlib import Path

from scraper.sql import sql_orm, sql_interface
from scraper.sql.sql_orm import Cazyme, CazyFamily


@pytest.fixture
def output_dir(test_dir):
    path = Path("test_outputs")
    path = path / "test_sql_outputs"
    return path


@pytest.fixture
def database_args(output_dir):
    args_dict = {
        "args": Namespace(
            output=output_dir,
        )
    }
    return args_dict


@pytest.fixture
def existing_db_args(output_dir):
    db_path = output_dir / "cazy_scrape_2021-02-07--15-35-13.db"
    args_dict = {
        "args": Namespace(
            output=output_dir,
            database=db_path
        )
    }
    return args_dict


# Unit tests for sql_orm


def test_regex_search(db_session):
    """Test a regular expression search can be performed successfully."""

    cazy_class = "PL"
    class_query = db_session.query(Cazyme.cazyme_id).\
        join(CazyFamily, Cazyme.families).\
        filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
        all()

    assert len(class_query) == 24


# Unit tests for sql_interface


def test_adding_a_new_protein(db_session):
    """Test adding a new protein to the local database."""
    time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    sql_interface.add_protein_to_db(
        "test_cazyme",
        "test_fam",
        "test_genus test_species",
        time_stamp,
        db_session,
        ec_numbers=["EC1.2.3.4"],
        genbank_accessions=["Genbank1", "Genbank2"],
        uniprot_accessions=["Uni1", "Uni2"],
        pdb_accessions=["PDB1", "PDB2"]
    )


def test_add_data_to_an_existing_record_in_db(db_session):
    """Test adding data to an existing CAZyme in the local database."""

    sql_interface.add_protein_to_db(
        "JCM19301_1832",
        "PL28",
        "pallidilutea JCM 19301",
        "GAL67220.1",
        db_session,
        ec_numbers=["EC4.2.2.-"],
        genbank_accessions=["Genbank1", "Genbank2"],
        uniprot_accessions=["Uni1", "Uni2"],
    )


def test_genbank_no_cazymes(db_session, monkeypatch):
    """test adding a protein to a database, GenBank accession is found with no linked CAZymes."""

    def mock_adding_a_new_protein(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_new_protein_to_db", mock_adding_a_new_protein)

    existing_genbank_with_no_cazyme = "test_genbank_no_cazyme"

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        existing_genbank_with_no_cazyme,
        db_session,
        ec_numbers=["EC4.2.2.-"],
        genbank_accessions=["Genbank1", "Genbank2"],
        uniprot_accessions=["Uni1", "Uni2"],
    )


def test_one_genbank_multiple_cazymes(db_session, monkeypatch):
    """test adding protein to db when genbank is found linked to muliple CAZymes."""

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_data_to_protein_record", mock_add_protein_to_db)

    genbank_with_multiple_cazymes = "one_genbank_multi_cazymes"

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        genbank_with_multiple_cazymes,
        db_session,
    )


def test_multiple_genbanks_multiple_cazymes(db_session, monkeypatch):
    """test adding protein to db when finding multiple identical CAZymes and GenBank accesisons."""

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_data_to_protein_record", mock_add_protein_to_db)

    identical_genbank_accession = "identical_accession"

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        identical_genbank_accession,
        db_session,
    )


def test_multiple_genbanks_one_cazyme(db_session, monkeypatch):
    """test adding protien to db when identical GenBank accessions with one CAZyme link."""

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_data_to_protein_record", mock_add_protein_to_db)

    multiple_accession = "multiple_accession"

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        multiple_accession,
        db_session,
    )


def test_multiple_genbanks_no_cazymes(db_session, monkeypatch):
    """test adding protien to db when identical GenBank accessions are linked to no CAZymes."""

    identical_genbanks_no_cazymes = "multi_accession_no_caz"

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_new_protein_to_db", mock_add_protein_to_db)

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        identical_genbanks_no_cazymes,
        db_session,
    )
