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

    sql_interface.add_protein_to_db(
        "test_cazyme",
        "test_fam",
        "test_genus test_species",
        "test_primary_genbank",
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
