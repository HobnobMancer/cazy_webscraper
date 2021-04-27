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
"""Tests the module sql which builds and interacts with an SQL database.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest
import shutil
import sys

from argparse import Namespace
from datetime import datetime
from pathlib import Path

from sqlalchemy.orm.exc import ObjectDeletedError

from scraper.utilities import file_io
from scraper.sql import sql_orm
from scraper.sql.sql_orm import (
    Cazyme,
    CazyFamily,
    Genbank,
    Taxonomy,
    Kingdom,
    EC,
    Uniprot,
    Pdb,
    Log,
    Cazymes_Genbanks,
)


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


def test_calling_class_instances():
    """Test calling Class instances strs and reprs."""
    cazyme = Cazyme(cazyme_name="test_cazyme")
    str(cazyme)
    repr(cazyme)

    tax = Taxonomy(genus="genus", species="species")
    str(tax)
    repr(tax)

    kng = Kingdom(kingdom="kingdom")
    str(kng)
    repr(kng)

    fam = CazyFamily(family="family")
    str(fam)
    repr(fam)

    fam = CazyFamily(family="family", subfamily="subfamily")
    str(fam)
    repr(fam)

    gb = Genbank(genbank_accession="accession")
    str(gb)
    repr(gb)

    cg = Cazymes_Genbanks(cazyme_id=1, genbank_id=2)
    str(cg)
    repr(cg)

    ec = EC(ec_number="EC1.2.3.4")
    str(ec)
    repr(ec)

    up = Uniprot(uniprot_accession="accession")
    str(up)
    repr(up)

    pdb = Pdb(pdb_accession="accession")
    str(pdb)
    repr(pdb)

    log = Log(date="date", time="time", classes="classes", families="families")
    str(log)
    repr(log)


def test_build_db():
    """Test building a database."""
    path_ = Path("tests")
    path_ = path_ / "test_outputs" / "test_outputs_sql" / "temp_dir"
    file_io.make_output_directory(path_, True, False)
    args = {"args": Namespace(output=path_)}
    sql_orm.build_db("time_stamp", args["args"])
    shutil.rmtree(path_)


def test_build_db_stdout():
    """Test building a database when output is stdout"""
    path_ = Path("tests")
    path_ = path_ / "test_outputs" / "test_outputs_sql" / "temp_dir"
    file_io.make_output_directory(path_, True, False)
    args = {"args": Namespace(output=sys.stdout)}
    sql_orm.build_db("time_stamp", args["args"])
    shutil.rmtree(path_)
