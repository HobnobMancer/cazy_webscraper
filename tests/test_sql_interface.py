#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Tests sql_interface.__init__.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


from argparse import Namespace
from datetime import datetime

import pytest

from sqlalchemy.exc import IntegrityError

from cazy_webscraper.sql import sql_interface, sql_orm


@pytest.fixture
def argsdict():
    args_data = {
        "args": Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
            cazy_synonyms=None,
            classes="PL,CE",
            config=None,
            citation=False,
            db_output="path",
            database=None,
            delete_old_relationships=False,
            ec_filter="EC3.2.1.3,EC2.5.4.6",
            force=False,
            families="CE1,CE2",
            genera="Aspergillus",
            kingdoms="Viruses",
            log=None,
            nodelete=False,
            nodelete_cache=False,
            nodelete_log=False,
            retries=10,
            sql_echo=False,
            subfamilies="PL1_1",
            species="soecies",
            strains="strains",
            timeout=45,
            validate=True,
            verbose=False,
            version=False,
        )
    }
    return args_data


def test_log_scrape(db_path, argsdict):
    """Test log_scrape_in_db()"""
    time_stamp = "time_Stamp"
    config_dict = {
        "classes": {'GH'},
        "GH": {"GH1"},
    }
    taxonomy_filters = {"genera": "Trichoderma", "species": {"1", "2"}, "strains": {"1"}}
    kingdoms = {"Bacteria"}
    ec_filter = {"32.15.1.2"}
    retrieved_annotations = "unit test"
    db = "unit test"

    db_connection = sql_orm.get_db_connection(db_path, False, False)

    with sql_orm.Session(bind=db_connection) as session:
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            taxonomy_filters,
            kingdoms,
            ec_filter,
            db,
            retrieved_annotations,
            session,
            argsdict["args"],
        )

        session.rollback()


def test_log_scrape_class_error(db_path, argsdict):
    """Test log_scrape_in_db()"""
    time_stamp = "time_Stamp"
    config_dict = {
        "GH": {"GH1"},
    }
    taxonomy_filters = {"genera": "Trichoderma", "species": {"1", "2"}, "strains": {"1"}}
    kingdoms = []
    ec_filter = []
    retrieved_annotations = "unit test"
    db = "unit test"

    db_connection = sql_orm.get_db_connection(db_path, False, False)

    with sql_orm.Session(bind=db_connection) as session:
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            taxonomy_filters,
            kingdoms,
            ec_filter,
            db,
            retrieved_annotations,
            session,
            argsdict["args"],
        )

        session.rollback()


def test_log_scrape_no_classes(db_path, argsdict):
    """Test log_scrape_in_db()"""
    time_stamp = "time_Stamp"
    config_dict = {
        "classes": [],
        "GH": {"GH1", "GH1"},
        "CE": set(),
    }
    taxonomy_filters = dict()
    kingdoms = {"Bacteria", "Viruses"}
    ec_filter = {"32.15.1.2"}
    retrieved_annotations = "unit test"
    db = "unit test"

    db_connection = sql_orm.get_db_connection(db_path, False, False)

    with sql_orm.Session(bind=db_connection) as session:
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            taxonomy_filters,
            kingdoms,
            ec_filter,
            db,
            retrieved_annotations,
            session,
            argsdict["args"],
        )

        session.rollback()


def test_insert_data(db_path, argsdict):
    """Test insert_data"""
    db_connection = sql_orm.get_db_connection(db_path, False, False)

    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    sql_interface.insert_data(db_connection, "Kingdoms", ["Kingdom"], [(time_stamp,)])

    with sql_orm.Session(bind=db_connection) as session:
        session.rollback()


def test_get_table(db_path, argsdict):
    """Test get_gbk_table_dict"""
    db_connection = sql_orm.get_db_connection(db_path, False, False)

    sql_interface.get_gbk_table_dict(db_connection)
