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
"""Configuration file for pytest files.
Contains fixtures used by multiple test files.
"""


import logging
import json

from datetime import datetime
from pathlib import Path

import pytest

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import Session
from sqlalchemy import create_engine


Base = declarative_base()


@pytest.fixture
def test_dir():
    return Path("tests/")


@pytest.fixture
def test_input_dir(test_dir):
    dir_path = test_dir / "test_inputs"
    return dir_path


@pytest.fixture
def null_logger():
    logger = logging.getLogger("Test_ac_number_retrieval logger")
    logger.addHandler(logging.NullHandler())
    return logger


@pytest.fixture
def cazy_home_url():
    return "http://www.cazy.org"


@pytest.fixture
def cazy_dictionary(test_input_dir):
    dict_path = test_input_dir / "cazy_dictionary.json"
    with open(dict_path, "r") as fh:
        cazy_dict = json.load(fh)
    return cazy_dict


@pytest.fixture
def time_stamp():
    time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")  # used in naming files
    return time_stamp


# Define fixtures for testing SQL ORM and interface


@pytest.fixture(scope="session")
def engine():
    db_path = "tests/test_inputs/test_inputs_sql/unit_test_db_2021-03-01--15-06-59.db"
    return create_engine(f"sqlite+pysqlite:///{db_path}", echo=False)


@pytest.fixture(scope="session")
def tables(engine):
    Base.metadata.create_all(engine)
    yield
    Base.metadata.drop_all(engine)


@pytest.fixture
def db_session(engine, tables):
    """Returns an sqlalchemy session, and after the test tears down everything properly."""
    connection = engine.connect()
    # begin the nested transaction
    transaction = connection.begin()
    # use the connection with the already started transaction
    session = Session(bind=connection)

    yield session

    session.close()
    # roll back the broader transaction
    transaction.rollback()
    # put back the connection to the connection pool
    connection.close()
