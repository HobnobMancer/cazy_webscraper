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


from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker


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


# define fixtures for testing SQL orm and interface, these are written here becuase the _addoption
# needs to be in the conftest.py, and it keeps the related fixtures together


def pytest_addoption(parser):
    db_path = "tests/test_outputs/test_outputs_sql/cazy_scrape_2021-02-08--10-58-51.db"
    parser.addoption(
        '--dburl',
        action='store', 
        default=f"sqlite+pysqlite:///{db_path}",
        help='url of the db for unit tests',
    )


@ pytest.fixture(scope='session')
def db_engine(request):
    """Yield a SQLAlchemy engine that is suppressed after the test session ends."""
    db_path = request.config.getoption("--dburl")
    engine_ = create_engine(db_path, echo=True)
    yield engine_
    engine_.dispose()


@pytest.fixture(scope='session')
def db_session_factory(db_engine):
    """Return a SQLAlchemy scoped session factory."""
    return scoped_session(sessionmaker(bind=db_engine))


@pytest.fixture(scope='function')
def db_session(db_session_factory):
    """Yield a SQLAlchemy connection that is rolledback after the unit test is complete."""
    session_ = db_session_factory()
    yield session_
    session_.rollback()
    session_.close()
