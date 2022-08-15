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
import re

import pandas as pd

from datetime import datetime
from pathlib import Path

import pytest

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import Session
from sqlalchemy import create_engine


Base = declarative_base()


@pytest.fixture
def time_stamp():
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    return time_stamp

@pytest.fixture
def start_time():   
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)
    return start_time


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


@pytest.fixture
def mock_return_logger(*args, **kwards):
    return logging.getLogger('mock_logger')


@pytest.fixture
def mock_config_logger(*args, **kwargs):
    return


@pytest.fixture
def config_dict():
    configuration_dict = {
        'classes': [],
        "Glycoside Hydrolases (GHs)": ["GH3"],
        'GlycosylTransferases (GTs)': [],
        "Polysaccharide Lyases (PLs)": None,
        'Carbohydrate Esterases (CEs)': [],
        'Auxiliary Activities (AAs)': [],
        'Carbohydrate-Binding Modules (CBMs)': [],
    }
    return configuration_dict


@pytest.fixture()
def mock_return_none(*args, **kwargs):
    return


# Define fixtures for testing SQL ORM and interface


@pytest.fixture()
def db_path():
    return Path("tests/test_inputs/unit_test_database/unit_test_2022_08_11.db")


@pytest.fixture(scope="function")
def engine(db_path):
    return create_engine(f"sqlite+pysqlite:///{db_path}", echo=False)


@pytest.fixture(scope="function")
def tables(engine):
    Base.metadata.create_all(engine)
    yield
    Base.metadata.drop_all(engine)


@pytest.fixture
def db_connection(engine, tables):
    """Returns an sqlalchemy session, and after the test tears down everything properly."""
    connection = engine.connect()
    # begin the nested transaction
    transaction = connection.begin()

    yield transaction

    # roll back the broader transaction
    transaction.rollback()
    # put back the connection to the connection pool
    connection.close()
