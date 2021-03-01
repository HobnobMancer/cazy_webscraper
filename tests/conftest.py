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
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine


Base = declarative_base()
Session = sessionmaker()


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
def db_session():
    """Open session to local SQL database."""
    db_path = "tests/test_outputs/test_outputs_sql/cazy_scrape_2021-02-08--10-58-51.db"
    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=False)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)

    return Session()
