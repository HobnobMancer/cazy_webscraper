#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2026
# (c) University of Strathclyde 2026
# (c) James Hutton Institute 2026
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# ehobbs@ebi.ac.uk
#
# The MIT License
"""Shared pytest fixtures for the v3 (src/) test suite.

Every fixture here targets src/ - there is no v2 (legacy cazy_webscraper package) coverage in
this suite any more.
"""


import logging

from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

from src.sql import sql_orm


@pytest.fixture
def test_dir():
    return Path("tests/")


@pytest.fixture
def test_input_dir(test_dir):
    return test_dir / "test_inputs"


@pytest.fixture
def null_logger():
    logger = logging.getLogger("test_null_logger")
    logger.addHandler(logging.NullHandler())
    return logger


@pytest.fixture
def time_stamp():
    """A timestamp in the format src/scripts/*.py coordination functions expect."""
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


@pytest.fixture
def start_time():
    return pd.to_datetime(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


@pytest.fixture
def db_path(tmp_path):
    """Path to a fresh, empty v3-schema local CAZyme database.

    Built directly from the current SQLAlchemy models (src.sql.sql_orm), never from a
    bundled .db file - the only one bundled in test_inputs/ (unit_test_database/) is v2
    schema and cannot be opened by v3 code.
    """
    path = tmp_path / "test.db"
    sql_orm.get_db_connection(path, False, new=True)
    return path
