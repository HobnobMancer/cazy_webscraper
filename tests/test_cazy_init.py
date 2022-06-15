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
"""Tests cazy module

These test are intened to be run from the root of the repository using:
pytest -v
"""


import json
import pandas as pd

from argparse import Namespace, ArgumentParser
from datetime import datetime
from pathlib import Path

import pytest

from sqlalchemy.exc import IntegrityError
from saintBioutils.utilities import logger as saint_logger
from saintBioutils.utilities import file_io as saint_fileIO

from cazy_webscraper import cazy

@pytest.fixture
def cazy_file_path():
    return "tests/test_inputs/test_inputs_cazy/cazy_data.txt"


@pytest.fixture
def cazy_zip_path():
    _path = "tests/test_inputs/test_inputs_cazy/cazy_db_timestamp.zip"
    return _path


def test_get_cazy_file(cazy_file_path):
    argsdict = {"args": Namespace(
        retries=10,
        cazy_data=cazy_file_path,
    )}
    cazy.get_cazy_txt_file_data(
        "tests/test_inputs/test_inputs_cazy/",
        "time_stamp",
        argsdict['args'],
    )


def test_parsing_cazy_zip(monkeypatch):
    argsdict = {"args": Namespace(
        retries=10,
        cazy_data=None,
    )}

    def mock_download(*args, **kwards):
        return

    monkeypatch.setattr(cazy, "get_cazy_file", mock_download)

    cazy.get_cazy_txt_file_data(
        Path("tests/test_inputs/test_inputs_cazy/"),
        "time_stamp",
        argsdict['args'],
    )
