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
"""Tests the script cazy_webscraper/__init__.py which coordinates connection to db.

These test are intened to be run from the root of the repository using:
pytest -v
"""

from argparse import Namespace

from logging import getLogger
import logging
import os

from cazy_webscraper import cazy
from cazy_webscraper.utilities.parse_configuration.cazy_class_synonym_dict import cazy_synonym_dict
from saintBioutils.utilities import logger as saint_logger
import pytest
import sys

from cazy_webscraper import connect_to_new_db, connect_existing_db


def test_no_existing_db(time_stamp, start_time, mock_config_logger, monkeypatch):
    """"Test connect_exiting_db when db does not exist"""
    argsdict = {
        "args": Namespace(
            database="fakePath",
            verbose=False,
        )
    }

    monkeypatch.setattr(saint_logger, "config_logger", mock_config_logger)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        connect_existing_db(argsdict['args'], time_stamp, start_time)
    assert pytest_wrapped_e.type == SystemExit