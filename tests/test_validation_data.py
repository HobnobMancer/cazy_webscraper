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
"""Tests crawler.get_validation_data.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import logging
import pytest

from argparse import Namespace
from pathlib import Path

from cazy_webscraper.crawler import get_validation_data
from cazy_webscraper.crawler.get_validation_data import CazyClass


@pytest.fixture
def cazy_url():
    return "html/www.cazy.org"


@pytest.fixture
def cache_dir():
    return Path("tests/test_outputs/test_outputs_validation_data")


def test_get_val_data_no_classes(start_time, cazy_url, cache_dir, monkeypatch):
    """Test get_validation_data when no CAZy class URLs are retrieved"""

    argsdict = {
        "args": Namespace(
            force=True,
            nodelete_cache=True,
        )
    }
    # mock making output dir

    # mock get_cazy_classes, return none
    def mock_return_none(*args, **kwards):
        return

    monkeypatch.setattr(get_validation_data, "make_output_directory", mock_return_none)
    monkeypatch.setattr(get_validation_data, "get_cazy_classes", mock_return_none)

    get_validation_data.get_validation_data(
        cazy_url,
        set(),
        {},
        {},
        cache_dir,
        logging.getLogger(),
        start_time,
        argsdict['args'],
    )


def test_get_val_data_success(start_time, cazy_url, cache_dir, monkeypatch):
    """Test get_validation_data when successful"""

    argsdict = {
        "args": Namespace(
            force=True,
            nodelete_cache=True,
        )
    }
    # mock making output dir

    classes = [
        CazyClass("GH", "www.cazy.org", 0),
    ]

    config_dict = {"GH": ["GH1", "GH2"]}

    incorrect_urls = ["incorrect_url"]

    fam_pop = {"GH1": 1, "GH2": 2}

    err_message = ""

    # mock get_cazy_classes, return none
    def mock_return_none(*args, **kwards):
        return

    def mock_classes(*args, **kwards):
        return classes

    def mock_get_fam_pop(*args, **kwards):
        return fam_pop, err_message, incorrect_urls, None

    monkeypatch.setattr(get_validation_data, "make_output_directory", mock_return_none)
    monkeypatch.setattr(get_validation_data, "get_cazy_classes", mock_classes)
    monkeypatch.setattr(get_validation_data, "get_cazy_family_pops", mock_get_fam_pop)

    get_validation_data.get_validation_data(
        cazy_url,
        set(),
        {},
        config_dict,
        cache_dir,
        logging.getLogger(),
        start_time,
        argsdict['args'],
    )


def test_get_val_data_failed_fam(start_time, cazy_url, cache_dir, monkeypatch):
    """Test get_validation_data when no can't get family population"""
    argsdict = {
        "args": Namespace(
            force=True,
            nodelete_cache=True,
            retries=2,
        )
    }
    # mock making output dir

    classes = [
        CazyClass("GH", "www.cazy.org", 0),
    ]
    config_dict = {"GH": ["GH1", "GH2"]}
    incorrect_urls = None
    fam_pop = None
    err_message = "message"
    failed_fams = ["GH1", "GH2"]

    # mock get_cazy_classes, return none
    def mock_return_none(*args, **kwards):
        return

    def mock_classes(*args, **kwards):
        return classes

    def mock_get_fam_pop(*args, **kwards):
        return fam_pop, err_message, incorrect_urls, failed_fams

    monkeypatch.setattr(get_validation_data, "make_output_directory", mock_return_none)
    monkeypatch.setattr(get_validation_data, "get_cazy_classes", mock_classes)
    monkeypatch.setattr(get_validation_data, "get_cazy_family_pops", mock_get_fam_pop)

    get_validation_data.get_validation_data(
        cazy_url,
        set(),
        {},
        config_dict,
        cache_dir,
        logging.getLogger(),
        start_time,
        argsdict['args'],
    )