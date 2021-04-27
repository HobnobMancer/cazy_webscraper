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
"""Tests the script parse_local_pages which scrapes CAZyme data from a local CAZy page library.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import pytest
import sys

import pandas as pd

from argparse import Namespace
from datetime import datetime

from bs4 import BeautifulSoup

from scraper import cazy_webscraper, crawler, sql, utilities
from scraper import crawler
from scraper.crawler import Family
from scraper.crawler.cazy_html_pages import parse_local_pages
from scraper.utilities import file_io, parse_configuration, parsers


@pytest.fixture
def input_dir(test_dir):
    path_ = test_dir / "test_inputs" / "test_inputs_crawler" / "parse_local_pages"
    return path_


@pytest.fixture
def args_no_files(test_dir):
    path_ = test_dir / "test_inputs" / "test_inputs_crawler" / "get_cazy_pages"
    out = test_dir / "test_outputs" / "test_outputs_parse_local_pages"

    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
            scrape_files=path_,
            output=out,
        )
    }
    return argsdict


@pytest.fixture
def args(test_dir):
    path_ = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_all_inputs" / "test_get_pagination_urls"

    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
            scrape_files=path_,
            output=sys.stdout,
        )
    }
    return argsdict


@pytest.fixture
def args_parse(test_dir, input_dir):
    path_ = test_dir / "test_outputs" / "test_outputs_parse_local_pages"
    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
            scrape_files=input_dir,
            output=path_,
        )
    }
    return argsdict


# test parse_local_pages


def test_local_parse(monkeypatch, input_dir, test_dir, args_parse, cazy_home_url):

    path_ = test_dir / "test_outputs" / "test_outputs_parse_local_pages"

    file_io.make_output_directory(path_, True, False)

    no_cazyme_table = input_dir / "http___www_cazy_org_GH60_all_html.html"
    incorrect_format = input_dir / "incorrect_format.html"
    deleted_fam = input_dir / "http___www_cazy_org_GH21_all_html.html"
    empty_fam = input_dir / "http___www_cazy_org_GH61_all_html.html"

    def mock_get_html_paths(*args, **kwargs):
        return [no_cazyme_table, incorrect_format, deleted_fam, empty_fam]

    def mock_sql(*args, **kwargs):
        return {"sql": "error message"}

    monkeypatch.setattr(parse_local_pages, "get_html_files", mock_get_html_paths)
    monkeypatch.setattr(crawler, "row_to_protein", mock_sql)

    time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_local_pages.parse_local_pages(
            args=args_parse["args"],
            cazy_home=cazy_home_url,
            start_time=start_time,
            time_stamp=time_stamp,
            session="session",
            taxonomy_filters=set(),
            ec_filters=set(),
        )
    assert pytest_wrapped_e.type == SystemExit

    file_io.make_output_directory(path_, True, False)


# test get_html_files


def test_local_no_files(test_dir, args_no_files):
    """Test get_html_files when no files are found."""

    path_ = test_dir / "test_inputs" / "test_inputs_crawler" / "get_cazy_pages"

    file_io.make_output_directory(path_, True, False)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_local_pages.get_html_files(args_no_files["args"])
    assert pytest_wrapped_e.type == SystemExit


def test_local_files(test_dir, args):
    """Test get_html_files when no files are found."""

    parse_local_pages.get_html_files(args["args"])
