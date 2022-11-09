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
from bs4 import BeautifulSoup
from pathlib import Path
from requests.exceptions import MissingSchema

from cazy_webscraper.crawler import get_validation_data
from cazy_webscraper.crawler.get_validation_data import CazyClass


@pytest.fixture
def args():
    args = {"args": Namespace(
        retries=2,
        timeout=45,
    )}
    return args


@pytest.fixture
def args_subfam_true():
    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=45,
        )
    }
    return argsdict


@pytest.fixture
def args_subfam_false():
    argsdict = {
        "args": Namespace(
            subfamilies=False,
            retries=2,
            timeout=45,
        )
    }
    return argsdict


@pytest.fixture
def cazy_url():
    return "html/www.cazy.org"


@pytest.fixture
def cazy_home_page(input_dir):
    file_path = input_dir / "class_url_pages" / "cazy_homepage.html"
    return file_path


@pytest.fixture
def cache_dir():
    return Path("tests/test_outputs/test_outputs_validation_data")


@pytest.fixture
def cazy_class_page_no_urls(input_dir):
    file_path = input_dir / "family_url_pages" / "cazy_classpage_no_urls.html"
    return file_path


@pytest.fixture
def cazy_home_no_spip(input_dir):
    file_path = input_dir / "class_url_pages" / "cazy_homepage_no_spip_out.html"
    return file_path


@pytest.fixture
def cazy_home_no_urls(input_dir):
    file_path = input_dir / "class_url_pages" / "cazy_homepage_no_urls.html"
    return file_path


@pytest.fixture
def cazy_class_page(input_dir):
    file_path = input_dir / "family_url_pages" / "cazy_classpage.html"
    return file_path


@pytest.fixture
def family_urls(input_dir):
    file_path = input_dir / "test_family_urls.txt"

    with open(file_path, "r") as fh:
        fam_lines = fh.read().splitlines()

    fam_list = []
    for line in fam_lines:
        fam_list.append([line, 0])

    return fam_list


@pytest.fixture
def family_h3_element(cazy_class_page):
    with open(cazy_class_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    return [_ for _ in
            soup.find_all("h3", {"class": "spip"}) if
            str(_.contents[0]) == "Tables for Direct Access"][0]


@pytest.fixture
def input_dir(test_input_dir):
    dir_path = test_input_dir / "test_inputs_crawler"
    return dir_path


@pytest.fixture
def no_subfam_h3_element(input_dir):
    file_path = input_dir / "cazy_classpage_no_subfams.html"
    with open(file_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    return [_ for _ in
            soup.find_all("h3", {"class": "spip"}) if
            str(_.contents[0]) == "Tables for Direct Access"][0]


@pytest.fixture
def subfamily_urls(input_dir):
    file_path = input_dir / "subfamily_urls.txt"
    with open(file_path, "r") as fh:
        fam_string = fh.read()
    fam_string = fam_string[1:-1]
    fam_string = fam_string.replace("'", "")
    fam_list = fam_string.split(", ")
    return fam_list

###### Test functions


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


# test get_cazy_classes


def test_get_class_urls_fail(cazy_url, cazy_dictionary, monkeypatch, args, cache_dir, start_time):
    """Test get_cazy_class_urls home_page not returned"""

    def mock_get_home_page(*args, **kwargs):
        return None, "error"

    monkeypatch.setattr(get_validation_data, "get_page", mock_get_home_page)

    assert get_validation_data.get_cazy_classes(
            cazy_url,
            None,
            cazy_dictionary,
            cache_dir,
            start_time,
            args["args"],
            unit_test=True,
        ) is None


def test_get_class_urls_exclusions_none(
    cazy_url,
    cache_dir,
    cazy_home_page,
    cazy_dictionary,
    monkeypatch,
    start_time,
    args,
):
    """Test get_cazy_class_urls when excluded_classess is None."""
    with open(cazy_home_page, "r") as fp:
        home_page = BeautifulSoup(fp, features="lxml")

        def mock_get_home_page(*args, **kwargs):
            return [home_page, None]

    monkeypatch.setattr(get_validation_data, "get_page", mock_get_home_page)

    result = get_validation_data.get_cazy_classes(
            cazy_url,
            None,
            cazy_dictionary,
            cache_dir,
            start_time,
            args["args"],
            unit_test=True,
    )
    assert len(result) == 6


def test_get_class_urls_exclusions_given(
    cazy_url,
    cache_dir,
    cazy_home_page,
    cazy_dictionary,
    monkeypatch,
    start_time,
    args,
):
    """Test get_cazy_class_urls when excluded_classess is not None."""
    with open(cazy_home_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    exclusions = ["<strong>Glycoside Hydrolases (GHs)</strong>"]

    def mock_get_home_page(*args, **kwargs):
        return [soup, None]

    monkeypatch.setattr(get_validation_data, "get_page", mock_get_home_page)

    result = get_validation_data.get_cazy_classes(
        cazy_url,
        exclusions,
        cazy_dictionary,
        cache_dir,
        start_time,
        args["args"],
        unit_test=True,
    )

    assert len(result) == 5


def test_get_class_urls_attribute(
    cazy_url,
    cache_dir,
    cazy_home_no_spip,
    cazy_dictionary,
    monkeypatch,
    start_time,
    args,
):
    """Test get_cazy_class_urls when attribute error is raised."""
    with open(cazy_home_no_spip) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    exclusions = ["<strong>Glycoside Hydrolases (GHs)</strong>"]

    def mock_get_home_page(*args, **kwargs):
        return [soup, None]

    monkeypatch.setattr(get_validation_data, "get_page", mock_get_home_page)

    assert get_validation_data.get_cazy_classes(
        cazy_url,
        exclusions,
        cazy_dictionary,
        cache_dir,
        start_time,
        args["args"],
        unit_test=True,
        ) is None


def test_get_class_urls_no_urls(
    cazy_url,
    cache_dir,
    cazy_home_no_urls,
    cazy_dictionary,
    monkeypatch,
    start_time,
    args,
):
    """Test get_cazy_class_urls when no class urls are returned from the HTML webpage."""
    with open(cazy_home_no_urls) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    exclusions = ["<strong>Glycoside Hydrolases (GHs)</strong>"]

    def mock_get_home_page(*args, **kwargs):
        return [soup, None]

    monkeypatch.setattr(get_validation_data, "get_page", mock_get_home_page)

    assert get_validation_data.get_cazy_classes(
        cazy_url,
        exclusions,
        cazy_dictionary,
        cache_dir,
        start_time,
        args["args"],
        unit_test=True,
        ) is None


# test get_families_urls


def test_get_family_urls_no_urls(cazy_url, args_subfam_false, monkeypatch, cazy_class_page_no_urls):
    """Tests get_families_urls when no Family URls are returned."""
    with open(cazy_class_page_no_urls) as fp:
        page = BeautifulSoup(fp, features="lxml")

    fam, message, incorrect_urls = get_validation_data.get_families_urls(
        cazy_url,
        "Glycoside Hydrolases (GHs)",
        page,
        args_subfam_false["args"],
    )
    assert fam is None
    assert message == "Failed to retrieve URLs to CAZy families for Glycoside Hydrolases (GHs)"
    assert incorrect_urls == []


def test_get_family_urls_no_urls_no_subfam_true(
    cazy_url,
    args_subfam_true,
    monkeypatch,
    cazy_class_page_no_urls,
):
    """Tests get_families_urls when no Family URls are returned, and retrieving subfamiles
    but none are retrieved."""
    with open(cazy_class_page_no_urls) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return page, None

    def mock_subfams(*args, **kwargs):
        return

    monkeypatch.setattr(get_validation_data, "get_subfamily_links", mock_subfams)

    fam, message, incorrect_urls = get_validation_data.get_families_urls(
        cazy_url,
        "Glycoside Hydrolases (GHs)",
        page,
        args_subfam_true["args"],
    )

    assert fam is None
    assert incorrect_urls == []


def test_get_family_urls_no_urls_subfam_true(
    cazy_url,
    args_subfam_true,
    monkeypatch,
    cazy_class_page_no_urls,
):
    """Tests get_families_urls when no Family URls are returned, and retrieving subfamiles
    and are retrieved."""
    with open(cazy_class_page_no_urls) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_subfams(*args, **kwargs):
        return ["http://www.cazy.org/GH5_1.html"]

    monkeypatch.setattr(get_validation_data, "get_subfamily_links", mock_subfams)

    fam, message, incorrect_urls = get_validation_data.get_families_urls(
        cazy_url,
        "Glycoside Hydrolases (GHs)",
        page,
        args_subfam_true["args"],
    )


def test_get_family_urls_success(
    cazy_class_page,
    args_subfam_true,
    monkeypatch,
):
    """Test get_families_urls when successful, and subfamilies is True."""
    with open(cazy_class_page) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_get_subfams(*args, **kwargs):
        return []

    monkeypatch.setattr(get_validation_data, "get_subfamily_links", mock_get_subfams)

    fam, message, incorrect_urls = get_validation_data.get_families_urls(
        cazy_url,
        "Glycoside Hydrolases (GHs)",
        page,
        args_subfam_true["args"],
    )

    print(incorrect_urls)

    assert incorrect_urls == []


# test get_subfamily_links


def test_get_subfam_links_successul(family_h3_element, subfamily_urls):
    """Test get_subfamily_links when links are retrieved."""

    res = get_validation_data.get_subfamily_links(
        family_h3_element,
        "http://www.cazy.org",
    )
    print(subfamily_urls)
    assert res == subfamily_urls


def test_get_subfam_links_no_urls(no_subfam_h3_element):
    """Test get_subfamily_links when no urls are retrieved."""

    assert None is get_validation_data.get_subfamily_links(
        no_subfam_h3_element,
        "http://www.cazy.org",
    )


# test browser dectorator

def test_browser_decorator():
    """Test browser_decorator to ensure proper handling if unsuccessful."""
    args = {"args": Namespace(timeout=10)}
    result = get_validation_data.get_page('www.caz!!!!!!!!y.org', args["args"], max_tries=1)
    assert True == (result[0] is None) and (type(result[1]) is MissingSchema)


# test CAZyClass


def test_cazy_class():
    new_class = get_validation_data.CazyClass("GH", "url", 0)
    new_class = get_validation_data.CazyClass("GH", "url", 0, {})


# test get fam populations


def test_get_fam_pops_failed_connection(
    cazy_url,
    cache_dir,
    start_time,
    monkeypatch,
):
    """Unit test get_cazy_family_pops() when can't get page"""
    argsdict = {
        "args": Namespace(
            force=True,
            nodelete_cache=True,
            retries=2,
        )
    }

    def mock_get_page(*args, **kwards):
        return None, "error"

    monkeypatch.setattr(get_validation_data, "get_page", mock_get_page)

    result, error, incorrect_urls, failed_connections = get_validation_data.get_cazy_family_pops(
        "GH",
        "www.cazy_url",
        cazy_url,
        ["GH1"],
        cache_dir,
        start_time,
        argsdict["args"],
        unit_test=True,
    )

    assert result is None
    assert incorrect_urls == []
    assert failed_connections == []


def test_get_fam_pops_no_fams(
    cazy_url,
    cache_dir,
    start_time,
    cazy_class_page,
    monkeypatch,
):
    """Unit test get_cazy_family_pops() when no fam urls are retrieved"""
    with open(cazy_class_page) as fp:
        page = BeautifulSoup(fp, features="lxml")

    argsdict = {
        "args": Namespace(
            force=True,
            nodelete_cache=True,
            retries=2,
        )
    }

    def mock_get_page(*args, **kwards):
        return page, None

    def mock_get_fams(*args, **kwards):
        return None, "error", []

    monkeypatch.setattr(get_validation_data, "get_page", mock_get_page)
    monkeypatch.setattr(get_validation_data, "get_families_urls", mock_get_fams)

    result, error, incorrect_urls, failed_connections = get_validation_data.get_cazy_family_pops(
        "GH",
        "www.cazy_url",
        cazy_url,
        ["GH1"],
        cache_dir,
        start_time,
        argsdict["args"],
        unit_test=True,
    )

    assert result is None
    assert error == "error"
    assert incorrect_urls == []
    assert failed_connections == []


def test_get_fam_pops_success(
    cazy_url,
    cache_dir,
    start_time,
    cazy_class_page,
    monkeypatch,
):
    """Unit test get_cazy_family_pops() when no fam urls are retrieved"""
    with open(cazy_class_page) as fp:
        page = BeautifulSoup(fp, features="lxml")

    argsdict = {
        "args": Namespace(
            force=True,
            nodelete_cache=True,
            retries=2,
        )
    }

    def mock_get_page(*args, **kwards):
        return page, None

    def mock_get_fams(*args, **kwards):
        return ['http://www.cazy.org/GH50.html', 'http://www.cazy.org/GH1.html'], "error", []

    monkeypatch.setattr(get_validation_data, "get_page", mock_get_page)
    monkeypatch.setattr(get_validation_data, "get_families_urls", mock_get_fams)

    result, error, incorrect_urls, failed_connections = get_validation_data.get_cazy_family_pops(
        "GH",
        "www.cazy_url",
        cazy_url,
        ["GH1"],
        cache_dir,
        start_time,
        argsdict["args"],
        unit_test=True,
    )

    assert result == ['http://www.cazy.org/GH50.html', 'http://www.cazy.org/GH1.html']
    assert error == "error"
    assert incorrect_urls == []
    assert failed_connections == []
