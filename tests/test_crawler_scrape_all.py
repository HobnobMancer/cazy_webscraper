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
"""Tests the module for scrapping proteins from the 'all' pages of CAZy families.

The functions tested in the is script are located in scraper.crawler.parse_cazy_families.scrape_all

These test are intened to be run from the root of the repository using:
pytest -v
"""

import pytest

from argparse import Namespace

from bs4 import BeautifulSoup

from scraper import crawler
from scraper.crawler import Family
from scraper.crawler.parse_cazy_families import scrape_all
from scraper.sql import sql_interface


@pytest.fixture
def fam_url():
    url = "http://www.cazy.org/GH14_all.html"
    return url


@pytest.fixture
def args():
    dict_ = {'args': Namespace(retries=2, timeout=10)}
    return dict_


@pytest.fixture
def fam_template():
    fam = Family("famName", "CAZyClass", fam_url)
    return fam


@pytest.fixture
def input_dir(test_dir):
    path_ = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_all_inputs"
    return path_


@pytest.fixture
def protein_gen():
    result_list = [
        {"url": None, "format": None, "sql": None, "error": None},
        {"url": "www.cazy.org/GH1.html", "format": "no internet connection", "sql": None, "error": "error message"},
        {"url": None, "format": "sql error", "sql": "protein_name", "error": "error message"},
        {"url": "www.cazy.org/GH3.html", "format": "sql error", "sql": "protein_name", "error": "error message"},
    ]
    return (_ for _ in result_list)


# test parse_family_via_all_pages()


def test_all_parse_family_first_parse(cazy_home_url, monkeypatch):
    """Test parse_family_via_all_pages() when it's the first time."""

    fam = Family("testFam", "testClass", "testURL")

    def mock_parse_pagination_page(*args, **kwargs):
        return fam, None, None, None

    monkeypatch.setattr(scrape_all, "parse_first_paginiation_page", mock_parse_pagination_page)

    scrape_all.parse_family_via_all_pages(
        family=fam,
        cazy_home=cazy_home_url,
        taxonomy_filters=None,
        ec_filters=None,
        args=None,
        session="session",
    )


def test_all_parse_family_retry_paginiation(cazy_home_url, monkeypatch):
    """Test parse_family_via_all_pages() when need to retry parsing the first paginiation page."""

    fam = Family(
        "testFam",
        "testClass",
        "testURL",
        failed_pages={"url.com": ["first_pagination", 1]},
    )

    def mock_parse_pagination_page(*args, **kwargs):
        return fam, None, None, None

    monkeypatch.setattr(scrape_all, "parse_first_paginiation_page", mock_parse_pagination_page)

    scrape_all.parse_family_via_all_pages(
        family=fam,
        cazy_home=cazy_home_url,
        taxonomy_filters=None,
        ec_filters=None,
        args=None,
        session="session",
    )


def test_all_parse_family_one_rescrape(cazy_home_url, monkeypatch):
    """Test parse_family_via_all_pages when there's one page to rescrape."""

    family = Family("testFam", "testClass", "testURL", failed_pages={"url.com": 1})

    def mock_parse_proteins(*args, **kwargs):
        return family, None, None, None

    monkeypatch.setattr(scrape_all, "parse_proteins", mock_parse_proteins)

    scrape_all.parse_family_via_all_pages(
        family=family,
        cazy_home=cazy_home_url,
        taxonomy_filters=None,
        ec_filters=None,
        args=None,
        session="session",
    )


def test_all_parse_family_multiple_rescrapes(cazy_home_url, monkeypatch):
    """Test parse_family_via_all_pages when there are multiple fams to rescrape."""

    family = Family(
        "testFam",
        "testClass",
        "testURL",
        failed_pages={
            "url.com": 1,
            "www.com": 2,
            "com.com": 3
        },
    )

    def mock_parse_proteins(*args, **kwargs):
        return family, None, None, None

    monkeypatch.setattr(scrape_all, "parse_proteins", mock_parse_proteins)

    scrape_all.parse_family_via_all_pages(
        family=family,
        cazy_home=cazy_home_url,
        taxonomy_filters=None,
        ec_filters=None,
        args=None,
        session="session",
    )


# test parse_first_paginiation_page


def test_all_parse_first_pag_incorrect_format(cazy_home_url):
    """Test parse_first_paginiation_page when the URL is incorrectly formatted."""

    fam = Family("famName", "CAZyClass", "url")

    scrape_all.parse_first_paginiation_page(
        first_pagination_url="url.url.url",
        family=fam,
        failed_scrapes=[],
        sql_failures=[],
        format_failures=[],
        cazy_home=cazy_home_url,
        args="args",
        taxonomy_filters=set(),
        ec_filters=[],
        session="session",
    )


def test_all_parse_first_pag_no_protein_total(fam_url, cazy_home_url, args, monkeypatch):
    """Test parse_first_paginiation_page when no number of total proteins is retrurned."""
    fam = Family("famName", "CAZyClass", fam_url)

    def mock_get_pag_data(*args, **kwargs):
        return None, "total"

    def mock_parse_error(*args, **kwargs):
        return fam, None, None, None

    monkeypatch.setattr(scrape_all, "get_paginiation_data", mock_get_pag_data)
    monkeypatch.setattr(scrape_all, "parse_total_proteins_error", mock_parse_error)

    scrape_all.parse_first_paginiation_page(
        first_pagination_url=fam_url,
        family=fam,
        failed_scrapes=[],
        sql_failures=[],
        format_failures=[],
        cazy_home=cazy_home_url,
        args=args["args"],
        taxonomy_filters=set(),
        ec_filters=[],
        session="session",
    )


def test_all_parse_first_pag_success(fam_url, cazy_home_url, args, monkeypatch):
    """Test parse_first_paginiation_page when successful."""
    fam = Family("famName", "CAZyClass", fam_url)

    def mock_get_pag_data(*args, **kwargs):
        return None, 1256

    def mock_parse_proteins(*args, **kwargs):
        return fam, None, None, None

    monkeypatch.setattr(scrape_all, "get_paginiation_data", mock_get_pag_data)
    monkeypatch.setattr(scrape_all, "parse_proteins", mock_parse_proteins)

    scrape_all.parse_first_paginiation_page(
        first_pagination_url=fam_url,
        family=fam,
        failed_scrapes=[],
        sql_failures=[],
        format_failures=[],
        cazy_home=cazy_home_url,
        args=args["args"],
        taxonomy_filters=set(),
        ec_filters=[],
        session="session",
    )


# test get_pagination_data


def test_all_get_pag_no_page(fam_url, fam_template, cazy_home_url, args, monkeypatch):
    """Test get_pagination_data when no page is retrieved from CAZy."""

    def mock_get_page(*args, **kwargs):
        return None, "error"

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)

    res1, res2 = scrape_all.get_paginiation_data(
        first_pagination_url=fam_url,
        family=fam_template,
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )
    # assert type(res1) is dict
    # assert res2 is None


def test_all_get_pag_deleted_fam(fam_url, fam_template, cazy_home_url, args, monkeypatch):
    """Test get_pagination_data when parsing a 'deteled' family."""

    def mock_get_page(*args, **kwargs):
        return "page", None

    def mock_get_urls(*args, **kwargs):
        return ['url'], 'Deleted family!'

    def mock_add_sql(*args, **kwargs):
        return

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)
    monkeypatch.setattr(scrape_all, "get_paginiation_page_urls", mock_get_urls)
    monkeypatch.setattr(sql_interface, "add_deleted_cazy_family", mock_add_sql)

    res1, res2 = scrape_all.get_paginiation_data(
        first_pagination_url=fam_url,
        family=fam_template,
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )

    # assert type(res1) is dict
    # assert res2 is None


def test_all_get_pag_empty_fam(fam_url, fam_template, cazy_home_url, args, monkeypatch):
    """Test get_pagination_data when parsing an 'empty' family."""

    def mock_get_page(*args, **kwargs):
        return "page", None

    def mock_get_urls(*args, **kwargs):
        return ['url'], 'Empty family'

    def mock_add_sql(*args, **kwargs):
        return

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)
    monkeypatch.setattr(scrape_all, "get_paginiation_page_urls", mock_get_urls)
    monkeypatch.setattr(sql_interface, "add_deleted_cazy_family", mock_add_sql)

    res1, res2 = scrape_all.get_paginiation_data(
        first_pagination_url=fam_url,
        family=fam_template,
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )

    # assert type(res1) is dict
    # assert res2 is None


def test_all_get_pag_failed(fam_url, fam_template, cazy_home_url, args, monkeypatch):
    """Test get_pagination_data when failed to retrieve the total protein count."""

    def mock_get_page(*args, **kwargs):
        return "page", None

    def mock_get_urls(*args, **kwargs):
        return ['url'], 'Failed Retrieval'

    def mock_add_sql(*args, **kwargs):
        return

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)
    monkeypatch.setattr(scrape_all, "get_paginiation_page_urls", mock_get_urls)
    monkeypatch.setattr(sql_interface, "add_deleted_cazy_family", mock_add_sql)

    res1, res2 = scrape_all.get_paginiation_data(
        first_pagination_url=fam_url,
        family=fam_template,
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )

    # assert type(res1) is dict
    # assert res2 is None


def test_all_get_pag_zero(fam_url, fam_template, cazy_home_url, args, monkeypatch):
    """Test get_pagination_data when no page urls are retrieved."""

    def mock_get_page(*args, **kwargs):
        return "page", None

    def mock_get_urls(*args, **kwargs):
        return [], 200

    def mock_add_sql(*args, **kwargs):
        return

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)
    monkeypatch.setattr(scrape_all, "get_paginiation_page_urls", mock_get_urls)
    monkeypatch.setattr(sql_interface, "add_deleted_cazy_family", mock_add_sql)

    res1, res2 = scrape_all.get_paginiation_data(
        first_pagination_url=fam_url,
        family=fam_template,
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )

    # assert type(res1) is dict
    # assert res2 is None


def test_all_get_pag_success(fam_url, fam_template, cazy_home_url, args, monkeypatch):
    """Test get_pagination_data when all is successful."""

    def mock_get_page(*args, **kwargs):
        return "page", None

    def mock_get_urls(*args, **kwargs):
        return ['url'], 1000

    def mock_add_sql(*args, **kwargs):
        return

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)
    monkeypatch.setattr(scrape_all, "get_paginiation_page_urls", mock_get_urls)
    monkeypatch.setattr(sql_interface, "add_deleted_cazy_family", mock_add_sql)

    res1, res2 = scrape_all.get_paginiation_data(
        first_pagination_url=fam_url,
        family=fam_template,
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )

    # assert type(res1) is list
    # assert type(res2) is int


# test get_pagination_page_urls()


def test_all_get_urls_no_pag(fam_url, cazy_home_url, input_dir):
    """Test get_pagination_page_urls when family is not paginated."""
    test_input_path = input_dir / "test_get_pagination_urls" / "no_pagination_pag.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    scrape_all.get_paginiation_page_urls(
        fam_url,
        soup,
        cazy_home_url,
        "fam_template",
    )


def test_all_get_urls_pag(fam_url, cazy_home_url, input_dir):
    """Test get_pagination_page_urls when the family has pagination."""
    test_input_path = input_dir / "test_get_pagination_urls" / "pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    scrape_all.get_paginiation_page_urls(
        fam_url,
        soup,
        cazy_home_url,
        "fam_template",
    )


def test_all_get_urls_deleted_fam(fam_url, cazy_home_url, input_dir):
    """test get_paginiation_page_urls when parsing a deleted family."""
    test_input_path = input_dir / "test_get_pagination_urls" / "deleted_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    scrape_all.get_paginiation_page_urls(
        fam_url,
        soup,
        cazy_home_url,
        "fam_template",
    )


def test_all_get_urls_empty_error(fam_url, cazy_home_url, input_dir):
    """test get_paginiation_page_urls when parsing an empty family."""
    test_input_path = input_dir / "test_get_pagination_urls" / "empty_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    scrape_all.get_paginiation_page_urls(
        fam_url,
        soup,
        cazy_home_url,
        "fam_template",
    )


def test_all_get_urls_double_error(fam_url, cazy_home_url, input_dir):
    """test get_paginiation_page_urls when cannot retrieve the total protein count."""
    test_input_path = input_dir / "test_get_pagination_urls" / "no_protein_total.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    scrape_all.get_paginiation_page_urls(
        fam_url,
        soup,
        cazy_home_url,
        "fam_template",
    )


# test parse_total_proteins_error()


def test_all_parse_total_proteins_error_first_pag(fam_url, fam_template, args):
    """Test parse_total_proteins_error() when it's the first failure for the pagination page."""

    errors = {
        "url": fam_url,
        "format": "error message",
    }

    family, failed_scrapes, sql_failures, format_failures = scrape_all.parse_total_proteins_error(
        first_pagination_url=fam_url,
        family=fam_template,
        failed_scrapes=[],
        sql_failures=[],
        format_failures=[],
        errors=errors,
        args=args["args"],
    )

    # assert format_failures == ['http://www.cazy.org/GH14_all.html\terror message']
    # assert failed_scrapes == []
    # assert sql_failures == []


def test_all_parse_total_proteins_error_final_pag(fam_url, fam_template, args):
    """Test parse_total_proteins_error() when it's NOT the first failure for the pagination page."""

    fam = fam_template
    fam.failed_pages = {fam_url: ["first_pagination", 2]}

    errors = {
        "url": fam_url,
        "format": "error message",
    }

    family, failed_scrapes, sql_failures, format_failures = scrape_all.parse_total_proteins_error(
        first_pagination_url=fam_url,
        family=fam,
        failed_scrapes=[],
        sql_failures=[],
        format_failures=[],
        errors=errors,
        args=args["args"],
    )

    # assert format_failures == ['http://www.cazy.org/GH14_all.html\terror message']
    # assert failed_scrapes == [
        (
            'http://www.cazy.org/GH14_all.html\tCAZyClass\tFailed to connect to this '
            'page of proteins for famName, and raised the following error message:\n'
            'http://www.cazy.org/GH14_all.html'
        ),
    ]
    # assert sql_failures == []


# test parse_proteins()


def test_all_parse_proteins_s_all(protein_gen, fam_template, args, monkeypatch):
    """Test parse_proteins from scrape_all.py in the crawler module."""
    fam = fam_template
    fam.failed_pages = {"www.cazy.org/GH3.html": 3}

    def mock_parse_protein_tables(*args, **kwargs):
        return protein_gen

    monkeypatch.setattr(scrape_all, "parse_protein_table", mock_parse_protein_tables)

    family, failed_scrapes, sql_failures, format_failures = scrape_all.parse_proteins(
        protein_page_urls=[1, 2, 3, 4],
        total_proteins=100,
        family=fam,
        taxonomy_filters=set(),
        ec_filters=set(),
        session="session",
        args=args["args"],
        failed_scrapes=[],
        sql_failures=[],
        format_failures=[],
    )

    # assert failed_scrapes == ['w\tCAZyClass\tFailed to connect to this page of proteins for famName, and raised the following error message:\nw']
    # assert sql_failures == ['protein_name was not added to the database\tand raised the following error when atempting to do so:\nerror message', 'protein_name was not added to the database\tand raised the following error when atempting to do so:\nerror message']
    # assert format_failures == ['no internet connection was inconsistently formated\traised the following error:\nerror message', 'sql error was inconsistently formated\traised the following error:\nerror message', 'sql error was inconsistently formated\traised the following error:\nerror message']


# test parse_protein_table()


def test_all_parse_table_no_page(fam_url, args, monkeypatch):
    """Test parse_protein_table when no page is retrieved from CAZy."""

    def mock_get_page(*args, **kwargs):
        return None, "error message"

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)

    scrape_all.parse_protein_table(
        protein_page_url=fam_url,
        family_name="testFam",
        taxonomy_filters=set(),
        ec_filters=set(),
        session="session",
        args=args["args"],
    )


def test_all_parse_table_deleted_fam(fam_url, args, input_dir, monkeypatch):
    """Test parse_protein_table when it is a deleted family."""
    test_input_path = input_dir / "test_get_pagination_urls" / "deleted_fam.html"
    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)

    scrape_all.parse_protein_table(
        protein_page_url=fam_url,
        family_name="testFam",
        taxonomy_filters=set(),
        ec_filters=set(),
        session="session",
        args=args["args"],
    )


def test_all_parse_table_empty_fam(fam_url, args, input_dir, monkeypatch):
    """Test parse_protein_table when it is an empty family."""
    test_input_path = input_dir / "test_get_pagination_urls" / "empty_fam.html"
    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)

    scrape_all.parse_protein_table(
        protein_page_url=fam_url,
        family_name="testFam",
        taxonomy_filters=set(),
        ec_filters=set(),
        session="session",
        args=args["args"],
    )


def test_all_parse_table_no_protein_total(fam_url, args, input_dir, monkeypatch):
    """Test parse_protein_table when no protein total can be retrieved to check if an empty fam."""
    test_input_path = input_dir / "test_get_pagination_urls" / "no_protein_total.html"
    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)

    scrape_all.parse_protein_table(
        protein_page_url=fam_url,
        family_name="testFam",
        taxonomy_filters=set(),
        ec_filters=set(),
        session="session",
        args=args["args"],
    )


def test_all_parse_table_no_table(fam_url, args, input_dir, monkeypatch):
    """Test parse_protein_table when no CAZyme table is retrieved and there should be one."""

    test_input_path = input_dir / "test_parsing_proteins" / "no_cazyme_table.html"
    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)

    scrape_all.parse_protein_table(
        protein_page_url=fam_url,
        family_name="testFam",
        taxonomy_filters=set(),
        ec_filters=set(),
        session="session",
        args=args["args"],
    )


def test_all_parse_table_successful(fam_url, args, input_dir, protein_gen, monkeypatch):
    """Test parse_protein_table when all is successful."""

    test_input_path = input_dir / "test_get_pagination_urls" / "no_pagination_pag.html"
    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    def mock_row_to_protein(*args, **kwargs):
        return protein_gen

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "row_to_protein", mock_row_to_protein)

    scrape_all.parse_protein_table(
        protein_page_url=fam_url,
        family_name="testFam",
        taxonomy_filters=set(),
        ec_filters=set(),
        session="session",
        args=args["args"],
    )


# test parse_protein_table_dict()


def test_all_parse_page_dict_no_page(monkeypatch, fam_template, args):
    """Test scrape_all.parse_protein_table_dict() when no page is returned"""

    def mock_get_page(*args, **kwargs):
        return None, "error"

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)

    scrape_all.parse_protein_table_dict(
        family=fam_template,
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={},
        args=args["args"],
    )


def test_all_parse_page_dict_no_table(monkeypatch, fam_template, args, input_dir):
    """Test scrape_all.parse_protein_table_dict() when there is no CAZyme table"""
    test_input_path = input_dir / "test_parsing_proteins" / "no_table.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)

    scrape_all.parse_protein_table_dict(
        family=fam_template,
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={}, 
        args=args["args"],
    )


def test_all_parse_page_dict_deleted(monkeypatch, fam_template, args, input_dir):
    """Test scrape_all.parse_protein_table_dict() when it is a deleted family"""
    test_input_path = input_dir / "test_get_pagination_urls" / "deleted_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)

    scrape_all.parse_protein_table_dict(
        family=fam_template,
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={}, 
        args=args["args"],
    )


def test_all_parse_page_dict_empty(monkeypatch, fam_template, args, input_dir):
    """Test scrape_all.parse_protein_table_dict() when it is a empty family"""
    test_input_path = input_dir / "test_get_pagination_urls" / "empty_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)

    scrape_all.parse_protein_table_dict(
        family=fam_template,
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={}, 
        args=args["args"],
    )


def test_all_parse_page_dict_empty_table(monkeypatch, fam_template, args, input_dir):
    """Test scrape_all.parse_protein_table_dict() when there the CAZyme table is not populated with
    entries and it should be"""
    test_input_path = input_dir / "test_parsing_proteins" / "no_cazyme_table.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)

    scrape_all.parse_protein_table_dict(
        family=fam_template,
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={}, 
        args=args["args"],
    )


def test_all_parse_page_dict_success(monkeypatch, fam_template, args, input_dir):
    """Test scrape_all.parse_protein_table_dict() when all is successful"""
    test_input_path = input_dir / "test_get_pagination_urls" / "pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    def mock_row_to_dict(*args, **kwargs):
        return {"error": "error message", "protein": "protein name and or accession"}, {}

    monkeypatch.setattr(scrape_all, "get_page", mock_get_page)
    monkeypatch.setattr(scrape_all, "row_to_protein_in_dict", mock_row_to_dict)

    scrape_all.parse_protein_table_dict(
        family=fam_template,
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={}, 
        args=args["args"],
    )