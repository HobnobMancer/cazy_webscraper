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
from scraper.crawler.parse_cazy_families import scrape_by_kingdom
from scraper.sql import sql_interface


@pytest.fixture
def fam_url():
    url = "http://www.cazy.org/GH1_bacteria.html"
    return url


@pytest.fixture
def fam_template():
    fam = Family("famName", "CAZyClass", fam_url)
    return fam


@pytest.fixture
def args():
    dict_ = {'args': Namespace(retries=2, timeout=10)}
    return dict_


@pytest.fixture
def input_dir(test_dir):
    path_ = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_by_kingdom_inputs"
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



# test parse_family_by_kingdom()


def test_kgdm_parse_fam_first(cazy_home_url, monkeypatch):
    """Test parse_family_by_kingdom when parsing fam for the first time."""
    fam = Family("testFam", "testClass", "url.url")

    def mock_parse_pagination(*args, **kwargs):
        return fam, [], [], []

    monkeypatch.setattr(scrape_by_kingdom, "parse_kingdom_first_pagination_page", mock_parse_pagination)

    scrape_by_kingdom.parse_family_by_kingdom(
        family=fam,
        cazy_home=cazy_home_url,
        taxonomy_filters=set(),
        kingdoms=["Bacteria"],
        ec_filters=set(),
        args="args",
        session="session",
    )


def test_kgdm_parse_fam_reparse_pag(cazy_home_url, monkeypatch):
    """Test parse_family_by_kingdom when need to reparse the pagination page."""
    fam = Family(
        "testFam",
        "testClass",
        "url.url",
        {"Bacteria": {"url.com": ["first_paginiation", 1]}},
    )

    def mock_parse_pagination(*args, **kwargs):
        return fam, [], [], []

    monkeypatch.setattr(scrape_by_kingdom, "parse_kingdom_first_pagination_page", mock_parse_pagination)

    scrape_by_kingdom.parse_family_by_kingdom(
        family=fam,
        cazy_home=cazy_home_url,
        taxonomy_filters=set(),
        kingdoms=["Bacteria"],
        ec_filters=set(),
        args="args",
        session="session",
    )


def test_kgdm_parse_fam_reparse_pages(cazy_home_url, monkeypatch):
    """Test parse_family_by_kingdom when there's multiple pages to rescrape."""
    fam = Family(
        "testFam",
        "testClass",
        "url.url",
        {
            "Bacteria": {
                "url.com": 1,
                "urltest.com": 1,
                "test.com": 1,
            },
        },
    )

    def mock_parse_pagination(*args, **kwargs):
        return fam, [], []

    monkeypatch.setattr(scrape_by_kingdom, "parse_protein_pages", mock_parse_pagination)

    scrape_by_kingdom.parse_family_by_kingdom(
        family=fam,
        cazy_home=cazy_home_url,
        taxonomy_filters=set(),
        kingdoms=["Bacteria"],
        ec_filters=set(),
        args="args",
        session="session",
    )


# test parse_kingdom_first_pagination_page


def test_kgdm_first_pag_format(fam_template, cazy_home_url):
    """Test parse_kingdom_first_pagination_page when the URL is incorrectly formatted."""

    scrape_by_kingdom.parse_kingdom_first_pagination_page(
        first_pagination_url="URL",
        family=fam_template,
        kingdom="Bacteria",
        failed_scrapes=[],
        sql_failures=[],
        format_failures=[],
        cazy_home=cazy_home_url,
        taxonomy_filters=set(),
        ec_filters=set(),
        args="args",
        session="session",
    )


def test_kgdm_first_pag_no_protein_total(fam_template, fam_url, cazy_home_url, monkeypatch):
    """Test parse_kingdom_first_pagination_page when no protein total was retrieved"""

    def mock_parse_pag(*args, **kwargs):
        return [], None

    def mock_parse_errors(*args, **kwargs):
        return fam_template, [], [], []

    monkeypatch.setattr(scrape_by_kingdom, "get_kingdom_pagination_data", mock_parse_pag)
    monkeypatch.setattr(scrape_by_kingdom, "parse_kingdom_total_proteins_error", mock_parse_errors)

    scrape_by_kingdom.parse_kingdom_first_pagination_page(
        first_pagination_url=fam_url,
        family=fam_template,
        kingdom="Bacteria",
        failed_scrapes=[],
        sql_failures=[],
        format_failures=[],
        cazy_home=cazy_home_url,
        taxonomy_filters=set(),
        ec_filters=set(),
        args="args",
        session="session",
    )


def test_kgdm_first_pag_successful(fam_template, fam_url, cazy_home_url, monkeypatch):
    """Test parse_kingdom_first_pagination_page when all is successful."""

    def mock_parse_pag(*args, **kwargs):
        return ['url'], 100

    def mock_parse_tables(*args, **kwargs):
        return fam_template, [], [], []

    monkeypatch.setattr(scrape_by_kingdom, "get_kingdom_pagination_data", mock_parse_pag)
    monkeypatch.setattr(scrape_by_kingdom, "parse_protein_pages", mock_parse_tables)

    scrape_by_kingdom.parse_kingdom_first_pagination_page(
        first_pagination_url=fam_url,
        family=fam_template,
        kingdom="Bacteria",
        failed_scrapes=[],
        sql_failures=[],
        format_failures=[],
        cazy_home=cazy_home_url,
        taxonomy_filters=set(),
        ec_filters=set(),
        args="args",
        session="session",
    )


# test get_kingdom_pagination_data()


def test_kgdm_get_kingdom_pag_no_page(fam_url, fam_template, cazy_home_url, args, monkeypatch):
    """Test get_kingdom_pagination_data when no page is retrieved from CAZy."""

    def mock_get_page(*args, **kwargs):
        return None, "error message"

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.get_kingdom_pagination_data(
        first_pagination_url=fam_url,
        family=fam_template,
        kingdom="Bacteria",
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )


def test_kgdm_get_kingdom_pag_deteled(fam_url, fam_template, cazy_home_url, args, input_dir, monkeypatch):
    """Test get_kingdom_pagination_data when the family is a deleted family."""

    test_input_path = input_dir / "kb_pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, "error message"

    def mock_get_urls(*args, **kwargs):
        return ['url'], "Deleted family!"

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)
    monkeypatch.setattr(scrape_by_kingdom, "get_kingdom_page_urls", mock_get_urls)
    monkeypatch.setattr(sql_interface, "add_deleted_cazy_family", mock_sql)

    scrape_by_kingdom.get_kingdom_pagination_data(
        first_pagination_url=fam_url,
        family=fam_template,
        kingdom="Bacteria",
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )


def test_kgdm_get_kingdom_pag_empty(fam_url, fam_template, cazy_home_url, args, input_dir, monkeypatch):
    """Test get_kingdom_pagination_data when the family is an empty family."""

    test_input_path = input_dir / "kb_pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, "error message"

    def mock_get_urls(*args, **kwargs):
        return ['url'], "Empty family"

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)
    monkeypatch.setattr(scrape_by_kingdom, "get_kingdom_page_urls", mock_get_urls)
    monkeypatch.setattr(sql_interface, "add_deleted_cazy_family", mock_sql)

    scrape_by_kingdom.get_kingdom_pagination_data(
        first_pagination_url=fam_url,
        family=fam_template,
        kingdom="Bacteria",
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )


def test_kgdm_get_kingdom_pag_failed(fam_url, fam_template, cazy_home_url, args, input_dir, monkeypatch):
    """Test get_kingdom_pagination_data when pagination data is not retrieved."""

    test_input_path = input_dir / "kb_pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, "error message"

    def mock_get_urls(*args, **kwargs):
        return ['url'], 'Failed Retrieval'

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)
    monkeypatch.setattr(scrape_by_kingdom, "get_kingdom_page_urls", mock_get_urls)
    monkeypatch.setattr(sql_interface, "add_deleted_cazy_family", mock_sql)

    scrape_by_kingdom.get_kingdom_pagination_data(
        first_pagination_url=fam_url,
        family=fam_template,
        kingdom="Bacteria",
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )


def test_kgdm_get_kingdom_pag_no_urls(fam_url, fam_template, cazy_home_url, args, input_dir, monkeypatch):
    """Test get_kingdom_pagination_data when no urls are retrieved."""

    test_input_path = input_dir / "kb_pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, "error message"

    def mock_get_urls(*args, **kwargs):
        return [], 100

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)
    monkeypatch.setattr(scrape_by_kingdom, "get_kingdom_page_urls", mock_get_urls)
    monkeypatch.setattr(sql_interface, "add_deleted_cazy_family", mock_sql)

    scrape_by_kingdom.get_kingdom_pagination_data(
        first_pagination_url=fam_url,
        family=fam_template,
        kingdom="Bacteria",
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )


def test_kgdm_get_kingdom_pag_successful(fam_url, fam_template, cazy_home_url, args, input_dir, monkeypatch):
    """Test get_kingdom_pagination_data when nall is successful."""

    test_input_path = input_dir / "kb_pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, "error message"

    def mock_get_urls(*args, **kwargs):
        return ['url'], 100

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)
    monkeypatch.setattr(scrape_by_kingdom, "get_kingdom_page_urls", mock_get_urls)
    monkeypatch.setattr(sql_interface, "add_deleted_cazy_family", mock_sql)

    scrape_by_kingdom.get_kingdom_pagination_data(
        first_pagination_url=fam_url,
        family=fam_template,
        kingdom="Bacteria",
        cazy_home=cazy_home_url,
        args=args["args"],
        session="session",
    )


# test get_kingdom_page_urls()


def test_kgdm_get_urls_no_pag(fam_url, input_dir, cazy_home_url):
    """Test get_kingdom_page_urls when fam has no pagination."""

    test_input_path = input_dir / "kb_no_pag.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    scrape_by_kingdom.get_kingdom_page_urls(
        first_pagination_url=fam_url,
        first_pagination_page=soup,
        kingdom="Bacteria",
        cazy_home=cazy_home_url,
        family_name="testFam",
    )


def test_kgdm_get_urls_pag(fam_url, input_dir, cazy_home_url):
    """Test get_kingdom_page_urls when fam has has pagination."""

    test_input_path = input_dir / "kb_pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    scrape_by_kingdom.get_kingdom_page_urls(
        first_pagination_url=fam_url,
        first_pagination_page=soup,
        kingdom="Bacteria",
        cazy_home=cazy_home_url,
        family_name="testFam",
    )


def test_kgdm_get_urls_deleted(fam_url, test_dir, cazy_home_url):
    """Test get_kingdom_page_urls when fam is a deleted fam."""

    test_input_path = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_all_inputs" / "test_get_pagination_urls" / "deleted_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    scrape_by_kingdom.get_kingdom_page_urls(
        first_pagination_url=fam_url,
        first_pagination_page=soup,
        kingdom="Bacteria",
        cazy_home=cazy_home_url,
        family_name="testFam",
    )


def test_kgdm_get_urls_empty(fam_url, test_dir, cazy_home_url):
    """Test get_kingdom_page_urls when fam is an empty fam."""

    test_input_path = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_all_inputs" / "test_get_pagination_urls" / "empty_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    scrape_by_kingdom.get_kingdom_page_urls(
        first_pagination_url=fam_url,
        first_pagination_page=soup,
        kingdom="Bacteria",
        cazy_home=cazy_home_url,
        family_name="testFam",
    )


def test_kgdm_get_urls_no_count(fam_url, test_dir, cazy_home_url):
    """Test get_kingdom_page_urls when the protein total cannot be retrieved."""

    test_input_path = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_all_inputs" / "test_get_pagination_urls" / "no_protein_total.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    scrape_by_kingdom.get_kingdom_page_urls(
        first_pagination_url=fam_url,
        first_pagination_page=soup,
        kingdom="Bacteria",
        cazy_home=cazy_home_url,
        family_name="testFam",
    )


# test parse_kingdom_total_proteins_error()


def test_kgdm_parse_total_proteins_error_first_pag(fam_url, fam_template, args):
    """Test parse_total_proteins_error() when it's the first failure for the pagination page."""

    errors = {
        "url": fam_url,
        "format": "error message",
    }

    family, failed_scrapes, sql_failures, format_failures = scrape_by_kingdom.parse_kingdom_total_proteins_error(
        first_pagination_url=fam_url,
        family=fam_template,
        kingdom="Bacteria",
        failed_scrapes=[],
        sql_failures=[],
        format_failures=[],
        errors=errors,
        args=args["args"],
    )

    # assert format_failures == ['http://www.cazy.org/GH1_bacteria.html\terror message']
    # assert failed_scrapes == []
    # assert sql_failures == []


def test_kgdm_parse_total_proteins_error_final_pag(fam_url, fam_template, args):
    """Test parse_total_proteins_error() when it's NOT the first failure for the pagination page."""

    fam = fam_template
    fam.failed_pages = {"Bacteria": {fam_url: ["first_pagination", 2]}}

    errors = {
        "url": fam_url,
        "format": "error message",
    }

    family, failed_scrapes, sql_failures, format_failures = scrape_by_kingdom.parse_kingdom_total_proteins_error(
        first_pagination_url=fam_url,
        family=fam_template,
        kingdom="Bacteria",
        failed_scrapes=[],
        sql_failures=[],
        format_failures=[],
        errors=errors,
        args=args["args"],
    )

    # assert format_failures == ['http://www.cazy.org/GH1_bacteria.html\terror message']
    # assert failed_scrapes == ['http://www.cazy.org/GH1_bacteria.html\tCAZyClass\tBacteria\tfamName\tFailed to connect to this page of proteins for famName, and raised the following error message:\nhttp://www.cazy.org/GH1_bacteria.html']
    # assert sql_failures == []


# test parse_protein_pages()


def test_kgdm_all_parse_proteins(protein_gen, fam_template, args, monkeypatch):
    """Test parse_proteins from scrape_all.py in the crawler module."""
    fam = fam_template
    fam.failed_pages = {"Bacteria": {"www.cazy.org/GH3.html": 3}}

    def mock_parse_protein_pages(*args, **kwargs):
        return protein_gen

    monkeypatch.setattr(scrape_by_kingdom, "parse_kingdom_protein_tables", mock_parse_protein_pages)

    family, failed_scrapes, sql_failures = scrape_by_kingdom.parse_protein_pages(
        family=fam,
        kingdom="Bacteria",
        protein_page_urls=[1, 2, 3, 4],
        total_proteins=100,
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        sql_failures=[],
        session="session",
        args=args["args"],
    )

    # assert failed_scrapes == ['www.cazy.org/GH3.html\tCAZyClass\tFailed to connect to this page of proteins for famName\terror message']
    # assert sql_failures == ['protein_name was not added to the database, and raised the following error when atempting to do so:\nerror message', 'protein_name was not added to the database, and raised the following error when atempting to do so:\nerror message']


# test parse_kingdom_protein_tables()


def test_kgdm_parse_proteins_no_page(args, monkeypatch):
    """Test parse_kingdom_protein_tables when no page is returned."""
    def mock_get_page(*args, **kwargs):
        return None, "error"

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.parse_kingdom_protein_tables(
        protein_page_url=fam_url,
        family_name="testFam",
        kingdom="Bacteria",
        taxonomy_filters=set(),
        ec_filters=set(),
        args=args["args"],
        session="session",
    )


def test_kgdm_parse_proteins_deleted(test_dir, args, monkeypatch):
    """Test parse_kingdom_protein_tables when it is a deleted family."""
    test_input_path = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_all_inputs" / "test_get_pagination_urls" / "deleted_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.parse_kingdom_protein_tables(
        protein_page_url=fam_url,
        family_name="testFam",
        kingdom="Bacteria",
        taxonomy_filters=set(),
        ec_filters=set(),
        args=args["args"],
        session="session",
    )


def test_kgdm_parse_proteins_empty(test_dir, args, monkeypatch):
    """Test parse_kingdom_protein_tables when it is an empty family."""
    test_input_path = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_all_inputs" / "test_get_pagination_urls" / "empty_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.parse_kingdom_protein_tables(
        protein_page_url=fam_url,
        family_name="testFam",
        kingdom="Bacteria",
        taxonomy_filters=set(),
        ec_filters=set(),
        args=args["args"],
        session="session",
    )


def test_kgdm_parse_proteins_no_table(input_dir, args, monkeypatch):
    """Test parse_kingdom_protein_tables when there is no protein table."""
    test_input_path = input_dir / "no_table.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.parse_kingdom_protein_tables(
        protein_page_url=fam_url,
        family_name="testFam",
        kingdom="Bacteria",
        taxonomy_filters=set(),
        ec_filters=set(),
        args=args["args"],
        session="session",
    )


def test_kgdm_parse_proteins_empty_table(input_dir, args, monkeypatch):
    """Test parse_kingdom_protein_tables when the CAZyme table is not populated by should be."""
    test_input_path = input_dir / "no_table.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.parse_kingdom_protein_tables(
        protein_page_url=fam_url,
        family_name="testFam",
        kingdom="Bacteria",
        taxonomy_filters=set(),
        ec_filters=set(),
        args=args["args"],
        session="session",
    )


def test_kgdm_parse_proteins_successul(input_dir, args, monkeypatch):
    """Test parse_kingdom_protein_tables when the CAZyme table is not populated by should be."""
    test_input_path = input_dir / "kb_pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.parse_kingdom_protein_tables(
        protein_page_url=fam_url,
        family_name="testFam",
        kingdom="Bacteria",
        taxonomy_filters=set(),
        ec_filters=set(),
        args=args["args"],
        session="session",
    )


# test parse_protein_pages_dict()


def test_kngm_parse_page_dict_no_page(monkeypatch, fam_template, args):
    """Test parse_protein_pages_dict() when no page is returned"""

    def mock_get_page(*args, **kwargs):
        return None, "error"

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.parse_protein_pages_dict(
        family=fam_template,
        kingdom="Bacteria",
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={}, 
        args=args["args"],
    )


def test_kngm_parse_page_dict_no_table(monkeypatch, fam_template, args, input_dir):
    """Test parse_protein_pages_dict() when there is no CAZyme table"""
    test_input_path = input_dir / "no_table.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.parse_protein_pages_dict(
        family=fam_template,
        kingdom="Bacteria",
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={}, 
        args=args["args"],
    )


def test_kngm_parse_page_dict_deleted(monkeypatch, fam_template, args, test_dir):
    """Test parse_protein_pages_dict() when it is a deleted family"""
    test_input_path = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_all_inputs" / "test_get_pagination_urls" / "deleted_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.parse_protein_pages_dict(
        family=fam_template,
        kingdom="Bacteria",
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={}, 
        args=args["args"],
    )


def test_kngm_parse_page_dict_empty(monkeypatch, fam_template, args, test_dir):
    """Test parse_protein_pages_dict() when it is an empty family"""
    test_input_path = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_all_inputs" / "test_get_pagination_urls" / "empty_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.parse_protein_pages_dict(
        family=fam_template,
        kingdom="Bacteria",
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={}, 
        args=args["args"],
    )


def test_kngm_parse_page_dict_continued_incomplete(monkeypatch, fam_template, args, input_dir):
    """Test parse_protein_pages_dict() when the page is incompletely returned"""
    test_input_path = input_dir / "kb_incomplete.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page)

    scrape_by_kingdom.parse_protein_pages_dict(
        family=fam_template,
        kingdom="Bacteria",
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={}, 
        args=args["args"],
    )


def test_kngm_parse_page_dict_success(monkeypatch, fam_template, args, input_dir):
    """Test parse_protein_pages_dict() when all is successful"""
    test_input_path = input_dir / "kb_pag_page.html"

    output_1 = {"error": "error message", "protein": "protein name and or accession"}
    output_2 = {}

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page_(*args, **kwargs):
        return soup, "error"

    def mock_row_to_protein(*args, **kwargs):
        return output_1, output_2

    monkeypatch.setattr(scrape_by_kingdom, "get_page", mock_get_page_)
    monkeypatch.setattr(scrape_by_kingdom, "row_to_protein_in_dict", mock_row_to_protein)

    scrape_by_kingdom.parse_protein_pages_dict(
        family=fam_template,
        kingdom="Bacteria",
        protein_page_urls=["test_url"],
        taxonomy_filters=set(),
        ec_filters=set(),
        failed_scrapes=[],
        session={},
        args=args["args"],
    )

# test parse_url_error_kngdm()

