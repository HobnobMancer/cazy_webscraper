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
"""Tests the script get_cazy_pages which creates a local library of HTML pages from CAZy.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import pytest

import pandas as pd

from argparse import Namespace
from datetime import datetime

from bs4 import BeautifulSoup

from scraper import cazy_webscraper, crawler, sql, utilities
from scraper import crawler
from scraper.crawler import Family
from scraper.crawler.cazy_html_pages import get_cazy_pages
from scraper.utilities import file_io, parse_configuration, parsers


@pytest.fixture
def input_dir(test_dir):
    path_ = test_dir / "test_inputs" / "test_inputs_crawler" / "get_cazy_pages"
    return path_


@pytest.fixture
def output_dir(test_dir):
    path_ = test_dir / "test_outputs"
    return path_


@pytest.fixture
def html_dir(output_dir):
    path_ = output_dir / "test_outputs_get_cazy_pages" / "html_pages"
    return path_


@pytest.fixture
def logs_dir(output_dir):
    path_ = output_dir / "test_outputs_get_cazy_pages" / "test_logs"
    return path_


@pytest.fixture
def input_pages(test_dir):
    path_ = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_all_inputs" / "test_get_pagination_urls"
    return path_


@pytest.fixture
def pag_dir(test_dir):
    path_ = test_dir / "test_inputs" / "test_inputs_crawler" / "scrape_by_kingdom_inputs"
    return path_


@pytest.fixture
def args(logs_dir):
    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
            output=logs_dir,
        )
    }
    return argsdict


@pytest.fixture
def args_pages(html_dir):
    argsdict = {
        "args": Namespace(
            subfamilies=True,
            retries=2,
            timeout=5,
            output=html_dir,
        )
    }
    return argsdict


@pytest.fixture
def start_time():
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)
    return start_time


@pytest.fixture
def protein_gen():
    result_list = [
        {"url": None, "format": None, "sql": None, "error": None},
        {"url": "www.cazy.org/GH1.html", "error": "error message"},
        {"url": None, "error": "error message"},
        {"url": "www.cazy.org/GH3.html", "error": "error message"},
    ]
    return (_ for _ in result_list)


@pytest.fixture
def protein_gen_success():
    result_list = [
        {"url": None, "format": None, "sql": None, "error": None},
        {"url": None, "format": None, "sql": None, "error": None},
        {"url": None, "format": None, "sql": None, "error": None},
    ]
    return (_ for _ in result_list)


@pytest.fixture
def fam_url():
    url = "http://www.cazy.org/GH14_all.html"
    return url


@pytest.fixture
def fam_template():
    fam = Family("famName", "CAZyClass", fam_url)
    return fam


# test get_cazy_pages


def test_pages_first_class_parse_no_fam_urls(
    cazy_dictionary, cazy_home_url, logs_dir, args, start_time, monkeypatch,
):
    """Test get_cazy_pages when parsing the CAZy class for the first time 
    and no fam urls are retrieved."""

    file_io.make_output_directory(logs_dir, True, False)

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass("test_class", "test_class_url.html", 0)
        return [class1]

    def mock_get_families(*args, **kwargs):
        return None, "test error message", ["test_url1", "test_url2"]

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_families)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_cazy_pages.get_cazy_pages(
            args=args["args"],
            cazy_home=cazy_home_url,
            time_stamp="time_stamp",
            excluded_classes=None,
            cazy_dict=cazy_dictionary,
            config_dict=None,
            kingdoms=None,
            start_time=start_time,
            )
    assert pytest_wrapped_e.type == SystemExit

    file_io.make_output_directory(logs_dir, True, False)


def test_pages_reparse_class_parse_not_kingdom_no_config(
    cazy_dictionary, cazy_home_url, logs_dir, args, start_time, monkeypatch,
):
    """Test get_cazy_pages when reparsing the CAZy class"""

    file_io.make_output_directory(logs_dir, True, False)

    fam1 = crawler.Family("test_fam", "test_class", "test_url")

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass(
            name="test_class",
            url="test_class_url.html",
            tries=0,
            failed_families={fam1: 0}
        )
        return [class1]

    def mock_get_families(*args, **kwargs):
        return [fam1], "error message", ["in", "cor", "rect", "urls"]

    def mock_parse_family(*args, **kwargs):
        return fam1, True, ["fail1", "fail2"], ["format error"]

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_families)
    monkeypatch.setattr(get_cazy_pages, "parse_all_family", mock_parse_family)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_cazy_pages.get_cazy_pages(
            args=args["args"],
            cazy_home=cazy_home_url,
            time_stamp="time_stamp",
            excluded_classes=None,
            cazy_dict=cazy_dictionary,
            config_dict=None,
            kingdoms='all',
            start_time=start_time,
            )
    assert pytest_wrapped_e.type == SystemExit

    file_io.make_output_directory(logs_dir, True, False)


def test_pages_reparse_class_parse_kingdom_no_config(
    cazy_dictionary, cazy_home_url, logs_dir, args, start_time, monkeypatch,
):
    """Test get_cazy_pages when reparsing the CAZy class"""

    file_io.make_output_directory(logs_dir, True, False)

    fam1 = crawler.Family("test_fam", "test_class", "test_url")

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass(
            name="test_class",
            url="test_class_url.html",
            tries=0,
            failed_families={fam1: 0}
        )
        return [class1]

    def mock_get_families(*args, **kwargs):
        return [fam1], "error message", ["in", "cor", "rect", "urls"]

    def mock_parse_family(*args, **kwargs):
        return fam1, True, ["fail1", "fail2"], ["format error"]

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_families)
    monkeypatch.setattr(get_cazy_pages, "parse_family_by_kingdom", mock_parse_family)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_cazy_pages.get_cazy_pages(
            args=args["args"],
            cazy_home=cazy_home_url,
            time_stamp="time_stamp",
            excluded_classes=None,
            cazy_dict=cazy_dictionary,
            config_dict=None,
            kingdoms=["Bacteria"],
            start_time=start_time,
            )
    assert pytest_wrapped_e.type == SystemExit

    file_io.make_output_directory(logs_dir, True, False)


# def test_pages_reparse_class_parse_not_kingdom_config(
#     cazy_dictionary, cazy_home_url, logs_dir, args, start_time, monkeypatch,
# ):
#     """Test get_cazy_pages when reparsing the CAZy class"""

#     file_io.make_output_directory(logs_dir, True, False)

#     fam1 = crawler.Family("GH3_1", "test_class", "test_url")

#     config_dict = {"Glycoside Hydrolases": ["GH3"]}

#     def mock_get_classes(*args, **kwargs):
#         class1 = crawler.CazyClass(
#             name="Glycoside Hydrolases",
#             url="test_class_url.html",
#             tries=0,
#             failed_families={fam1: 0}
#         )
#         return [class1]

#     def mock_get_families(*args, **kwargs):
#         return [fam1], "error message", ["in", "cor", "rect", "urls"]

#     def mock_parse_family(*args, **kwargs):
#         return fam1, True, ["fail1", "fail2"], ["format error"]

#     monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
#     monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_families)
#     monkeypatch.setattr(get_cazy_pages, "parse_family_by_kingdom", mock_parse_family)

#     with pytest.raises(SystemExit) as pytest_wrapped_e:
#         get_cazy_pages.get_cazy_pages(
#             args=args["args"],
#             cazy_home=cazy_home_url,
#             time_stamp="time_stamp",
#             excluded_classes=None,
#             cazy_dict=cazy_dictionary,
#             config_dict=config_dict,
#             kingdoms="all",
#             start_time=start_time,
#             )
#     assert pytest_wrapped_e.type == SystemExit

#     file_io.make_output_directory(logs_dir, True, False)


def test_pages_reparse_class_parse_kingdom_config(
    cazy_dictionary, cazy_home_url, logs_dir, start_time, monkeypatch,
):
    """Test get_cazy_pages when reparsing the CAZy class"""

    args = {
        "args": Namespace(
            subfamilies=False,
            retries=2,
            timeout=5,
            output=logs_dir,
        )
    }

    file_io.make_output_directory(logs_dir, True, False)

    fam1 = crawler.Family("GH3_1", "test_class", "test_url")

    config_dict = {"Glycoside Hydrolases": ["GH3"]}

    def mock_get_classes(*args, **kwargs):
        class1 = crawler.CazyClass(
            name="Glycoside Hydrolases",
            url="test_class_url.html",
            tries=0,
            failed_families={fam1: 0}
        )
        return [class1]

    def mock_get_families(*args, **kwargs):
        return [fam1], "error message", ["in", "cor", "rect", "urls"]

    def mock_parse_family(*args, **kwargs):
        return fam1, True, ["fail1", "fail2"], ["format error"]

    monkeypatch.setattr(crawler, "get_cazy_classes", mock_get_classes)
    monkeypatch.setattr(crawler, "get_cazy_family_urls", mock_get_families)
    monkeypatch.setattr(get_cazy_pages, "parse_family_by_kingdom", mock_parse_family)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_cazy_pages.get_cazy_pages(
            args=args["args"],
            cazy_home=cazy_home_url,
            time_stamp="time_stamp",
            excluded_classes=None,
            cazy_dict=cazy_dictionary,
            config_dict=config_dict,
            kingdoms=["Bacteria"],
            start_time=start_time,
            )
    assert pytest_wrapped_e.type == SystemExit

    file_io.make_output_directory(logs_dir, True, False)


# test parse_all_family()


def test_pages_parse_all_rescrape_first_pag(args, protein_gen, cazy_home_url, monkeypatch):
    """Test parse_all_family when rescraping a family and 
    need to rescrape the first paginiation page."""

    fam = Family(
        "testFam",
        "testClass",
        "testURL",
        failed_pages={"url.com": ["first_pagination", 1]},
    )

    def mock_parse_pagination_page(*args, **kwargs):
        return ['urls', 'urls']

    def mock_get_pages(*args, **kwargs):
        return protein_gen

    monkeypatch.setattr(get_cazy_pages, "get_paginiation_pages", mock_parse_pagination_page)
    monkeypatch.setattr(get_cazy_pages, "get_html_page", mock_get_pages)

    get_cazy_pages.parse_all_family(
        family=fam,
        cazy_home=cazy_home_url,
        args=args["args"],
    )


def test_pages_parse_all_format_error(args, protein_gen, cazy_home_url, monkeypatch):
    """Test parse_all_family when family url is not formatted correctly."""

    fam = Family(
        "testFam",
        "testClass",
        "testURL",
    )

    def mock_parse_pagination_page(*args, **kwargs):
        return [1, 2, 3, 4, 5, 6]

    def mock_get_pages(*args, **kwargs):
        return protein_gen

    monkeypatch.setattr(get_cazy_pages, "get_paginiation_pages", mock_parse_pagination_page)
    monkeypatch.setattr(get_cazy_pages, "get_html_page", mock_get_pages)

    get_cazy_pages.parse_all_family(
        family=fam,
        cazy_home=cazy_home_url,
        args=args["args"],
    )


def test_pages_parse_all_first_fam_parse(fam_url, args, protein_gen, cazy_home_url, monkeypatch):
    """Test parse_all_family when parsing the fam for the first time."""

    fam = Family(
        "testFam",
        "testClass",
        fam_url,
    )

    def mock_parse_pagination_page(*args, **kwargs):
        return [1, 2, 3, 4, 5, 6]

    def mock_get_pages(*args, **kwargs):
        return protein_gen

    monkeypatch.setattr(get_cazy_pages, "get_paginiation_pages", mock_parse_pagination_page)
    monkeypatch.setattr(get_cazy_pages, "get_html_page", mock_get_pages)

    get_cazy_pages.parse_all_family(
        family=fam,
        cazy_home=cazy_home_url,
        args=args["args"],
    )


# test get_pagination_pages()


def test_pages_get_pag_no_page(args, cazy_home_url, fam_url, fam_template, monkeypatch):
    """Test get_pagination_pages when no page is returned."""

    def mock_get_page(*args, **kwargs):
        return None, "error message"

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    res1, res2 = get_cazy_pages.get_pagination_pages(
        first_pagination_url=fam_url,
        family=fam_template,
        cazy_home=cazy_home_url,
        args=args["args"],
    )

    assert res1 == {'url': 'http://www.cazy.org/GH14_all.html\tCAZyClass\tFailed to connect to first pagination page for famName, therefore could not retrieve URLs to all pagination pages\terror message', 'format': None}
    assert res2 is None


def test_pages_get_page_successful(args, cazy_home_url, fam_url, fam_template, monkeypatch):
    """Test get_pagination_pages when successful"""

    def mock_get_page(*args, **kwargs):
        return "page", "error message"

    def mock_get_urls(*args, **kwargs):
        return ["url"]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(get_cazy_pages, "get_pagination_page_urls", mock_get_urls)

    res1 = get_cazy_pages.get_pagination_pages(
        first_pagination_url=fam_url,
        family=fam_template,
        cazy_home=cazy_home_url,
        args=args["args"],
    )

    assert res1 == ["url"]


# test get_pagination_page_urls()


def test_pages_get_urls_pag(fam_url, input_pages, cazy_home_url):
    """Test get_pagination_page_urls() when there is pagination."""

    test_input_path = input_pages / "pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    res = get_cazy_pages.get_pagination_page_urls(fam_url, soup, cazy_home_url, "testFam")
    assert len(res) == 40


def test_pages_get_urls_no_pag(fam_url, input_pages, cazy_home_url):
    """Test get_pagination_page_urls() when there is pagination."""

    test_input_path = input_pages / "no_pagination_pag.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    res = get_cazy_pages.get_pagination_page_urls(fam_url, soup, cazy_home_url, "testFam")
    assert res == ['http://www.cazy.org/GH14_all.html']


# test get_html_page()


def test_pages_get_html_pages_no_page(fam_url, args_pages, html_dir, monkeypatch):
    """Test get_html_pages() when no page is returned."""
    file_io.make_output_directory(html_dir, True, False)

    def mock_get_pages(*args, **kwargs):
        return None, "error"

    monkeypatch.setattr(crawler, "get_page", mock_get_pages)

    get_cazy_pages.get_html_page(fam_url, "testFam", args_pages["args"])

    file_io.make_output_directory(html_dir, True, False)


def test_pages_get_html_pages_deleted_fam(fam_url, args_pages, html_dir, input_pages, monkeypatch):
    """Test get_html_pages() when no page is returned."""
    file_io.make_output_directory(html_dir, True, False)

    test_input_path = input_pages / "deleted_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_pages(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(crawler, "get_page", mock_get_pages)

    get_cazy_pages.get_html_page(fam_url, "testFam", args_pages["args"])

    file_io.make_output_directory(html_dir, True, False)


def test_pages_get_html_pages_no_protein_count(fam_url, args_pages, html_dir, input_pages, monkeypatch):
    """Test get_html_pages() when no protien count can be retrieved."""
    file_io.make_output_directory(html_dir, True, False)

    test_input_path = input_pages / "no_protein_total.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_pages(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(crawler, "get_page", mock_get_pages)

    get_cazy_pages.get_html_page(fam_url, "testFam", args_pages["args"])

    file_io.make_output_directory(html_dir, True, False)


def test_pages_get_html_pages_empty(fam_url, args_pages, html_dir, input_pages, monkeypatch):
    """Test get_html_pages() when the family is empty."""
    file_io.make_output_directory(html_dir, True, False)

    test_input_path = input_pages / "empty_fam.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_pages(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(crawler, "get_page", mock_get_pages)

    get_cazy_pages.get_html_page(fam_url, "testFam", args_pages["args"])

    file_io.make_output_directory(html_dir, True, False)


def test_pages_get_html_pages_no_table(fam_url, args_pages, html_dir, input_pages, monkeypatch):
    """Test get_html_pages() when the CAZyme table is not populated."""
    file_io.make_output_directory(html_dir, True, False)

    test_input_path = str(input_pages).replace("test_get_pagination_urls", "test_parsing_proteins/no_cazyme_table.html")

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_pages(*args, **kwargs):
        return soup, None

    def mock_retry_connection(*args, **kwargs):
        return None

    monkeypatch.setattr(crawler, "get_page", mock_get_pages)
    monkeypatch.setattr(get_cazy_pages, "retry_failed_connections", mock_retry_connection)

    get_cazy_pages.get_html_page(fam_url, "testFam", args_pages["args"])

    file_io.make_output_directory(html_dir, True, False)


# test retry_failed_connections()


def test_pages_retry_failed_connections(fam_url, args_pages, html_dir, input_pages, monkeypatch):
    """Test retry_failed_connections()"""

    file_io.make_output_directory(html_dir, True, False)

    test_input_path = input_pages / "pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_pages(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(crawler, "get_page", mock_get_pages)

    get_cazy_pages.retry_failed_connections(
        fam_url,
        "TestFam",
        args_pages["args"],
        "filename.html"
    )

    file_io.make_output_directory(html_dir, True, False)


# test parse_family_by_kingdom()


def test_pages_kingdom_reparse_pagination(args, protein_gen, cazy_home_url, monkeypatch):
    """Test parse_family_by_kingdom() when reparsing a family and need to reparse the pagination."""

    test_fam = Family("famName", "CAZyClass", "http://www.cazy.org/GH14.html")
    test_fam.failed_pages = {"http://www.cazy.org/GH14_all.html": 1}

    def mock_get_pag(*args, **kwargs):
        return ["http://www.cazy.org/GH14_all.html"]

    def mock_get_pages(*args, **kwargs):
        return protein_gen

    monkeypatch.setattr(get_cazy_pages, "get_pagination_pages_kingdom", mock_get_pag)
    monkeypatch.setattr(get_cazy_pages, "get_html_page", mock_get_pages)

    get_cazy_pages.parse_family_by_kingdom(
        family=test_fam,
        cazy_home=cazy_home_url,
        args=args["args"],
        kingdoms=["Bacteria"],
    )


def test_pages_kingdom_first_parse_format_error(args, protein_gen, cazy_home_url, monkeypatch):
    """Test parse_family_by_kingdom() when parsing a fam for the first time and 
    URL is incorrectly formated."""

    test_fam = Family("famName", "CAZyClass", "http://www.caasdasdaszy.org/GH14.html")

    def mock_get_pag(*args, **kwargs):
        return ["http://www.cazy.org/GH14_all.html"]

    def mock_get_pages(*args, **kwargs):
        return protein_gen

    monkeypatch.setattr(get_cazy_pages, "get_pagination_pages_kingdom", mock_get_pag)
    monkeypatch.setattr(get_cazy_pages, "get_html_page", mock_get_pages)

    get_cazy_pages.parse_family_by_kingdom(
        family=test_fam,
        cazy_home=cazy_home_url,
        args=args["args"],
        kingdoms=["Bacteria"],
    )


def test_pages_kingdom_successful(args, protein_gen_success, cazy_home_url, monkeypatch):
    """Test parse_family_by_kingdom() when all is successful."""

    test_fam = Family("famName", "CAZyClass", "http://www.cazy.org/GH14.html")

    def mock_get_pag(*args, **kwargs):
        return ["http://www.cazy.org/GH14_all.html"]

    def mock_get_pages(*args, **kwargs):
        return protein_gen_success

    monkeypatch.setattr(get_cazy_pages, "get_pagination_pages_kingdom", mock_get_pag)
    monkeypatch.setattr(get_cazy_pages, "get_html_page", mock_get_pages)

    get_cazy_pages.parse_family_by_kingdom(
        family=test_fam,
        cazy_home=cazy_home_url,
        args=args["args"],
        kingdoms=["Bacteria"],
    )


# test get_pagination_pages_kingdom


def test_pages_king_pag_no_page(monkeypatch, fam_url, fam_template, cazy_home_url, args):
    """Test get_pagination_pages_kingdom when no page is returned."""

    def mock_get_page(*args, **kwargs):
        return None, "error"

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    get_cazy_pages.get_pagination_pages_kingdom(
        first_pagination_url=fam_url,
        family=fam_template,
        cazy_home=cazy_home_url,
        args=args["args"],
        kingdom="Bacteria"
    )


def test_pages_king_pag_successful(monkeypatch, fam_url, fam_template, cazy_home_url, args):
    """Test get_pagination_pages_kingdom when successful."""

    def mock_get_page(*args, **kwargs):
        return "page", "error"

    def mock_pag(*args, **kwargs):
        return [1, 2]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(get_cazy_pages, "get_tax_page_urls", mock_pag)

    get_cazy_pages.get_pagination_pages_kingdom(
        first_pagination_url=fam_url,
        family=fam_template,
        cazy_home=cazy_home_url,
        args=args["args"],
        kingdom="Bacteria"
    )


# test get_tax_pages_urls()


def test_pages_get_tax_pag_no_pag(pag_dir, fam_url, cazy_home_url):
    """Test get_tax_page_urls when there is no pagination."""

    test_input_path = pag_dir / "kb_no_pag.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    get_cazy_pages.get_tax_page_urls(fam_url, soup, cazy_home_url, "fam")


def test_pages_get_tax_pag_pag(pag_dir, fam_url, cazy_home_url):
    """Test get_tax_page_urls when there is pagination."""

    test_input_path = pag_dir / "kb_pag_page.html"

    with open(test_input_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    get_cazy_pages.get_tax_page_urls(fam_url, soup, cazy_home_url, "fam")
