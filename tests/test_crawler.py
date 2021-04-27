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
"""Tests the script scraper.crawler submodule which crawls through the CAZy database.

These test are intened to be run from the root of the repository using:
pytest -v
"""

import types
import pytest

from argparse import Namespace
from requests.exceptions import MissingSchema

from bs4 import BeautifulSoup

from scraper import crawler
from scraper.sql import sql_interface


@pytest.fixture
def args():
    args = {"args": Namespace(
        retries=2,
        timeout=45,
    )}
    return args


@pytest.fixture
def input_dir(test_input_dir):
    dir_path = test_input_dir /"test_inputs_crawler"
    return dir_path


@pytest.fixture
def cazy_home_page(input_dir):
    file_path = input_dir / "class_url_pages" / "cazy_homepage.html"
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


# HTML pages for testing retrieval of URLs to CAZy family pages from a CAZy class summary page


@pytest.fixture
def cazy_class_page(input_dir):
    file_path = input_dir / "family_url_pages" / "cazy_classpage.html"
    return file_path


@pytest.fixture
def cazy_class_page_no_urls(input_dir):
    file_path = input_dir / "family_url_pages" / "cazy_classpage_no_urls.html"
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
def subfamily_urls(input_dir):
    file_path = input_dir / "subfamily_urls.txt"
    with open(file_path, "r") as fh:
        fam_string = fh.read()
    fam_string = fam_string[1:-1]
    fam_string = fam_string.replace("'", "")
    fam_list = fam_string.split(", ")
    return fam_list


@pytest.fixture
def no_subfam_h3_element(input_dir):
    file_path = input_dir / "cazy_classpage_no_subfams.html"
    with open(file_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    return [_ for _ in
            soup.find_all("h3", {"class": "spip"}) if
            str(_.contents[0]) == "Tables for Direct Access"][0]


@pytest.fixture
def gh147_page(input_dir):
    file_path = input_dir / "cazy_gh147_page.html"
    return file_path


@pytest.fixture
def protein_without_ec(input_dir):
    file_path = input_dir / "row_to_protein" / "protein_without_ec.html"
    return file_path


@pytest.fixture
def protein_with_ec(input_dir):
    file_path = input_dir / "row_to_protein" / "protein_with_ec.html"
    return file_path


@pytest.fixture
def protein_no_primary_gbks(input_dir):
    file_path = input_dir / "row_to_protein" / "protein_no_primary_genbanks.html"
    return file_path


@pytest.fixture
def protein_no_gbks(input_dir):
    file_path = input_dir / "row_to_protein" / "protein_no_gbks.html"
    return file_path


@pytest.fixture
def protein_multiple_primary_gbks(input_dir):
    file_path = input_dir / "row_to_protein" / "protein_multiple_primary_gbks.html"
    return file_path


# test classes


def test_cazy_class():
    """Test building a CAZy class instance."""
    new_class = crawler.CazyClass("GH", "class_url", 0)
    new_class_failed_families = crawler.CazyClass("GH", "class_url", 0, {"GH1": "GH1_url"})

    assert new_class.name == "GH"
    assert new_class_failed_families.name == "GH"
    new_class
    repr(new_class)


def test_family_get_name():
    """Tests get family name for Family."""
    fam_1 = crawler.Family("GH1", "Glycoside_Hydrolases(GH)", "class_url")
    fam_2 = crawler.Family("GH2", "Glycoside_Hydrolases(GH)", "class_url", {"url":0})

    assert "GH1" == fam_1.name
    assert "GH2" == fam_2.name
    fam_1
    repr(fam_1)


# test get_cazy_class_urls


def test_get_class_urls_fail(cazy_home_url, cazy_dictionary, monkeypatch, args):
    """Test get_cazy_class_urls home_page not returned"""

    def mock_get_home_page(*args, **kwargs):
        return None, "error"

    monkeypatch.setattr(crawler, "get_page", mock_get_home_page)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        crawler.get_cazy_classes(cazy_home_url, None, cazy_dictionary, args["args"])
    assert pytest_wrapped_e.type == SystemExit


def test_get_class_urls_exclusions_none(
    cazy_home_url,
    cazy_home_page,
    cazy_dictionary,
    monkeypatch,
    args,
):
    """Test get_cazy_class_urls when excluded_classess is None."""
    with open(cazy_home_page, "r") as fp:
        home_page = BeautifulSoup(fp, features="lxml")

        def mock_get_home_page(*args, **kwargs):
            return [home_page, None]

    monkeypatch.setattr(crawler, "get_page", mock_get_home_page)

    result = crawler.get_cazy_classes(
        cazy_home_url,
        None,
        cazy_dictionary,
        args["args"],
    )
    assert len(result) == 6


def test_get_class_urls_exclusions_given(
    cazy_home_url,
    cazy_home_page,
    monkeypatch,
    cazy_dictionary,
    args,
):
    """Test get_cazy_class_urls when excluded_classess is not None."""
    with open(cazy_home_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    exclusions = ["<strong>Glycoside Hydrolases (GHs)</strong>"]

    def mock_get_home_page(*args, **kwargs):
        return [soup, None]

    monkeypatch.setattr(crawler, "get_page", mock_get_home_page)

    result = crawler.get_cazy_classes(
        cazy_home_url,
        exclusions,
        cazy_dictionary,
        args["args"],
    )

    assert len(result) == 5


def test_get_class_urls_attribute(
    cazy_home_url,
    cazy_home_no_spip,
    monkeypatch,
    cazy_dictionary,
    args,
):
    """Test get_cazy_class_urls when attribute error is raised."""
    with open(cazy_home_no_spip) as fp:
        soup = fp.read()

    exclusions = ["<strong>Glycoside Hydrolases (GHs)</strong>"]

    def mock_get_home_page(*args, **kwargs):
        return [soup, None]

    monkeypatch.setattr(crawler, "get_page", mock_get_home_page)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        crawler.get_cazy_classes(
            cazy_home_url,
            exclusions,
            cazy_dictionary,
            args["args"],
        )
    assert pytest_wrapped_e.type == SystemExit


def test_get_class_urls_no_urls(
    cazy_home_url,
    cazy_home_no_urls,
    monkeypatch,
    cazy_dictionary,
    args,
):
    """Test get_cazy_class_urls when no class urls are returned from the HTML webpage."""
    with open(cazy_home_no_urls) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    exclusions = ["<strong>Glycoside Hydrolases (GHs)</strong>"]

    def mock_get_home_page(*args, **kwargs):
        return [soup, None]

    monkeypatch.setattr(crawler, "get_page", mock_get_home_page)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        crawler.get_cazy_classes(
            cazy_home_url,
            exclusions,
            cazy_dictionary,
            args["args"],
        )
    assert pytest_wrapped_e.type == SystemExit


# test get_cazy_family_urls


def test_get_family_urls_no_page(args_subfam_true, monkeypatch, cazy_home_url):
    """Test get_cazy_family_urls when no page is returned."""

    def mock_get_page(*args, **kwargs):
        return None, "error"

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    assert crawler.get_cazy_family_urls(
        "class_url",
        "class_name",
        cazy_home_url,
        args_subfam_true["args"],
    ) == (None, 'error', None)


def test_get_family_urls_no_urls(args_subfam_false, monkeypatch, cazy_class_page_no_urls):
    """Tests get_cazy_family_urls when no Family URls are returned."""
    with open(cazy_class_page_no_urls) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return page, None

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    fam, message, incorrect_urls = crawler.get_cazy_family_urls(
        "http://www.cazy.org/Glycoside-Hydrolases.html",
        "Glycoside Hydrolases (GHs)",
        "http://www.cazy.org/",
        args_subfam_false["args"],
    )
    assert fam is None
    assert message == "Failed to retrieve URLs to CAZy families for Glycoside Hydrolases (GHs)\n"
    assert incorrect_urls is None


def test_get_family_urls_no_urls_no_subfam_true(
    args_subfam_true,
    monkeypatch,
    cazy_class_page_no_urls,
):
    """Tests get_cazy_family_urls when no Family URls are returned, and retrieving subfamiles
    but none are retrieved."""
    with open(cazy_class_page_no_urls) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return page, None

    def mock_subfams(*args, **kwargs):
        return

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_subfamily_links", mock_subfams)

    fam, message, incorrect_urls = crawler.get_cazy_family_urls(
        "http://www.cazy.org/Glycoside-Hydrolases.html",
        "Glycoside Hydrolases (GHs)",
        "http://www.cazy.org/",
        args_subfam_true["args"],
    )

    assert fam is None
    assert message == (
        "Could not retrieve family and subfamily URLs from class page of "
        "Glycoside Hydrolases (GHs)\n"
    )
    assert incorrect_urls is None


def test_get_family_urls_no_urls_subfam_true(
    args_subfam_true,
    monkeypatch,
    cazy_class_page_no_urls,
):
    """Tests get_cazy_family_urls when no Family URls are returned, and retrieving subfamiles
    and are retrieved."""
    with open(cazy_class_page_no_urls) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return page, None

    def mock_subfams(*args, **kwargs):
        return ["http://www.cazy.org/GH5_1.html"]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_subfamily_links", mock_subfams)

    fam, message, incorrect_urls = crawler.get_cazy_family_urls(
        "http://www.cazy.org/Glycoside-Hydrolases.html",
        "Glycoside Hydrolases (GHs)",
        "http://www.cazy.org/",
        args_subfam_true["args"],
    )


def test_get_family_urls_success(
    cazy_class_page,
    args_subfam_true,
    family_urls,
    monkeypatch,
):
    """Test get_cazy_family_urls when successful, and subfamilies is True."""
    with open(cazy_class_page) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return page, None

    def mock_get_subfams(*args, **kwargs):
        return []

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_subfamily_links", mock_get_subfams)

    fam, message, incorrect_urls = crawler.get_cazy_family_urls(
        "http://www.cazy.org/Glycoside-Hydrolases.html",
        "Glycoside Hydrolases (GHs)",
        "http://www.cazy.org/",
        args_subfam_true["args"],
    )

    assert incorrect_urls == [
        (
            "http://www.cazy.org//GH1testtesttesttesttest.html\tGlycoside Hydrolases (GHs)\t"
            "Format of the URL is incorrect\t'NoneType' object has no attribute 'group'"
        ),
        (
            "http://www.cazy.org//GH2aaa.html\tGlycoside Hydrolases (GHs)\t"
            "Format of the URL is incorrect\t'NoneType' object has no attribute 'group'"
        ),
        (
            "http://www.cazy.org//GHtesttesttest3.html\tGlycoside Hydrolases (GHs)\t"
            "Format of the URL is incorrect\t'NoneType' object has no attribute 'group'"
        ),
    ]


# test get_subfamily_links


def test_get_subfam_links_successul(family_h3_element, subfamily_urls):
    """Test get_subfamily_links when links are retrieved."""

    res = crawler.get_subfamily_links(
        family_h3_element,
        "http://www.cazy.org",
    )
    print(subfamily_urls)
    assert res == subfamily_urls


def test_get_subfam_links_no_urls(no_subfam_h3_element, ):
    """Test get_subfamily_links when no urls are retrieved."""

    assert None is crawler.get_subfamily_links(
        no_subfam_h3_element,
        "http://www.cazy.org",
    )


# test row_to_protein()


def test_row_to_protein_no_ecs(protein_without_ec, monkeypatch, args_subfam_true):
    """Test row_to_protein when no EC#s are listed."""
    with open(protein_without_ec) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        kingdom='bacteria',
        ec_filters=["EC1.2.3.4"],
        session="session",
        args=args_subfam_true["args"],
    )


def test_row_to_protein_ec(protein_with_ec, monkeypatch, args_subfam_true):
    """Test row_to_protein when EC#s are listed."""
    with open(protein_with_ec) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        kingdom='bacteria',
        ec_filters=["1.2.3.4"],
        session="session",
        args=args_subfam_true["args"],
    )


def test_row_to_protein_no_genbanks(protein_no_gbks, monkeypatch, args_subfam_true):
    """Test row_to_protein when a protein has no primary or non-primary GenBank accessions."""
    with open(protein_no_gbks) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        kingdom='bacteria',
        ec_filters=["3.2.1.-"],
        session="session",
        args=args_subfam_true["args"],
    )


def test_row_to_protein_no_primary_gbk_with_nonprimary(
    protein_no_primary_gbks,
    monkeypatch,
    args_subfam_true,
):
    """Test row_to_protein when a protein has no primary GenBank accession but has non-primary
    GenBank accessions."""
    with open(protein_no_primary_gbks,) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        kingdom='bacteria',
        ec_filters=["3.2.1.-"],
        session="session",
        args=args_subfam_true["args"],
    )


def test_row_to_protein_with_multipl_primaries(
    protein_multiple_primary_gbks,
    monkeypatch,
    args_subfam_true,
):
    """Test row_to_protein when protein has multiple primary GenBank and UniProt accessions."""
    with open(protein_multiple_primary_gbks) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        kingdom='bacteria',
        ec_filters=["3.2.1.-"],
        session="session",
        args=args_subfam_true["args"],
    )


def test_row_to_protein_sql_error(
    protein_multiple_primary_gbks,
    monkeypatch,
    args_subfam_true,
):
    """Test row_to_protein when SQL raises an error."""
    with open(protein_multiple_primary_gbks) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        raise TypeError

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        kingdom='bacteria',
        ec_filters=["3.2.1.-"],
        session="session",
        args=args_subfam_true["args"],
    )


# test row_to_protein dict


def test_row_to_protein_in_dict_no_ecs(protein_without_ec, monkeypatch, args_subfam_true):
    """Test row_to_protein when no EC#s are listed."""
    with open(protein_without_ec) as fp:
        row = BeautifulSoup(fp, features="lxml")

    crawler.row_to_protein_in_dict(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        ec_filters=["EC1.2.3.4"],
        session={},
    )


def test_row_to_protein_in_dict_ec(protein_with_ec, monkeypatch, args_subfam_true):
    """Test row_to_protein when EC#s are listed."""
    with open(protein_with_ec) as fp:
        row = BeautifulSoup(fp, features="lxml")

    crawler.row_to_protein_in_dict(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        ec_filters=["1.2.3.4"],
        session={},
    )


def test_row_to_protein_in_dict_no_genbanks(protein_no_gbks, monkeypatch, args_subfam_true):
    """Test row_to_protein when a protein has no primary or non-primary GenBank accessions."""
    with open(protein_no_gbks) as fp:
        row = BeautifulSoup(fp, features="lxml")

    crawler.row_to_protein_in_dict(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        ec_filters=["3.2.1.-"],
        session={},
    )


def test_row_to_protein_in_dict_no_primary_gbk_with_nonprimary(
    protein_no_primary_gbks,
    monkeypatch,
    args_subfam_true,
):
    """Test row_to_protein when a protein has no primary GenBank accession but has non-primary
    GenBank accessions."""
    with open(protein_no_primary_gbks,) as fp:
        row = BeautifulSoup(fp, features="lxml")

    crawler.row_to_protein_in_dict(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        ec_filters=["3.2.1.-"],
        session={},
    )


def test_row_to_protein_in_dict_with_multipl_primaries(
    protein_multiple_primary_gbks,
    monkeypatch,
    args_subfam_true,
):
    """Test row_to_protein when protein has multiple primary GenBank and UniProt accessions."""
    with open(protein_multiple_primary_gbks) as fp:
        row = BeautifulSoup(fp, features="lxml")

    crawler.row_to_protein_in_dict(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        ec_filters=["3.2.1.-"],
        session={},
    )


def test_row_to_protein_in_dict_sql_error(
    protein_multiple_primary_gbks,
    monkeypatch,
    args_subfam_true,
):
    """Test row_to_protein when SQL raises an error."""
    with open(protein_multiple_primary_gbks) as fp:
        row = BeautifulSoup(fp, features="lxml")

    crawler.row_to_protein_in_dict(
        row=row,
        family_name="GH147",
        taxonomy_filters=set(["Bacteroides caccae"]),
        ec_filters=["3.2.1.-"],
        session={},
    )


# test get_all_accessions


def test_get_all_accessions(gh147_page):
    with open(gh147_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    cazyme_table = soup.select("table")[1]
    for row in cazyme_table.select("tr"):
        try:
            if (row.attrs["class"] == ['royaume']) and (row.text.strip() != 'Top'):
                continue
        except KeyError:
            pass

        if ('class' not in row.attrs) and ('id' not in row.attrs):  # row contains protein data
            gbk_bs_elements = [
                _ for _ in row.select("td")[3].contents if getattr(_, "name", None) != "br"
            ]
            uni_bs_elements = [
                _ for _ in row.select("td")[4].contents if getattr(_, "name", None) != "br"
            ]
            pdb_bs_elements = [
                _ for _ in row.select("td")[5].contents if getattr(_, "name", None) != "br"
            ]

            gbk_nonprimary = crawler.get_all_accessions(gbk_bs_elements)
            uni_nonprimary = crawler.get_all_accessions(uni_bs_elements)
            pdb_accessions = crawler.get_all_accessions(pdb_bs_elements)


# browser decorator and get_page


def test_browser_decorator():
    """Test browser_decorator to ensure proper handling if unsuccessful."""
    args = {"args": Namespace(timeout=10)}
    result = crawler.get_page('www.caz!!!!!!!!y.org', args["args"], max_tries=1)
    assert True == (result[0] is None) and (type(result[1]) is MissingSchema)
