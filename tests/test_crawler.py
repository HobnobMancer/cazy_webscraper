#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License

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
def input_dir(test_input_dir):
    dir_path = test_input_dir / "test_inputs_crawler"
    return dir_path


@pytest.fixture
def cazy_home_page(input_dir):
    file_path = input_dir / "cazy_homepage.html"
    return file_path


@pytest.fixture
def cazy_home_no_spip(input_dir):
    file_path = input_dir / "cazy_homepage_no_spip_out.html"
    return file_path


@pytest.fixture
def args_datasplit_family():
    argsdict = {
        "args": Namespace(
            data_split="family",
            subfamilies=True,
        )
    }
    return argsdict


@pytest.fixture
def cazy_class_page(input_dir):
    file_path = input_dir / "cazy_classpage.html"
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
def protein_gen():
    result_list = [
        {"url": None, "error": None, "sql": None},
        {"url": "www.cazy.org/GH1.html", "error": "no internet connection", "sql": None},
        {"url": None, "error": "sql error", "sql": "protein_name"},
    ]
    return (_ for _ in result_list)


@pytest.fixture
def gh1_page(input_dir):
    file_path = input_dir / "cazy_gh1_page.html"
    return file_path


@pytest.fixture
def gh147_page(input_dir):
    file_path = input_dir / "cazy_gh147_page.html"
    return file_path


@pytest.fixture
def no_links_page(input_dir):
    file_path = input_dir / "cazy_no_links.html"
    return file_path


@pytest.fixture
def protein_without_ec(input_dir):
    file_path = input_dir / "protein_without_ec.html"
    return file_path


@pytest.fixture
def protein_with_ec(input_dir):
    file_path = input_dir / "protein_with_ec.html"
    return file_path


@pytest.fixture
def protein_with_gb_synonyms(input_dir):
    file_path = input_dir / "protein_with_genbank_synonyms.html"
    return file_path


@pytest.fixture
def protein_with_no_gb(input_dir):
    file_path = input_dir / "protein_no_genbank.html"
    return file_path


@pytest.fixture
def protein_with_no_uniprot_no_pdb(input_dir):
    file_path = input_dir / "protein_no_uniprot_no_pdb.html"
    return file_path


# test classes


def test_cazy_class():
    """Test building a CAZy class instance."""
    new_class = crawler.CazyClass("GH", "class_url", 0)
    new_class_failed_families = crawler.CazyClass("GH", "class_url", 0, {"GH1": "GH1_url"})

    assert new_class.name == "GH"
    assert new_class_failed_families.name == "GH"


def test_family_get_name():
    """Tests get family name for Family."""
    family = crawler.Family("GH1", "Glycoside_Hydrolases(GH)", "class_url")

    exepected_name = "GH1"

    assert exepected_name == family.name


# test get_cazy_class_urls


def test_get_class_urls_exclusions_none(
    cazy_home_url,
    cazy_home_page,
    cazy_dictionary,
    monkeypatch,
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
        1,
        cazy_dictionary,
    )
    assert len(result) == 6


def test_get_class_urls_fail(cazy_home_url, cazy_dictionary, monkeypatch):
    """Test get_cazy_class_urls home_page not returned"""

    def mock_get_home_page(*args, **kwargs):
        return [None, "error"]

    monkeypatch.setattr(crawler, "get_page", mock_get_home_page)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        crawler.get_cazy_classes(cazy_home_url, None, 1, cazy_dictionary)
    assert pytest_wrapped_e.type == SystemExit


def test_get_class_urls_exclusions_given(
    cazy_home_url,
    cazy_home_page,
    monkeypatch,
    cazy_dictionary,
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
        1,
        cazy_dictionary,
    )

    assert len(result) == 5


def test_get_class_urls_attribute(
    cazy_home_url,
    cazy_home_no_spip,
    monkeypatch,
    cazy_dictionary,
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
            1,
            cazy_dictionary,
        )
    assert pytest_wrapped_e.type == SystemExit


# test get_cazy_family_urls


def test_get_family_urls_fail(args_datasplit_family, monkeypatch):
    """Test get_cazy_family_urls when no page is returned."""

    def mock_get_page(*args, **kwargs):
        return [None, "error"]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    assert crawler.get_cazy_family_urls(
        "class_url",
        "cazy_home_url",
        "class_name",
        args_datasplit_family["args"],
    ) == (None, 'error', None)


def test_get_family_urls_success(
    cazy_class_page,
    args_datasplit_family,
    family_urls,
    monkeypatch,
):
    """Test get_cazy_family_urls when successful, and subfamilies is True."""
    with open(cazy_class_page) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return [page, None]

    def mock_get_subfams(*args, **kwargs):
        return []

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_subfamily_links", mock_get_subfams)

    crawler.get_cazy_family_urls(
        "http://www.cazy.org/Glycoside-Hydrolases.html",
        "Glycoside Hydrolases (GHs)",
        "http://www.cazy.org/",
        args_datasplit_family["args"],
    )


# test get_subfamily_links


def test_get_subfam_links_len_0(family_h3_element, subfamily_urls, ):
    """Test get_subfamily_links when no links are retrieved."""

    assert subfamily_urls == crawler.get_subfamily_links(
        family_h3_element,
        "http://www.cazy.org",
    )


def test_get_subfam_links_urls(no_subfam_h3_element, ):
    """Test get_subfamily_links when urls are retrieved."""

    assert None is crawler.get_subfamily_links(
        no_subfam_h3_element,
        "http://www.cazy.org",
    )


# test parse_family()


def test_parse_family_incorrect_url():
    """Tests parse_family() when the first page URL is in the incorrect format."""

    test_fam = crawler.Family(
        "GHs3ad",
        "Glycoside Hydrolases (GHs)",
        "http://www.caaazy.org/GHaaaaaa1.html",
    )

    crawler.parse_family(
        test_fam,
        "http://www.cazy.org/",
        1,
        "sessions",
    )


def test_parse_family_no_page(monkeypatch):
    """Test parse_family() no page is returned."""

    def mock_get_page(*args, **kwargs):
        return

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    test_fam = crawler.Family(
        "GH3",
        "Glycoside Hydrolases (GHs)",
        "http://www.cazssssy.org/GH3.html",
    )

    crawler.parse_family(
        test_fam,
        "http://www.cazy.org/",
        1,
        "sessions",
    )


def test_parse_family_no_page_urls(monkeypatch):
    """Tests parse_family() when protein page URLs are retrieved"""

    def mock_get_page(*args, **kwargs):
        return "mock_page", None

    def mock_page_urls(*args, **kwargs):
        return [], 0

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_protein_page_urls", mock_page_urls)

    test_fam = crawler.Family("GH3", "Glycoside Hydrolases (GHs)", "http://www.cazy.org/GH1.html")

    crawler.parse_family(
        test_fam,
        "http://www.cazy.org",
        2,
        "session",
    )


def test_parse_family_success(protein_gen, monkeypatch):
    """Test parse_family() when successful."""

    def mock_get_page(*args, **kwargs):
        return "mock_page", None

    def mock_page_urls(*args, **kwargs):
        return [], 0

    def mock_parse_proteins(*args, **kwargs):
        return protein_gen

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_protein_page_urls", mock_page_urls)
    monkeypatch.setattr(crawler, "parse_proteins", mock_parse_proteins)

    test_family = crawler.Family("GH1", "Glycoside Hydrolases (GH)", "www.cazy.org/GH1.html")

    crawler.parse_family(
        test_family,
        "http://www.cazy.org/GH1.html",
        2,
        "session",
    )


# test get_protein_page_urls()


def test_get_protein_page_urls_no_links(gh147_page):
    """Test get_protein_page_urls() on page with no links."""
    with open(gh147_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    crawler.get_protein_page_urls("http://www.cazy.org/GH147_all.html", soup, "http://www.cazy.org")


def test_get_protein_page_urls_no_pag(gh147_page):
    """Test get_protein_page_urls() on page with no pagination of proteins."""
    with open(gh147_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    crawler.get_protein_page_urls("http://www.cazy.org/GH147_all.html", soup, "http://www.cazy.org")


def test_get_protein_page_urls_pag(gh1_page):
    """Test get_protein_page_urls() on page with pagination of proteins."""
    with open(gh1_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    crawler.get_protein_page_urls("http://www.cazy.org/GH1_all.html", soup, "http://www.cazy.org")


# test parse_proteins


def test_parse_proteins_none(monkeypatch):
    """Test parse_proteins when no page is returned."""

    def mock_get_page(*args, **kwargs):
        return [None, "error"]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    assert True is isinstance(
        crawler.parse_proteins("protein_url", "family", "session"),
        types.GeneratorType,
    )


def test_parse_proteins(gh147_page, monkeypatch):
    """Tests parse_proteins when successful."""
    with open(gh147_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    crawler.parse_proteins("protein_url", "family", "session")


# test row_to_protein()


def test_row_to_protein_no_ecs(protein_without_ec, monkeypatch):
    """Test row_to_protein when no EC#s are listed."""
    with open(protein_without_ec) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", "session")


def test_row_to_protein_ec(protein_with_ec, monkeypatch):
    """Test row_to_protein when EC#s are listed."""
    with open(protein_with_ec) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", "session")


def test_row_to_protein_gb_synoymns(protein_with_gb_synonyms, monkeypatch):
    """Test row_to_protein when protein has GenBank synonyms."""
    with open(protein_with_gb_synonyms) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", "session")


def test_row_to_protein_no_gb(protein_with_no_gb, monkeypatch):
    """Test row_to_protein when protein has no GenBank accession."""
    with open(protein_with_no_gb) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", "session")


def test_row_to_protein_no_uniprot_no_pdb(protein_with_no_uniprot_no_pdb, monkeypatch):
    """Test row_to_protein when protein has no UniProt or PDB accessions."""
    with open(protein_with_no_uniprot_no_pdb) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", "session")


# browser decorator and get_page


def test_browser_decorator():
    """Test browser_decorator to ensure proper handling if unsuccessful."""

    result = crawler.get_page('www.caz!!!!!!!!y.org')
    assert True == (result[0] is None) and (type(result[1]) is MissingSchema)
