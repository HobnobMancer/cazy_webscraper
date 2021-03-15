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
def args():
    args = {"args": Namespace(
        retries=2,
        timeout=45,
    )}
    return args


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
def cazy_class_page(input_dir):
    file_path = input_dir / "cazy_classpage.html"
    return file_path


@pytest.fixture
def cazy_class_page_no_fams(input_dir):
    file_path = input_dir / "cazy_classpage_no_fams.html"
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
def gh1_page_no_count(input_dir):
    file_path = input_dir / "cazy_gh1_page_no_count.html"
    return file_path


@pytest.fixture
def gh147_page(input_dir):
    file_path = input_dir / "cazy_gh147_page.html"
    return file_path


@pytest.fixture
def pag_page(input_dir):
    file_page = input_dir / "pagination_page.html"
    return file_page


@pytest.fixture
def pag_page_no_kingdom(input_dir):
    file_page = input_dir / "pagination_page_no_kingdom.html"
    return file_page


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
def protein_no_primary_gbks(input_dir):
    file_path = input_dir / "protein_no_primary_genbanks.html"
    return file_path


@pytest.fixture
def protein_with_no_gb(input_dir):
    file_path = input_dir / "protein_no_genbank.html"
    return file_path


@pytest.fixture
def protein_with_multi_primary(input_dir):
    file_path = input_dir / "protein_with_multi_primary.html"
    return file_path


@pytest.fixture
def protein_with_no_uniprot_no_pdb(input_dir):
    file_path = input_dir / "protein_no_uniprot_no_pdb.html"
    return file_path


@pytest.fixture
def test_fam():
    test_fam = crawler.Family("GH3", "Glycoside Hydrolases (GHs)", "http://www.cazy.org/GH1.html")
    return test_fam


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


def test_get_class_urls_fail(cazy_home_url, cazy_dictionary, monkeypatch, args):
    """Test get_cazy_class_urls home_page not returned"""

    def mock_get_home_page(*args, **kwargs):
        return [None, "error"]

    monkeypatch.setattr(crawler, "get_page", mock_get_home_page)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        crawler.get_cazy_classes(cazy_home_url, None, cazy_dictionary, args["args"])
    assert pytest_wrapped_e.type == SystemExit


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


# test get_cazy_family_urls


def test_get_family_urls_fail(args_subfam_true, monkeypatch, cazy_home_url):
    """Test get_cazy_family_urls when no page is returned."""

    def mock_get_page(*args, **kwargs):
        return [None, "error"]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    assert crawler.get_cazy_family_urls(
        "class_url",
        "class_name",
        cazy_home_url,
        args_subfam_true["args"],
    ) == (None, 'error', None)


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
        return [page, None]

    def mock_get_subfams(*args, **kwargs):
        return []

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_subfamily_links", mock_get_subfams)

    crawler.get_cazy_family_urls(
        "http://www.cazy.org/Glycoside-Hydrolases.html",
        "Glycoside Hydrolases (GHs)",
        "http://www.cazy.org/",
        args_subfam_true["args"],
    )


def test_get_family_urls_no_urls_no_subfams(
    cazy_class_page_no_fams,
    args_subfam_false,
    monkeypatch,
):
    """Test get_cazy_family_urls when no family URLs are retrieved and args.subfamilies is False."""
    with open(cazy_class_page_no_fams) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return [page, None]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    crawler.get_cazy_family_urls(
        "http://www.cazy.org/Glycoside-Hydrolases.html",
        "Glycoside Hydrolases (GHs)",
        "http://www.cazy.org/",
        args_subfam_false["args"],
    )


def test_get_family_urls_no_urls_subfams(
    cazy_class_page_no_fams,
    args_subfam_true,
    monkeypatch,
):
    """Test get_cazy_family_urls when no family URLs are retrieved and args.subfamilies is True."""
    with open(cazy_class_page_no_fams) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return [page, None]

    def mock_subfamilies(*args, **kwargs):
        return [
            "http://www.cazy.org/GH5_1.html",
            "http://www.cazy.org/GH5_2.html",
            "http://www.cazy.org/GH5_3.html",
        ]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_subfamily_links", mock_subfamilies)

    crawler.get_cazy_family_urls(
        "http://www.cazy.org/Glycoside-Hydrolases.html",
        "Glycoside Hydrolases (GHs)",
        "http://www.cazy.org/",
        args_subfam_true["args"],
    )


def test_get_family_urls_no_urls_subfams_none(
    cazy_class_page_no_fams,
    args_subfam_true,
    monkeypatch,
):
    """Test get_cazy_family_urls when no family or subfamily URLs are retrieved."""
    with open(cazy_class_page_no_fams) as fp:
        page = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return [page, None]

    def mock_subfamilies(*args, **kwargs):
        return None

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_subfamily_links", mock_subfamilies)

    crawler.get_cazy_family_urls(
        "http://www.cazy.org/Glycoside-Hydrolases.html",
        "Glycoside Hydrolases (GHs)",
        "http://www.cazy.org/",
        args_subfam_true["args"],
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


def test_parse_family_incorrect_url(args):
    """Tests parse_family() when the first page URL is in the incorrect format."""

    test_fam = crawler.Family(
        "GHs3ad",
        "Glycoside Hydrolases (GHs)",
        "http://www.caaazy.org/GHaaaaaa1.html",
    )

    crawler.parse_family(
        test_fam,
        "http://www.cazy.org/",
        None,
        ['bacteria'],
        args['args'],
        "sessions",
    )


def test_parse_family_no_page(monkeypatch, args):
    """Test parse_family() no page is returned."""

    def mock_get_page(*args, **kwargs):
        return None, "error_message"

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    test_fam = crawler.Family(
        "GH3",
        "Glycoside Hydrolases (GHs)",
        "http://www.cazssssy.org/GH3.html",
    )

    crawler.parse_family(
        test_fam,
        "http://www.cazy.org/",
        None,
        ['bacteria'],
        args['args'],
        "sessions",
    )


def test_parse_family_no_page_urls(monkeypatch, args):
    """Tests parse_family() when protein page URLs are retrieved"""

    def mock_get_page(*args, **kwargs):
        return "mock_page", None

    def mock_page_urls(*args, **kwargs):
        return [], 0

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_tax_page_urls", mock_page_urls)

    test_fam = crawler.Family("GH3", "Glycoside Hydrolases (GHs)", "http://www.cazy.org/GH1.html")

    crawler.parse_family(
        test_fam,
        "http://www.cazy.org",
        None,
        ['bacteria'],
        args['args'],
        "session",
    )


def test_parse_family_success(protein_gen, monkeypatch, args):
    """Test parse_family() when successful."""

    def mock_get_page(*args, **kwargs):
        return "mock_page", None

    def mock_page_urls(*args, **kwargs):
        return [], 0

    def mock_parse_proteins(*args, **kwargs):
        return [
            {"url": None, "error": None, "sql": None},
            {"url": None, "error": None, "sql": None},
        ]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_tax_page_urls", mock_page_urls)
    monkeypatch.setattr(crawler, "parse_proteins", mock_parse_proteins)

    test_family = crawler.Family("GH1", "Glycoside Hydrolases (GH)", "www.cazy.org/GH1.html")

    crawler.parse_family(
        test_family,
        "http://www.cazy.org/GH1.html",
        None,
        ['bacteria'],
        args['args'],
        "session",
    )


def test_parse_family_sql_url_errors(protein_gen, monkeypatch, args):
    """Test parse_family() when URL and SQL errors are raised when parsing a protein."""

    def mock_get_page(*args, **kwargs):
        return "mock_page", None

    def mock_page_urls(*args, **kwargs):
        return [], 0

    def mock_parse_proteins(*args, **kwargs):
        return [
            {"url": 'www.url_address', "error": 'error message', "sql": None},
            {"url": None, "error": 'error message', "sql": 'protein and sql'},
        ]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_tax_page_urls", mock_page_urls)
    monkeypatch.setattr(crawler, "parse_proteins", mock_parse_proteins)

    test_family = crawler.Family("GH1", "Glycoside Hydrolases (GH)", "www.cazy.org/GH1.html")

    crawler.parse_family(
        test_family,
        "http://www.cazy.org/GH1.html",
        None,
        ['bacteria'],
        args['args'],
        "session",
    )


def test_parse_family_previous_failed_pages(protein_gen, monkeypatch, args):
    """Test parse_family() when the family has previously failed families."""

    def mock_get_page(*args, **kwargs):
        return "mock_page", None

    def mock_page_urls(*args, **kwargs):
        return [], 0

    def mock_parse_proteins(*args, **kwargs):
        return [
            {"url": None, "error": None, "sql": None},
            {"url": None, "error": None, "sql": None},
        ]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_tax_page_urls", mock_page_urls)
    monkeypatch.setattr(crawler, "parse_proteins", mock_parse_proteins)

    test_family = crawler.Family(
        "GH1",
        "Glycoside Hydrolases (GH)",
        "www.cazy.org/GH1.html",
        {"pageURL1": 1, "pageURL2":2}
    )

    crawler.parse_family(
        test_family,
        "http://www.cazy.org/GH1.html",
        None,
        ['bacteria'],
        args['args'],
        "session",
    )


# test get_tax_page_urls()


def test_get_tax_page_urls_no_links(gh147_page):
    """Test get_tax_page_urls() on page with no links."""
    with open(gh147_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    crawler.get_tax_page_urls(
        "http://www.cazy.org/GH147_bacteria.html",
        soup,
        'bacteria',
        "http://www.cazy.org",
        'GH147',
    )


def test_get_tax_page_urls_no_pag(gh147_page):
    """Test get_tax_page_urls() on page with no pagination of proteins."""
    with open(gh147_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    crawler.get_tax_page_urls(
        "http://www.cazy.org/GH147_bacteria.html",
        soup,
        'bacteria',
        "http://www.cazy.org",
        'GH147',
    )


def test_get_tax_page_urls_page(pag_page):
    """Test get_tax_page_urls() on page with pagination of proteins."""
    with open(pag_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    crawler.get_tax_page_urls(
        "http://www.cazy.org/GH1_bacteria.html",
        soup,
        'bacteria',
        "http://www.cazy.org",
        'GH147',
    )


def test_get_tax_pages_no_kingdom(pag_page_no_kingdom):
    """Test get_tax_page_urls() when there are no proteins for the current Kingdom."""
    with open(pag_page_no_kingdom) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    crawler.get_tax_page_urls(
        "http://www.cazy.org/GH1_bacteria.html",
        soup,
        'bacteria',
        "http://www.cazy.org",
        'GH147',
    )


# test parse_proteins


def test_parse_proteins_none(monkeypatch, args):
    """Test parse_proteins when no page is returned."""

    def mock_get_page(*args, **kwargs):
        return None, "error message"

    def mock_row(*args, **kwargs):
        return {"url": None, "error": None, "sql": None}

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "row_to_protein", mock_row)

    crawler.parse_proteins("protein_page_url", "GH1", None, 'bacteria', args['args'], "session")


def test_parse_proteins(gh147_page, monkeypatch, args):
    """Tests parse_proteins when successful."""
    with open(gh147_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    def mock_row(*args, **kwargs):
        return {"url": None, "error": None, "sql": None}

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "row_to_protein", mock_row)

    crawler.parse_proteins("protein_url", "family", None, 'bacteria', args['args'], "session")


def test_parse_proteins_full(args, db_session, gh147_page, monkeypatch):
    """Test parsing proteins and adding to a database."""
    with open(gh147_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_page(*args, **kwargs):
        return soup, None

    monkeypatch.setattr(crawler, "get_page", mock_page)

    crawler.parse_proteins("protein_url", "family", None, 'bacteria', args['args'], db_session)


# test row_to_protein()


def test_row_to_protein_no_ecs(protein_without_ec, monkeypatch, args_subfam_true):
    """Test row_to_protein when no EC#s are listed."""
    with open(protein_without_ec) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", set(["Bacteroides caccae"]), 'bacteria', "session", args_subfam_true["args"])


def test_row_to_protein_ec(protein_with_ec, monkeypatch, args_subfam_true):
    """Test row_to_protein when EC#s are listed."""
    with open(protein_with_ec) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", set(["Bacteroides"]), 'bacteria', "session", args_subfam_true["args"])


def test_row_to_protein_gb_synoymns(protein_with_gb_synonyms, monkeypatch, args_subfam_true):
    """Test row_to_protein when protein has GenBank synonyms."""
    with open(protein_with_gb_synonyms) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", None, 'bacteria', "session", args_subfam_true["args"])


def test_row_to_protein_no_primary_gbk(protein_no_primary_gbks, monkeypatch, args_subfam_true):
    """Test row_to_protein when protein has GenBank synonyms."""
    with open(protein_no_primary_gbks,) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", None, 'bacteria', "session", args_subfam_true["args"])


def test_row_to_protein_with_multipl_primaries(protein_no_primary_gbks, monkeypatch, args_subfam_true):
    """Test row_to_protein when protein has multiple primary GenBank and UniProt accessions."""
    with open(protein_no_primary_gbks,) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", None, 'bacteria', "session", args_subfam_true["args"])



def test_row_to_protein_no_gb(protein_with_no_gb, monkeypatch, args_subfam_true):
    """Test row_to_protein when protein has no GenBank accession."""
    with open(protein_with_no_gb) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", None, 'bacteria', "session", args_subfam_true["args"])


def test_row_to_protein_no_gb_sql_error(protein_with_multi_primary, monkeypatch, args_subfam_true):
    """Test row_to_protein when protein has no GenBank accession and SQL raises an error."""
    with open(protein_with_multi_primary) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        raise TypeError

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", None, 'bacteria', "session", args_subfam_true["args"])


def test_row_to_protein_no_uniprot_no_pdb(protein_with_no_uniprot_no_pdb, monkeypatch, args_subfam_true):
    """Test row_to_protein when protein has no UniProt or PDB accessions."""
    with open(protein_with_no_uniprot_no_pdb) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", set(), 'bacteria', "session", args_subfam_true["args"])


def test_row_to_protein_no_uniprot_no_pdb_sql_error(protein_with_no_uniprot_no_pdb, monkeypatch, args_subfam_true):
    """Test row_to_protein when protein has no UniProt or PDB accessions and SQL raises an error."""
    with open(protein_with_no_uniprot_no_pdb) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        raise TypeError

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", None, 'bacteria', "session", args_subfam_true["args"])


def test_row_to_protein_gb_synoymns_raise_error(protein_with_gb_synonyms, monkeypatch, args_subfam_true):
    """Test row_to_protein when protein has GenBank synonyms, and SQL raises an error."""
    with open(protein_with_gb_synonyms) as fp:
        row = BeautifulSoup(fp, features="lxml")

    def mock_sql(*args, **kwargs):
        raise TypeError

    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_sql)

    crawler.row_to_protein(row, "GH147", None, 'bacteria', "session", args_subfam_true["args"])


# browser decorator and get_page


def test_browser_decorator():
    """Test browser_decorator to ensure proper handling if unsuccessful."""
    args = {"args": Namespace(timeout=10)}
    result = crawler.get_page('www.caz!!!!!!!!y.org', args["args"], max_tries=1)
    assert True == (result[0] is None) and (type(result[1]) is MissingSchema)


# unit tests for parsing HTML pages labelled under 'all'


def test_parse_fam_kingdom_all(monkeypatch, cazy_home_url, args_subfam_false, test_fam):
    """Test parse_family() when kingdoms = 'all'."""

    def mock_parse_all(*args, **kwargs):
        return test_fam, False, "failed scrapes", "sql_failures"

    monkeypatch.setattr(crawler, "parse_family_via_all_pages", mock_parse_all)

    crawler.parse_family(
        test_fam,
        cazy_home_url,
        None,
        'all',
        args_subfam_false['args'],
        'session',
    )


def test_parsing_all_incorrect_url(cazy_home_url, args_subfam_false):
    """Test parse_family_via_all_pages() when the URL format is wrong"""

    test_fam = crawler.Family("GH3", "Glycoside Hydrolases (GHs)", "http://www.cssazy.org/GH1.html")

    crawler.parse_family_via_all_pages(
        test_fam,
        cazy_home_url,
        None,
        args_subfam_false['args'],
        'session'
    )


def test_parsing_all_no_first_paginiation(test_fam, cazy_home_url, args_subfam_true, monkeypatch):
    """Test parse_family_via_all_pages() when no page retured for the first pagination page."""

    def mock_get_page(*args, **kwargs):
        return None, "error message"

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    crawler.parse_family_via_all_pages(
        test_fam,
        cazy_home_url,
        None,
        args_subfam_true['args'],
        'session'
    )


def test_parsing_all_no_paginiation_urls(
    test_fam,
    cazy_home_url,
    args_subfam_true,
    monkeypatch,
    gh1_page,
):
    """Test parse_family_via_all_pages() when no paginiation page URLs returned."""

    with open(gh1_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_pages(*args, **kwargs):
        return soup, None

    def mock_get_urls(*args, **kwargs):
        return [], 0

    monkeypatch.setattr(crawler, "get_page", mock_get_pages)
    monkeypatch.setattr(crawler, "get_paginiation_page_urls", mock_get_urls)

    crawler.parse_family_via_all_pages(
        test_fam,
        cazy_home_url,
        None,
        args_subfam_true['args'],
        'session'
    )


def test_successful_parse_all(test_fam, gh1_page, monkeypatch, args_subfam_true, cazy_home_url):
    """Test parse_family_via_all_pages() when parse is successful."""

    with open(gh1_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    def mock_page_urls(*args, **kwargs):
        return ["url1", "url2"], 0

    def mock_parse_proteins(*args, **kwargs):
        return [
            {"url": None, "error": None, "sql": None},
            {"url": None, "error": None, "sql": None},
        ]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_paginiation_page_urls", mock_page_urls)
    monkeypatch.setattr(crawler, "parse_proteins_from_all", mock_parse_proteins)

    crawler.parse_family(
        test_fam,
        cazy_home_url,
        None,
        'all',
        args_subfam_true['args'],
        "session",
    )


def test_successful_parse_retry(test_fam, gh1_page, monkeypatch, args_subfam_true, cazy_home_url):
    """Test parse_family_via_all_pages() when not the first scrape attempt."""

    test_fam.failed_pages = {"http://www.cazy.org/GH1.html": 1}

    with open(gh1_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    def mock_page_urls(*args, **kwargs):
        return ["url1", "url2"], 0

    def mock_parse_proteins(*args, **kwargs):
        return [
            {"url": None, "error": None, "sql": None},
            {"url": None, "error": None, "sql": None},
        ]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_paginiation_page_urls", mock_page_urls)
    monkeypatch.setattr(crawler, "parse_proteins_from_all", mock_parse_proteins)

    crawler.parse_family(
        test_fam,
        cazy_home_url,
        None,
        'all',
        args_subfam_true['args'],
        "session",
    )


def test_get_paginiation_urls(monkeypatch, cazy_home_url, gh1_page):
    """Test get_paginiation_page_urls() when urls successfully retrieved."""

    with open(gh1_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    crawler.get_paginiation_page_urls(
        "http://www.cazy.org/GH1_all.html",
        soup,
        cazy_home_url,
        "GH147"
    )


def test_get_paginiation_no_total(monkeypatch, cazy_home_url, gh1_page_no_count):
    """Test get_paginiation_page_urls() when number of proteins unsuccessfully retrieved."""

    with open(gh1_page_no_count) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    crawler.get_paginiation_page_urls(
        "http://www.cazy.org/GH1_all.html",
        soup,
        cazy_home_url,
        "GH147"
    )


def test_get_paginiation_no_last(monkeypatch, cazy_home_url, gh147_page):
    """Test get_paginiation_page_urls() when urls successfully retrieved."""

    with open(gh147_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    crawler.get_paginiation_page_urls(
        "http://www.cazy.org/GH147_all.html",
        soup,
        cazy_home_url,
        "GH147",
    )


def test_parse_proteins_all_no_page(monkeypatch, args_subfam_true):
    """Test parse_proteins_from_all() when no connection is made to CAZy."""

    def mock_get_page(*args, **kwargs):
        return None, "error message"

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    crawler.parse_proteins_from_all(
        "page_url",
        "GH147",
        None,
        "session",
        args_subfam_true["args"],
    )


def test_parse_proteins_from_all_successfully(monkeypatch, args_subfam_true, gh147_page):
    """Test parse_proteins_from_all() when successful."""

    with open(gh147_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup, None

    def mock_adding_to_db(*args, **kwargs):
        return {"url": None, "error": None, "sql": None}

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(sql_interface, "add_protein_to_db", mock_adding_to_db)

    crawler.parse_proteins_from_all(
        "http://www.cazy.org/GH147_all.html",
        "GH147",
        None,
        "session",
        args_subfam_true['args'],
    )
