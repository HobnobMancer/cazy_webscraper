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
def protein_list():
    return [
        crawler.Protein("protein_name", "GH1", "1.2.3.4", "organism"),
        crawler.Protein(
            "protein",
            "GH1",
            "",
            "organism",
            {"GenBank": ["link1"], "UniProt": ["link2"], "PDB": ["link3"]},
        ),
    ]


@pytest.fixture
def protein_gen(protein_list):
    return (_ for _ in protein_list)


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


# test the classes Protein and Family


def test_protein_get_protein_dict():
    """Test the Protein class __str__ and __repr__."""
    protein = crawler.Protein(
        "protein_name",
        "GH1",
        ["1.2.3.4"],
        "organism",
        {"GenBank": ["link1", "linka", "linkb"], "UniProt": ["link2"], "PDB/3D": ["link3"]},
    )

    protein_dict = protein.get_protein_dict()

    expected_protein_dict = {
        'Protein_name': ['protein_name'],
        'CAZy_family': ['GH1'],
        'EC#': ['1.2.3.4'],
        'Source_organism': ['organism'],
        'GenBank': ['link1,\nlinka,\nlinkb'],
        'UniProt': ['link2'],
        'PDB/3D': ['link3'],
    }

    assert expected_protein_dict == protein_dict


def test_family_get_name():
    """Tests get family name for Family."""
    family = crawler.Family("GH1", "Glycoside_Hydrolases(GH)")

    family_name = family.get_family_name()

    exepected_name = "GH1"

    assert exepected_name == family_name


# test get_cazy_class_urls


def test_get_class_urls_exclusions_none(cazy_home_url, cazy_home_page, null_logger, monkeypatch):
    """Test get_cazy_class_urls when excluded_classess is None."""
    with open(cazy_home_page, "r") as fp:
        home_page = BeautifulSoup(fp, features="lxml")

        def mock_get_home_page(*args, **kwargs):
            return [home_page, None]

    monkeypatch.setattr(crawler, "get_page", mock_get_home_page)

    expected_result = [
        'http://www.cazy.org/Glycoside-Hydrolases.html',
        'http://www.cazy.org/GlycosylTransferases.html',
        'http://www.cazy.org/Polysaccharide-Lyases.html',
        'http://www.cazy.org/Carbohydrate-Esterases.html',
        'http://www.cazy.org/Auxiliary-Activities.html',
        'http://www.cazy.org/Carbohydrate-Binding-Modules.html',
    ]

    assert expected_result == crawler.get_cazy_class_urls(
        cazy_home_url,
        None,
        1,
        null_logger,
    )


def test_get_class_urls_fail(cazy_home_url, null_logger, monkeypatch):
    """Test get_cazy_class_urls home_page not returned"""

    def mock_get_home_page(*args, **kwargs):
        return [None, "error"]

    monkeypatch.setattr(crawler, "get_page", mock_get_home_page)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        crawler.get_cazy_class_urls(cazy_home_url, None, 1, null_logger)
    assert pytest_wrapped_e.type == SystemExit


def test_get_class_urls_exclusions_given(
    cazy_home_url,
    cazy_home_page,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_class_urls when excluded_classess is not None."""
    with open(cazy_home_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    exclusions = ["<strong>Glycoside Hydrolases (GHs)</strong>"]

    def mock_get_home_page(*args, **kwargs):
        return [soup, None]

    monkeypatch.setattr(crawler, "get_page", mock_get_home_page)

    expected_result = [
        'http://www.cazy.org/GlycosylTransferases.html',
        'http://www.cazy.org/Polysaccharide-Lyases.html',
        'http://www.cazy.org/Carbohydrate-Esterases.html',
        'http://www.cazy.org/Auxiliary-Activities.html',
        'http://www.cazy.org/Carbohydrate-Binding-Modules.html',
    ]

    assert expected_result == crawler.get_cazy_class_urls(
        cazy_home_url,
        exclusions,
        1,
        null_logger,
    )


def test_get_class_urls_attribute(
    cazy_home_url,
    cazy_home_no_spip,
    null_logger,
    monkeypatch,
):
    """Test get_cazy_class_urls when attribute error is raised."""
    with open(cazy_home_no_spip) as fp:
        soup = fp.read()

    exclusions = ["<strong>Glycoside Hydrolases (GHs)</strong>"]

    def mock_get_home_page(*args, **kwargs):
        return [soup, None]

    monkeypatch.setattr(crawler, "get_page", mock_get_home_page)

    assert None is crawler.get_cazy_class_urls(
        cazy_home_url,
        exclusions,
        1,
        null_logger,
    )


# test get_cazy_family_urls


def test_get_family_urls_fail(args_datasplit_family, null_logger, monkeypatch):
    """Test get_cazy_family_urls when no page is returned."""

    def mock_get_page(*args, **kwargs):
        return [None, "error"]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    assert None is crawler.get_cazy_family_urls(
        "class_url",
        "cazy_home_url",
        "class_name",
        args_datasplit_family["args"],
        null_logger,
    )


def test_get_family_urls_success(
    cazy_class_page,
    args_datasplit_family,
    family_urls,
    null_logger,
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

    assert family_urls == crawler.get_cazy_family_urls(
        "class_url",
        "http://www.cazy.org",
        "class_name",
        args_datasplit_family["args"],
        null_logger,
    )


# test get_subfamily_links


def test_get_subfam_links_len_0(family_h3_element, subfamily_urls, null_logger):
    """Test get_subfamily_links when no links are retrieved."""

    assert subfamily_urls == crawler.get_subfamily_links(
        family_h3_element,
        "http://www.cazy.org",
        null_logger,
    )


def test_get_subfam_links_urls(no_subfam_h3_element, null_logger):
    """Test get_subfamily_links when urls are retrieved."""

    assert None is crawler.get_subfamily_links(
        no_subfam_h3_element,
        "http://www.cazy.org",
        null_logger,
    )


# test parse_family()


def test_parse_family_fam_name(null_logger):
    """Tests parse_family() when family name is incorrectly formated"""

    crawler.parse_family(
        "http://www.cazy.org/GH_test_1.html",
        "GH_test_1",
        "http://www.cazy.org",
        null_logger,
    )


def test_parse_family_url(null_logger):
    """Test parse_family() whwn the url is incorrectly formatted."""

    crawler.parse_family(
        "http://www.cazy.org/GH_test_1.html",
        "GH1",
        "http://www.cazy.org",
        null_logger,
    )


def test_parse_family_no_page(null_logger, monkeypatch):
    """Tests parse_family() when no page retrieved for first paginiation page of proteins"""

    def mock_get_page(*args, **kwargs):
        return [None, "error message"]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    crawler.parse_family(
        "http://www.cazy.org/GH1.html",
        "GH1",
        "http://www.cazy.org",
        null_logger,
    )


def test_parse_family_success(protein_gen, null_logger, monkeypatch):
    """Test parse_family() when successful."""

    def mock_get_page(*args, **kwargs):
        return ["first paginiation page", None]

    def mock_get_urls(*args, **kwargs):
        return ["http://www.cazy.org/GH1_all.html"]

    def mock_parse_proteins(*args, **kwargs):
        return protein_gen

    monkeypatch.setattr(crawler, "get_page", mock_get_page)
    monkeypatch.setattr(crawler, "get_protein_page_urls", mock_get_urls)
    monkeypatch.setattr(crawler, "parse_proteins", mock_parse_proteins)

    crawler.parse_family(
        "http://www.cazy.org/GH1.html",
        "GH1",
        "http://www.cazy.org",
        null_logger,
    )


# test get_protein_page_urls()


def test_get_protein_page_urls_no_links(no_links_page):
    """Test get_protein_page_urls() on page with no links."""
    with open(no_links_page) as fp:
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


def test_parse_proteins_none(null_logger, monkeypatch):
    """Test parse_proteins when no page is returned."""

    def mock_get_page(*args, **kwargs):
        return [None, "error"]

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    assert True is isinstance(
        crawler.parse_proteins("protein_url", "family", null_logger),
        types.GeneratorType,
    )


def test_parse_proteins(gh147_page, null_logger, monkeypatch):
    """Tests parse_proteins when successful."""
    with open(gh147_page) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwargs):
        return soup

    monkeypatch.setattr(crawler, "get_page", mock_get_page)

    crawler.parse_proteins("protein_url", "family", null_logger)


# test row_to_protein()


def test_row_to_protein_no_ecs(protein_without_ec):
    """Test row_to_protein when no EC#s are listed."""
    with open(protein_without_ec) as fp:
        row = BeautifulSoup(fp, features="lxml")

    assert True is isinstance(
        crawler.row_to_protein(row, "GH147"),
        crawler.Protein
    )


def test_row_to_protein_ec(protein_with_ec):
    """Test row_to_protein when EC#s are listed."""
    with open(protein_with_ec) as fp:
        row = BeautifulSoup(fp, features="lxml")

    assert True is isinstance(
        crawler.row_to_protein(row, "GH147"),
        crawler.Protein
    )


def test_row_to_protein_gb_synoymns(protein_with_gb_synonyms):
    """Test row_to_protein when protein has GenBank synonyms."""
    with open(protein_with_gb_synonyms) as fp:
        row = BeautifulSoup(fp, features="lxml")

    assert True is isinstance(
        crawler.row_to_protein(row, "GH147"),
        crawler.Protein
    )


def test_row_to_protein_no_gb(protein_with_no_gb):
    """Test row_to_protein when protein has no GenBank accession."""
    with open(protein_with_no_gb) as fp:
        row = BeautifulSoup(fp, features="lxml")

    assert True is isinstance(
        crawler.row_to_protein(row, "GH147"),
        crawler.Protein
    )


# browser decorator and get_page


@pytest.mark.skip(reason="make trial testing quicker")
def test_browser_decorator():
    """Test browser_decorator to ensure proper handling if unsuccessful."""

    result = crawler.get_page('www.caz!!!!!!!!y.org')
    assert True == (result[0] is None) and (type(result[1]) is MissingSchema)
