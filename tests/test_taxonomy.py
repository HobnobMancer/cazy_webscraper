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
"""Tests the module taxonomy

These test are intened to be run from the root of the repository using:
pytest -v
"""


import logging
import pytest

from argparse import Namespace
from collections import namedtuple

from cazy_webscraper.ncbi.taxonomy import multiple_taxa


@pytest.fixture
def cazy_data():
    TaxData = namedtuple('Tax', ['kingdom', 'organism'])
    data = {
        "gbk1": {"kingdom": {'Bacteria'}, "taxonomy": {TaxData('king1', 'sp2'), TaxData('king2', 'sp2')}, "families": {'GH1', 'GH2'}},
        "gbk2": {"kingdom": {'Bacteria'}, "taxonomy": {TaxData('king2', 'sp2'), TaxData('king2', 'sp2')}, "families": {'GH1', 'GH2'}},
        "gbk3": {"kingdom": {'Bacteria'}, "taxonomy": {TaxData('king3', 'sp2'), TaxData('king2', 'sp2')}, "families": {'GH1', 'GH2'}},
        "gbk4": {"kingdom": {'Bacteria'}, "taxonomy": {TaxData('king1', 'sp2'), TaxData('king2', 'sp2')}, "families": {'GH1', 'GH2'}},
    }
    return data


def test_id_multi_tax(cazy_data):
    """Test identify_multiplpe_taxa"""
    logger = logging.getLogger(__name__)

    returned_list = multiple_taxa.identify_multiple_taxa(cazy_data, logger)

    assert len(returned_list) == 3


def test_select_first_organism(cazy_data):
    """Test select_first_organism()"""
    logger = logging.getLogger(__name__)

    multiple_taxa.select_first_organism(cazy_data, ['gbk1', 'gbk2', 'gbk3', 'gbk4'], logger)


def test_multi_taxa_invalid(cazy_data, monkeypatch):
    """Test replace_multiple_tax_with_invalid_ids"""
    gbk_accs = ['gbk1', 'gbk2', 'gbk3', 'gbk4']
    logger = logging.getLogger(__name__)

    def mock_replace_multiple(*args, **kwards):
        return cazy_data, False

    monkeypatch.setattr(multiple_taxa, "replace_multiple_tax", mock_replace_multiple)

    multiple_taxa.replace_multiple_tax_with_invalid_ids(cazy_data, gbk_accs, logger, "args")


def test_multi_taxa(cazy_data, monkeypatch):
    """Test replace_multiple_tax_with_invalid_ids"""
    gbk_accs = ['gbk1', 'gbk2', 'gbk3', 'gbk4']
    logger = logging.getLogger(__name__)

    def mock_replace_multiple(*args, **kwards):
        return cazy_data, True

    monkeypatch.setattr(multiple_taxa, "replace_multiple_tax", mock_replace_multiple)

    multiple_taxa.replace_multiple_tax_with_invalid_ids(cazy_data, gbk_accs, logger, "args")


def test_get_ncbi_tax(monkeypatch):
    argsdict = {"args": Namespace(
        retries=10,
        skip_ncbi_tax=True,
        ncbi_batch_size=200,
    )}

    entrez_result = "tests/test_inputs/test_inputs_ncbi_tax/entrezProt.xml"

    with open(entrez_result, "rb") as fh:
        result = fh

        def mock_entrez_tax_call(*args, **kwargs):
            """Mocks call to Entrez."""
            return result

        monkeypatch.setattr(multiple_taxa, "entrez_retry", mock_entrez_tax_call)

        output = multiple_taxa.get_ncbi_tax(
            {"WebEnv": 1, "QueryKey": 2},
            {'CAA35997.1': {'kingdom': {'kingdom'}, 'taxonomy': {'Bos taurus'}}},
            logging.getLogger(__name__),
            argsdict['args'],
        )
        assert output == {'CAA35997.1': {'kingdom': 'Eukaryota', 'taxonomy': {'Bos taurus'}, 'organism': 'Bos taurus'}}


def test_get_ncbi_tax_fails(monkeypatch):
    argsdict = {"args": Namespace(
        retries=10,
        skip_ncbi_tax=True,
        ncbi_batch_size=200,
    )}

    def mock_entrez_tax_call(*args, **kwargs):
        """Mocks call to Entrez."""
        return

    monkeypatch.setattr(multiple_taxa, "entrez_retry", mock_entrez_tax_call)

    multiple_taxa.get_ncbi_tax(
        {"WebEnv": 1, "QueryKey": 2},
        {'CAA35997.1': {'kingdom': {'kingdom'}, 'organism': {'organism'}}},
        logging.getLogger(__name__),
        argsdict['args'],
    )


def test_failed_replace_tax(monkeypatch):
    """Test replace_multiple_tax when failes to conenct to NCBI.Entrez"""
    argsdict = {"args": Namespace(
        retries=10,
        skip_ncbi_tax=True,
        ncbi_batch_size=200,
    )}

    def mock_return_none(*args, **kwards):
        return None

    def mock_select_org(*args, **kwards):
        return {}

    monkeypatch.setattr(multiple_taxa, "entrez_retry", mock_return_none)
    monkeypatch.setattr(multiple_taxa, "select_first_organism", mock_select_org)

    multiple_taxa.replace_multiple_tax({}, "gk1", logging.getLogger(__name__), argsdict["args"], [])


def test_failed_replace_tax_run_time(monkeypatch):
    """Test replace_multiple_tax when failes to conenct to NCBI.Entrez"""
    argsdict = {"args": Namespace(
        retries=10,
        skip_ncbi_tax=True,
        ncbi_batch_size=200,
    )}

    def mock_error(*args, **kwards):
        raise RuntimeError

    def mock_select_org(*args, **kwards):
        return {}

    def mock_replace_tax(*args, **kwards):
        return {}, False

    monkeypatch.setattr(multiple_taxa, "entrez_retry", mock_error)
    monkeypatch.setattr(multiple_taxa, "select_first_organism", mock_select_org)
    monkeypatch.setattr(multiple_taxa, "replace_multiple_tax_with_invalid_ids", mock_replace_tax)

    multiple_taxa.replace_multiple_tax({}, "gk1", logging.getLogger(__name__), argsdict["args"], [])
