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
"""Tests the __init__.py from the expand module.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace

from cazy_webscraper import ncbi


def test_epost(monkeypatch):
    """Get unit test input from https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=pubmed&id=11237011,12466850"""
    argsdict = {"args": Namespace(
        retries=10,
    )}

    efetch_result = "tests/test_inputs/test_inputs_ncbi_tax/epost.xml"

    with open(efetch_result, "rb") as fh:
        result = fh

        def mock_entrez_tax_call(*args, **kwargs):
            """Mocks call to Entrez."""
            return result

        monkeypatch.setattr(ncbi, "entrez_retry", mock_entrez_tax_call)

        output = ncbi.post_ids(['11237011', '12466850'], 'Pubmed', argsdict['args'])
        assert len(output) == 2


def test_epost_str(monkeypatch):
    """Get unit test input from https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=pubmed&id=11237011,12466850"""
    argsdict = {"args": Namespace(
        retries=10,
    )}

    efetch_result = "tests/test_inputs/test_inputs_ncbi_tax/epost.xml"

    with open(efetch_result, "rb") as fh:
        result = fh

        def mock_entrez_tax_call(*args, **kwargs):
            """Mocks call to Entrez."""
            return result

        monkeypatch.setattr(ncbi, "entrez_retry", mock_entrez_tax_call)

        output = ncbi.post_ids('11237011', 'Pubmed', argsdict['args'])
        assert len(output) == 2


def test_epost_fail(monkeypatch):
    """Get unit test input from https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=pubmed&id=11237011,12466850"""
    argsdict = {"args": Namespace(
        retries=10,
    )}

    def mock_entrez_tax_call(*args, **kwargs):
        """Mocks call to Entrez."""
        return

    monkeypatch.setattr(ncbi, "entrez_retry", mock_entrez_tax_call)

    output = ncbi.post_ids('11237011', 'Pubmed', argsdict['args'])
    assert len(output) == 2
