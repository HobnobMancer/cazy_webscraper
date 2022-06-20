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
"""Tests expand.taxonomy.get_ncbi_taxs.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


from argparse import Namespace

from Bio import Entrez

from cazy_webscraper.ncbi import genomes


def test_get_nuccore_ids(monkeypatch):
    """Retrieve mock result with https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?db=nuccore&dbfrom=protein&id=BCS34995.1&linkname=protein_nuccore"""
    argsdict = {"args": Namespace(
        retries=10,
    )}

    efetch_result = "tests/test_inputs/test_inputs_ncbi_genomes/elink_prot_nuccore.xml"

    def mock_post_ids(*args, **kwards):
        return "qk", "we"

    with open(efetch_result, "rb") as fh:
        result = fh

        def mock_entrez_tax_call(*args, **kwargs):
            """Mocks call to Entrez to retrieve nuccore ids."""
            return result

        monkeypatch.setattr(genomes, "entrez_retry", mock_entrez_tax_call)
        monkeypatch.setattr(genomes, "post_ids", mock_post_ids)

        output = genomes.get_nuccore_ids(['BCS34995.1'], {}, argsdict['args'])
        assert output == ({'2131947417'}, {})


def test_get_nuccore_ids_fail(monkeypatch):
    argsdict = {"args": Namespace(
        retries=10,
    )}

    def mock_post_ids(*args, **kwards):
        return "qk", "we"

    def mock_entrez_tax_call(*args, **kwargs):
        """Mocks call to Entrez to retrieve nuccore ids."""
        return

    monkeypatch.setattr(genomes, "entrez_retry", mock_entrez_tax_call)
    monkeypatch.setattr(genomes, "post_ids", mock_post_ids)

    genomes.get_nuccore_ids(['BCS34995.1'], {}, argsdict['args'])
