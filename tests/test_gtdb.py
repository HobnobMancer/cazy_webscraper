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
"""Tests the get_gtdb_taxs.py from the expand module.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest
from pathlib import Path

from cazy_webscraper.expand.gtdb import get_gtdb_tax


@pytest.fixture
def archaea_data_file_path():
    _path = Path("tests/test_inputs/test_inputs_gtdb/ar53_taxonomy.tsv.gz")
    return _path

@pytest.fixture
def test_genomes():
    _genomes = [
        'GCF_000979375.1',
        'GCF_11111111.1',
        'GCA_002506415.1',
    ]
    return _genomes


def test_parse_gtdb_datafile(archaea_data_file_path, test_genomes):
    result = get_gtdb_tax.get_lineage_data(archaea_data_file_path, test_genomes)
    assert result == {
        'GCF_000979375.1': {'lineage': 'd__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei', 'release': 'ar53_taxonomy.tsv'},
        'GCA_002506415.1': {'lineage': 'd__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei', 'release': 'ar53_taxonomy.tsv'}
    }
