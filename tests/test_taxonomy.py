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
"""Tests the module taxonomy

These test are intened to be run from the root of the repository using:
pytest -v
"""


import logging
import pytest

from cazy_webscraper import taxonomy


@pytest.fixture
def cazy_data():
    data = {
        "gbk1": {"kingdom": {'Bacteria'}, "organism": {"sp 1", "sp 2"}, "families": {'GH1', 'GH2'}},
        "gbk2": {"kingdom": {'Bacteria'}, "organism": {"sp 1", "sp 2"}, "families": {'GH1', 'GH2'}},
        "gbk3": {"kingdom": {'Bacteria'}, "organism": {"sp 1", "sp 2"}, "families": {'GH1', 'GH2'}},
        "gbk4": {"kingdom": {'Bacteria'}, "organism": {"sp 1", "sp 2"}, "families": {'GH1', 'GH2'}},
    }
    return data


def test_id_multi_tax(cazy_data):
    """Test identify_multiplpe_taxa"""
    logger = logging.getLogger(__name__)

    returned_list = taxonomy.identify_multiple_taxa(cazy_data, logger)

    assert len(returned_list) == 4


def test_select_first_organism(cazy_data):
    """Test select_first_organism()"""
    logger = logging.getLogger(__name__)

    taxonomy.select_first_organism(cazy_data, ['gbk1', 'gbk2', 'gbk3', 'gbk4'], logger)


def test_multi_taxa_invalid(cazy_data, monkeypatch):
    """Test replace_multiple_tax_with_invalid_ids"""
    gbk_accs = ['gbk1', 'gbk2', 'gbk3', 'gbk4']
    logger = logging.getLogger(__name__)

    def mock_replace_multiple(*args, **kwards):
        return cazy_data, False
    
    monkeypatch.setattr(taxonomy, "replace_multiple_tax", mock_replace_multiple)

    taxonomy.replace_multiple_tax_with_invalid_ids(cazy_data, gbk_accs, logger, "args")


def test_multi_taxa(cazy_data, monkeypatch):
    """Test replace_multiple_tax_with_invalid_ids"""
    gbk_accs = ['gbk1', 'gbk2', 'gbk3', 'gbk4']
    logger = logging.getLogger(__name__)

    def mock_replace_multiple(*args, **kwards):
        return cazy_data, True
    
    monkeypatch.setattr(taxonomy, "replace_multiple_tax", mock_replace_multiple)

    taxonomy.replace_multiple_tax_with_invalid_ids(cazy_data, gbk_accs, logger, "args")
    