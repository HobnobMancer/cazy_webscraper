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
"""Tests adding cazyme data via the sql.sql_interface.add_data module

These test are intened to be run from the root of the repository using:
pytest -v
"""


import json
import pandas as pd

from argparse import Namespace, ArgumentParser
from datetime import datetime
from pathlib import Path

import pytest

from sqlalchemy.exc import IntegrityError
from saintBioutils.utilities import logger as saint_logger
from saintBioutils.utilities import file_io as saint_fileIO

from cazy_webscraper.sql.sql_interface.add_data import add_cazyme_data


def test_add_kingdoms(monkeypatch):
    cazy_tax_data = {'Bacteria': {}, 'Viruses': {}}
    connection = None

    def mock_return_none(*args, **kwards):
        return

    def mock_kingdom_table(*args, **kwards):
        return {}

    monkeypatch.setattr(add_cazyme_data, "insert_data", mock_return_none)
    monkeypatch.setattr(add_cazyme_data, "get_kingdom_table_dict", mock_kingdom_table)

    add_cazyme_data.add_kingdoms(cazy_tax_data, connection)


def test_add_organisms(monkeypatch):
    tax_dict = {
        'Bacteria': {'Genus species'},
        'Eukaryota': {'Genus1 species1 strain1'},
    }
    connection = None

    def mock_kingdom_table(*args, **kwards):
        return {'Bacteria': 1, 'Eukaryota': 2}
    
    def mock_tax_table(*args, **kwards):
        return {}

    def mock_return_none(*args, **kwards):
        return

    monkeypatch.setattr(add_cazyme_data, "get_kingdom_table_dict", mock_kingdom_table)
    monkeypatch.setattr(add_cazyme_data, "get_taxs_table_dict", mock_tax_table)
    monkeypatch.setattr(add_cazyme_data, "insert_data", mock_return_none)

    add_cazyme_data.add_source_organisms(tax_dict, connection)


def test_add_fams(monkeypatch):
    cazy_data = {
        'gbk1': {'kingdom': 'k', 'organism': 'o', 'families': {
                'fam': {'subfam'},
                'fam1': {None},
            }
        }
    }
    connection = None

    def mock_fam_table(*args, **kwards):
        return {}

    def mock_return_none(*args, **kwards):
        return

    monkeypatch.setattr(add_cazyme_data, "get_fams_table_dict", mock_fam_table)
    monkeypatch.setattr(add_cazyme_data, "insert_data", mock_return_none)

    add_cazyme_data.add_cazy_families(cazy_data, connection)


def test_add_gbks(monkeypatch):

    def mock_gbk_table(*args, **kwards):
        return {}
    
    def mock_tax_table(*args, **kwards):
        return {
            'genus species': {'tax_id': 1, 'kingdom_id': 1},
        }

    def mock_return_none(*args, **kwards):
        return

    cazy_data = {
        'gbk1': {'kingdom': 'k', 'organism': 'genus species', 'families': {
                'fam': {'subfam'},
                'fam1': {None},
            }
        }
    }
    connection = None

    monkeypatch.setattr(add_cazyme_data, "get_gbk_table_dict", mock_gbk_table)
    monkeypatch.setattr(add_cazyme_data, "get_taxs_table_dict", mock_tax_table)
    monkeypatch.setattr(add_cazyme_data, "insert_data", mock_return_none)

    add_cazyme_data.add_genbanks(cazy_data, connection)


def test_add_gbks_fam_rel(monkeypatch):

    def mock_gbk_table(*args, **kwards):
        return {'gbk1': {'tax_id': 1, 'gbk_id': 1}}

    def mock_fam_table(*args, **kwards):
        return {'fam subfam': 1, 'fam1 _': 2}

    def mock_gbk_fam_table(*args, **kwards):
        return {}, []

    def mock_return_none(*args, **kwards):
        return

    cazy_data = {
        'gbk1': {'kingdom': 'k', 'organism': ('genus species',), 'families': {
                'fam': {'subfam'},
                'fam1': {None},
            }
        }
    }
    connection = None

    monkeypatch.setattr(add_cazyme_data, "get_gbk_table_dict", mock_gbk_table)
    monkeypatch.setattr(add_cazyme_data, "get_fams_table_dict", mock_fam_table)
    monkeypatch.setattr(add_cazyme_data, "get_gbk_fam_table_dict", mock_gbk_fam_table)
    monkeypatch.setattr(add_cazyme_data, "insert_data", mock_return_none)

    add_cazyme_data.add_genbank_fam_relationships(cazy_data, connection, 'args')
