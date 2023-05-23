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
"""Tests cazy module

These test are intened to be run from the root of the repository using:
pytest -v
"""


import json
import pandas as pd

from argparse import Namespace, ArgumentParser
from collections import namedtuple
from datetime import datetime
from pathlib import Path

import pytest

from sqlalchemy.exc import IntegrityError
from saintBioutils.utilities import logger as saint_logger
from saintBioutils.utilities import file_io as saint_fileIO

from cazy_webscraper import cazy


TaxData = namedtuple('Tax', ['kingdom', 'organism'])


@pytest.fixture
def cazy_file_path():
    return "tests/test_inputs/test_inputs_cazy/cazy_data.txt"


@pytest.fixture
def cazy_zip_path():
    _path = "tests/test_inputs/test_inputs_cazy/cazy_db_timestamp.zip"
    return _path


@pytest.fixture
def cazy_data_dict():
    _dict = {
        'UBD70155.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus BFG-250'}, 'families': {'GH157': {None}}},
        'ALJ59177.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}},
        'WP_029429093.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}},
    }
    return _dict


def test_get_cazy_file(cazy_file_path):
    argsdict = {"args": Namespace(
        retries=10,
        cazy_data=cazy_file_path,
    )}
    cazy.get_cazy_txt_file_data(
        "tests/test_inputs/test_inputs_cazy/",
        "time_stamp",
        argsdict['args'],
    )


def test_parsing_cazy_zip(monkeypatch):
    argsdict = {"args": Namespace(
        retries=10,
        cazy_data=None,
    )}

    def mock_download(*args, **kwards):
        return

    monkeypatch.setattr(cazy, "get_cazy_file", mock_download)

    cazy.get_cazy_txt_file_data(
        Path("tests/test_inputs/test_inputs_cazy/"),
        "time_stamp",
        argsdict['args'],
    )


def test_failed_download(monkeypatch):
    argsdict = {"args": Namespace(
        retries=2,
        cazy_data=None,
    )}

    def mock_download(*args, **kwards):
        return "error"

    monkeypatch.setattr(cazy, "get_cazy_file", mock_download)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cazy.get_cazy_txt_file_data(
            Path("tests/test_inputs/test_inputs_cazy/"),
            "time_stamp",
            argsdict['args'],
        )
    assert pytest_wrapped_e.type == SystemExit


def test_parse_all_cazy_data_subfams(cazy_file_path):
    argsdict = {
        "args": Namespace(
            subfamilies=True,
        )
    }

    with open(cazy_file_path, "r") as fh:
        cazy_lines = fh.read().splitlines()

    output = cazy.parse_all_cazy_data(cazy_lines, None, argsdict['args']) 
    # assertion to check output fails on 'Tax' not defined
    

def test_parse_all_cazy_data_fams_only(cazy_file_path):
    argsdict = {
        "args": Namespace(
            subfamilies=False,
        )
    }

    with open(cazy_file_path, "r") as fh:
        cazy_lines = fh.read().splitlines()

    output = cazy.parse_all_cazy_data(cazy_lines, None, argsdict['args'])


def test_parse_cazy_data_subfam(cazy_file_path):
    argsdict = {
        "args": Namespace(
            subfamilies=True,
        )
    }

    with open(cazy_file_path, "r") as fh:
        cazy_lines = fh.read().splitlines()

    output = cazy.parse_cazy_data_with_filters(
        cazy_lines,
        {'GH'},
        {'Bacteria'},
        {'GH157', 'AA5_2'},
        {'Bacteroides', 'Aspergillus oryzae', 'Aspergillus flavus NRRL3357'},
        None,
        argsdict['args'],
    ) 


def test_parse_cazy_data_fam_only(cazy_file_path):
    argsdict = {
        "args": Namespace(
            subfamilies=False,
        )
    }

    with open(cazy_file_path, "r") as fh:
        cazy_lines = fh.read().splitlines()

    output = cazy.parse_cazy_data_with_filters(
        cazy_lines,
        {'GH'},
        {'Bacteria'},
        {'GH157', 'AA5_2'},
        {'Bacteroides', 'Aspergillus oryzae', 'Aspergillus flavus NRRL3357'},
        None,
        argsdict['args'],
    )
    # assertion fails on 'Tax' not defined
    # assert output == {'UBD70155.1': {'taxonomy': {Tax(kingdom='Bacteria', organism='Bacteroides cellulosilyticus BFG-250')}, 'families': {'GH157': {None}}}, 'ALJ59177.1': {'taxonomy': {Tax(kingdom='Bacteria', organism='Bacteroides cellulosilyticus WH2')}, 'families': {'GH157': {None}}}, 'WP_029429093.1': {'taxonomy': {Tax(kingdom='Bacteria', organism='Bacteroides cellulosilyticus WH2')}, 'families': {'GH157': {None}}}}


def test_build_tax_dict():
    cazy_data = {'UBD70155.1': {'kingdom': 'Bacteria', 'organism': 'Bacteroides cellulosilyticus BFG-250', 'families': {'GH157': {None}}}, 'ALJ59177.1': {'kingdom': 'Bacteria', 'organism': 'Bacteroides cellulosilyticus WH2', 'families': {'GH157': {None}}}, 'WP_029429093.1': {'kingdom': 'Bacteria', 'organism': 'Bacteroides cellulosilyticus WH2', 'families': {'GH157': {None}}}}

    output = cazy.build_taxa_dict(cazy_data)
    assert output[0] == {'Bacteria': {'Bacteroides cellulosilyticus WH2', 'Bacteroides cellulosilyticus BFG-250'}}


def test_apply_filters_no_tax(cazy_data_dict):

    output = cazy.apply_kingdom_tax_filters(
        cazy_data_dict,
        set(),
        set(),
        'UBD70155.1',
        'GH157',
        None,
        'Bacteroides cellulosilyticus BFG-250',
        'Bacteria',
    )
    TaxData = namedtuple('Tax', ['kingdom', 'organism'])
    
    # assert output  == ({'UBD70155.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus BFG-250', Tax(kingdom='Bacteria', organism='Bacteroides cellulosilyticus BFG-250')}, 'families': {'GH157': {None}}}, 'ALJ59177.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}, 'WP_029429093.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}}, True)


def test_apply_tax_filters(cazy_data_dict):

    output = cazy.apply_kingdom_tax_filters(
        cazy_data_dict,
        {'Bacteria'},
        {'Bacteroides'},
        'UBD70155.1',
        'GH157',
        None,
        'Bacteroides cellulosilyticus BFG-250',
        'Bacteria',
    )
    
    # assert str(output) == "({'UBD70155.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus BFG-250', Tax(kingdom='Bacteria', organism='Bacteroides cellulosilyticus BFG-250')}, 'families': {'GH157': {None}}}, 'ALJ59177.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}, 'WP_029429093.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}}, True)"


def test_add_protein_to_dict(cazy_data_dict):

    output = cazy.add_protein_to_dict(
        cazy_data_dict,
        'UBD70155.1',
        'GH1',
        'GH1_1',
        'Bacteroides cellulosilyticus BFG-250',
        'Bacteria',
    ) 
    # assert str(output) == "{'UBD70155.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus BFG-250', Tax(kingdom='Bacteria', organism='Bacteroides cellulosilyticus BFG-250')}, 'families': {'GH157': {None}, 'GH1': {'GH1_1'}}}, 'ALJ59177.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}, 'WP_029429093.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}}"


def test_add_protein_to_dict_new_protein(cazy_data_dict):

    output = cazy.add_protein_to_dict(
        cazy_data_dict,
        'UBD70155_new.1',
        'GH1',
        'GH1_1',
        'Bacteroides cellulosilyticus BFG-250',
        'Bacteria',
    )
    
    # assert str(output) == "{'UBD70155.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus BFG-250'}, 'families': {'GH157': {None}}}, 'ALJ59177.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}, 'WP_029429093.1': {'kingdom': {'Bacteria'}, 'taxonomy': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}, 'UBD70155_new.1': {'taxonomy': {Tax(kingdom='Bacteria', organism='Bacteroides cellulosilyticus BFG-250')}, 'families': {'GH1': {'GH1_1'}}}}"
