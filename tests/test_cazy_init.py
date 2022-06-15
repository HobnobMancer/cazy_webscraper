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
from datetime import datetime
from pathlib import Path

import pytest

from sqlalchemy.exc import IntegrityError
from saintBioutils.utilities import logger as saint_logger
from saintBioutils.utilities import file_io as saint_fileIO

from cazy_webscraper import cazy

@pytest.fixture
def cazy_file_path():
    return "tests/test_inputs/test_inputs_cazy/cazy_data.txt"


@pytest.fixture
def cazy_zip_path():
    _path = "tests/test_inputs/test_inputs_cazy/cazy_db_timestamp.zip"
    return _path


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


def test_parse_all_cazy_data(cazy_file_path):
    with open(cazy_file_path, "r") as fh:
        cazy_lines = fh.read().splitlines()

    assert cazy.parse_all_cazy_data(cazy_lines, None) == {'UBD70155.1': {'kingdom': {'Bacteria'}, 'organism': {'Bacteroides cellulosilyticus BFG-250'}, 'families': {'GH157': {None}}}, 'ALJ59177.1': {'kingdom': {'Bacteria'}, 'organism': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}, 'WP_029429093.1': {'kingdom': {'Bacteria'}, 'organism': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}, 'BAX82587.1': {'kingdom': {'Bacteria'}, 'organism': {'Labilibaculum antarcticum SPP2'}, 'families': {'GH157': {None}}}, 'SDR68055.1': {'kingdom': {'Bacteria'}, 'organism': {'Polaribacter sp. KT25b'}, 'families': {'GH157': {None}}}, 'QXP79022.1': {'kingdom': {'Bacteria'}, 'organism': {'Winogradskyella sp. HaHa_3_26'}, 'families': {'GH157': {None}}}, 'QNK77973.1': {'kingdom': {'Bacteria'}, 'organism': {'Winogradskyella sp. PAMC22761'}, 'families': {'GH157': {None}}}, 'BCJ46567.1': {'kingdom': {'Bacteria'}, 'organism': {'Actinoplanes ianthinogenes NBRC 13996'}, 'families': {'AA5': {'AA5_2'}, 'GH1': {None}}}, 'ATO81963.1': {'kingdom': {'Bacteria'}, 'organism': {'Actinoplanes sp. SE50'}, 'families': {'AA5': {'AA5_2'}}}, 'AEV83893.1': {'kingdom': {'Bacteria'}, 'organism': {'Actinoplanes sp. SE50/110'}, 'families': {'AA5': {'AA5_2'}}}, 'SLL99371.1': {'kingdom': {'Bacteria'}, 'organism': {'Actinoplanes sp. SE50/110'}, 'families': {'AA5': {'AA5_2'}}}, 'CRK61066.1': {'kingdom': {'Bacteria'}, 'organism': {'Alloactinosynnema sp. L-07'}, 'families': {'AA5': {'AA5_2'}}}, 'AUZ88908.1': {'kingdom': {'Bacteria'}, 'organism': {'Arthrobacter agilis UMCV2'}, 'families': {'AA5': {'AA5_2'}}}, 'UKA55327.1': {'kingdom': {'Bacteria'}, 'organism': {'Arthrobacter sp. FW305-BF8'}, 'families': {'AA5': {'AA5_2'}}}, 'QMW46863.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus flavus AF13'}, 'families': {'AA5': {'AA5_2'}}}, 'QMW38907.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus flavus AF13'}, 'families': {'AA5': {'AA5_2'}}}, 'UCK59454.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus flavus CA14'}, 'families': {'AA5': {'AA5_2'}}}, 'UDD63818.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus flavus CA14'}, 'families': {'AA5': {'AA5_2'}}}, 'QRD87005.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus flavus NRRL 3357'}, 'families': {'AA5': {'AA5_2'}}}, 'QRD91111.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus flavus NRRL 3357'}, 'families': {'AA5': {'AA5_2'}}}, 'QMW26827.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus flavus NRRL3357'}, 'families': {'AA5': {'AA5_2'}}}, 'QMW34793.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus flavus NRRL3357'}, 'families': {'AA5': {'AA5_2'}}}, 'BAE64583.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus oryzae RIB40'}, 'families': {'AA5': {'AA5_2'}}}, 'BAE56565.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus oryzae RIB40'}, 'families': {'AA5': {'AA5_2'}}}, 'BCS28434.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus puulaauensis MK2'}, 'families': {'AA5': {'AA5_2'}}}, 'BCS18659.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus puulaauensis MK2'}, 'families': {'AA5': {'AA5_2'}}}, 'BAW27603.1': {'kingdom': {'Eukaryota'}, 'organism': {'Aspergillus stellatus NBRC 32302'}, 'families': {'AA5': {'AA5_2'}}}}


def test_parse_cazy_data(cazy_file_path):
    with open(cazy_file_path, "r") as fh:
        cazy_lines = fh.read().splitlines()

    assert cazy.parse_cazy_data_with_filters(
        cazy_lines,
        {'GH'},
        {'Bacteria'},
        {'GH157', 'AA5_2'},
        {'Bacteroides', 'Aspergillus oryzae', 'Aspergillus flavus NRRL3357'},
        None,
    ) == {'UBD70155.1': {'kingdom': {'Bacteria'}, 'organism': {'Bacteroides cellulosilyticus BFG-250'}, 'families': {'GH157': {None}}}, 'ALJ59177.1': {'kingdom': {'Bacteria'}, 'organism': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}, 'WP_029429093.1': {'kingdom': {'Bacteria'}, 'organism': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}}


def test_build_tax_dict():
    cazy_data = {'UBD70155.1': {'kingdom': {'Bacteria'}, 'organism': {'Bacteroides cellulosilyticus BFG-250'}, 'families': {'GH157': {None}}}, 'ALJ59177.1': {'kingdom': {'Bacteria'}, 'organism': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}, 'WP_029429093.1': {'kingdom': {'Bacteria'}, 'organism': {'Bacteroides cellulosilyticus WH2'}, 'families': {'GH157': {None}}}}
    assert {
        'Bacteria': {
            'Bacteroides cellulosilyticus BFG-250', 'Bacteroides cellulosilyticus WH2',
        }} == cazy.build_taxa_dict(cazy_data)
