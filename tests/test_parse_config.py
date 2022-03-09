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
"""Tests the module utilities.parse_configuration which parses the users configuration.

These test are intened to be run from the root of the repository using:
pytest -v
"""

from cazy_webscraper.utilities import parse_configuration


import pytest

from argparse import Namespace


def test_parse_cmd_fams():
    args_dict = {
        "args": Namespace(
            families="GH1,GT2,PL3,CE19,AA10,CBM4,QQQ444,CBMQQQ"
        )
    }

    config_dict = {
        "Glycoside Hydrolases (GHs)": set(),
        "GlycosylTransferases (GTs)": set(),
        "Polysaccharide Lyases (PLs)": set(),
        "Carbohydrate Esterases (CEs)": set(),
        "Auxiliary Activities (AAs)": set(),
        "Carbohydrate-Binding Modules (CBMs)": set(),
    }

    parse_configuration.get_cmd_defined_families(config_dict, args_dict["args"])


def test_excluded_classes(cazy_dictionary):
    config_dict = {
        "classes": set(),
        "Glycoside Hydrolases (GHs)": ["GH1", "GH2"],
    }

    parse_configuration.get_excluded_classes(config_dict, cazy_dictionary)


def test_excluded_classes_none(cazy_dictionary):
    config_dict = cazy_dictionary
    config_dict["classes"] = set()

    parse_configuration.get_excluded_classes(config_dict, cazy_dictionary)


def test_convert_empty():
    tax_dict = {'genera': {'Aspergillus', 'Trichoderma'}, 'species': set()}

    parse_configuration.convert_empty_sets_to_none(tax_dict)


def test_get_filter_set():
    tax_dict = {'genera': {'Aspergillus', 'Trichoderma'}, 'species': None}

    parse_configuration.get_filter_set(tax_dict)
