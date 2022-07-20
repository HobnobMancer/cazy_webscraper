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
"""Tests the module utilities.parse_configuration which parses the users configuration.

These test are intened to be run from the root of the repository using:
pytest -v
"""

from distutils.command.config import config
from cazy_webscraper.utilities import parse_configuration


import pytest

from argparse import Namespace


@pytest.fixture
def config_dict_blank():
    config_dict = {
        "classes": set(),
        "Glycoside Hydrolases (GHs)": set(),
        "GlycosylTransferases (GTs)": set(),
        "Polysaccharide Lyases (PLs)": set(),
        "Carbohydrate Esterases (CEs)": set(),
        "Auxiliary Activities (AAs)": set(),
        "Carbohydrate-Binding Modules (CBMs)": set(),
    }
    return config_dict


@pytest.fixture
def yaml_path():
    path_ = "tests/test_inputs/test_inputs_parse_configuration/config_file.yaml"
    return path_


@pytest.fixture
def tax_dict():
    tax_dict = {'genera': set(), 'species': set(), 'strains': set()}
    return tax_dict


def test_class_synomymns_error():
    args_dict = {
        "args": Namespace(
            cazy_synonyms ="testtesttest",
        )
    }

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_configuration.get_cazy_class_synonym_dict(args_dict['args'])
    assert pytest_wrapped_e.type == SystemExit


def test_class_synonyms():
    args_dict = {
        "args": Namespace(
            cazy_synonyms="tests/test_inputs/cazy_dictionary.json",
        )
    }
    parse_configuration.get_cazy_class_synonym_dict(args_dict['args'])


def test_yaml_config_no_yaml(config_dict_blank, cazy_dictionary, tax_dict):
    args_dict = {
        "args": Namespace(
            config="path",
        )
    }
    
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_configuration.get_yaml_configuration(
            config_dict_blank,
            cazy_dictionary,
            tax_dict,
            set(),
            args_dict['args'],
        )
    assert pytest_wrapped_e.type == SystemExit


def test_yaml_config(config_dict_blank, cazy_dictionary, tax_dict, monkeypatch, yaml_path):
    def mock_get_yaml_classes(*args, **kwards):
        return ['GH', 'PL']
    
    args_dict = {
        "args": Namespace(
            config=yaml_path,
        )
    }
    
    monkeypatch.setattr(parse_configuration, "get_yaml_cazy_classes", mock_get_yaml_classes)
    
    parse_configuration.get_yaml_configuration(
        config_dict_blank,
        cazy_dictionary,
        tax_dict,
        set(),
        args_dict['args'],
    )


def test_yaml_classes_keyerror(cazy_dictionary):
    config_dict = {}
    parse_configuration.get_yaml_cazy_classes(config_dict, cazy_dictionary)


def test_yaml_classes_none(cazy_dictionary):
    config_dict = {'classes': None}
    parse_configuration.get_yaml_cazy_classes(config_dict, cazy_dictionary)


def test_yaml_classes(cazy_dictionary, monkeypatch):
    def mock_get_classes(*args, **kwards):
        return ["GH", "PL"]

    config_dict = {'classes': ['PL']}

    monkeypatch.setattr(parse_configuration, "parse_user_cazy_classes", mock_get_classes)

    parse_configuration.get_yaml_cazy_classes(config_dict, cazy_dictionary)


def test_parse_classes(cazy_dictionary):
    parse_configuration.parse_user_cazy_classes(["GH", "Cat", "gh"], cazy_dictionary)


def test_parse_cmd_tax(monkeypatch, cazy_dictionary, config_dict_blank):

    def mock_get_classes(*args, **kwards):
        return ["GH", "PL"]

    def mock_get_fams(*args, **kwards):
        return config_dict_blank

    args_dict = {
        "args": Namespace(
            classes="GH,PL",
            families="GH1,GT2,PL3,CE19,AA10,CBM4,QQQ444,CBMQQQ",
            genera="Aspergillus,Trichoderma",
            species="Aspergillus niger,Genus species",
            strains="Genus species strain,genus species strain",
            kingdoms="bacteria,unclassified",
        )
    }

    tax_dict = {'genera': set(), 'species': set(), 'strains': set()}

    monkeypatch.setattr(parse_configuration, "parse_user_cazy_classes", mock_get_classes)
    monkeypatch.setattr(parse_configuration, "get_cmd_defined_families", mock_get_fams)

    parse_configuration.get_cmd_scrape_config(config_dict_blank, cazy_dictionary, tax_dict, set(), args_dict["args"])


def test_parse_cmd_fams(config_dict_blank):
    args_dict = {
        "args": Namespace(
            families="GH1,GT2,PL3,CE19,AA10,CBM4,QQQ444,CBMQQQ"
        )
    }

    parse_configuration.get_cmd_defined_families(config_dict_blank, args_dict["args"])


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


def test_ec_filters(yaml_path):
    args_dict = {
        "args": Namespace(
            config=yaml_path,
        )
    }
    ec_filters = "EC3.2.1.3,1.-.-.-,ec10.0.1.*"
    parse_configuration.get_ec_config(ec_filters, args_dict['args'])


def test_ec_err():
    args_dict = {
        "args": Namespace(
            config="path",
        )
    }
    ec_filters = "EC3.2.1.3,1.-.-.-,ec10.0.1.*"

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_configuration.get_ec_config(ec_filters, args_dict['args'])
    assert pytest_wrapped_e.type == SystemExit
    