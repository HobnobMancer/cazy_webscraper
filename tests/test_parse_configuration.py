#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License

"""Tests the submodule parse_configuration from utilities.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest
import sys

import pandas as pd

from argparse import Namespace

from scraper import sql
from scraper.utilities import parse_configuration


@pytest.fixture
def test_output_dir(test_dir):
    path = test_dir / "test_outputs" / "test_outputs_utilities"
    return path


@pytest.fixture
def making_output_dir(test_output_dir):
    path = test_output_dir / "testing_making_dir"
    return path


@pytest.fixture
def args_config_none():
    args_dict = {
        "args": Namespace(
            config=None,
            classes=None,
            families=None,
            genera=None,
            species=None,
            strains=None,
            kingdoms=None,
        )
    }
    return args_dict


@pytest.fixture
def parse_configuration_path():
    parse_configuration_path = parse_configuration.__file__
    return parse_configuration_path


@pytest.fixture
def args_config_fake(test_dir):
    path = test_dir / "fake_config.yaml"
    args_dict = {
        "args": Namespace(
            config=path
        )
    }
    return args_dict


@pytest.fixture
def input_dir(test_input_dir):
    path = test_input_dir / "test_inputs_parse_configuration"
    return path


@pytest.fixture
def config_file_path(input_dir):
    path = input_dir / "config_file.yaml"
    return path


@pytest.fixture
def args_config_file_cmd(config_file_path):
    args_dict = {
        "args": Namespace(
            config=config_file_path,
            classes="CE,AA",
            families="GH14,PL15",
            kingdoms="Bacteria",
        )
    }
    return args_dict


@pytest.fixture
def args_config_file(config_file_path):
    args_dict = {
        "args": Namespace(
            config=config_file_path,
            classes=None,
            families=None,
            kingdoms="Viruses,Archaea",
        )
    }
    return args_dict


@pytest.fixture
def args_config_cmd():
    args_dict = {
        "args": Namespace(
            config=None,
            classes="CE,AA",
            families="GH14,GT1,PL15,CE6,AA10,CBM50,DD21,CBMAA1",
            genera=None,
            species=None,
            strains=None,
            kingdoms="Archaea,Bacteria"
        )
    }
    return args_dict


@pytest.fixture
def empty_config_dict():
    config_dict = {
        'classes': [],
        'Glycoside Hydrolases (GHs)': [],
        'GlycosylTransferases (GTs)': [],
        'Polysaccharide Lyases (PLs)': [],
        'Carbohydrate Esterases (CEs)': [],
        'Auxiliary Activities (AAs)': [],
        'Carbohydrate-Binding Modules (CBMs)': [],
    }
    return config_dict


@pytest.fixture
def args_no_yaml():
    args_dict = {
        "args": Namespace(
            config="test/test/test_yaml.yaml",
            classes=None,
            families=None,
        )
    }
    return args_dict


@pytest.fixture
def testing_df():
    df_data = [["A", "B", "C"]]
    df = pd.DataFrame(df_data, columns=["C1", "C2", "C3"])
    return df


@pytest.fixture
def stdout_args(test_dir):
    args_dict = {
        "args": Namespace(
            output=sys.stdout,
        )
    }
    return args_dict


@pytest.fixture
def output_args(test_dir):
    path = test_dir / "test_outputs" / "test_outputs_parse_configuration"
    args_dict = {
        "args": Namespace(
            output=path,
        )
    }
    return args_dict


@pytest.fixture
def tax_no_yaml_args():
    args_dict = {
        "args": Namespace(
            config="fake_yaml",
            genera="Bacillus,Priestia",
            species="Arabidopsis thaliana",
            strains="Bathycoccus prasinos RCC1105,Cyanidioschyzon merolae strain 10D",
        )
    }
    return args_dict


@pytest.fixture
def tax_args(config_file_path):
    args_dict = {
        "args": Namespace(
            config=config_file_path,
            genera="Bacillus,Priestia",
            species="Arabidopsis thaliana",
            strains="Bathycoccus prasinos RCC1105,Cyanidioschyzon merolae strain 10D",
        )
    }
    return args_dict


@pytest.fixture
def raw_config_dict():
    raw_config_dict = {
        "classes": ["Carbohydrate Esterases (CEs)"],
        "Glycoside Hydrolases (GHs)": ["GH1"],
        "GlycosylTransferases (GTs)": ["GT1"],
        "Polysaccharide Lyases (PLs)": ["PL2", "PL3"],
        "Carbohydrate Esterases (CEs)": [],
        "Auxiliary Activities (AAs)": ["AA1"],
        "Carbohydrate-Binding Modules (CBMs)": ["CBM2"],
        "genera": ["Trichoderma"],
        "species": ["Aspergillus niger"],
        "strains": ["Acidianus ambivalens LEI 10"],
        "kingdoms": ["Archaea", "Bacteria"],
    }
    return raw_config_dict


# test parse_configuration()


def test_filenotfound_config(parse_configuration_path):
    """Test parse_config when the path to the configuraiton file is wrong."""
    args = {"args":Namespace(config="fake_page")}

    with pytest.raises(SystemExit) as pytest_wrapped_err:
        parse_configuration.parse_configuration(parse_configuration_path, args["args"])
    assert pytest_wrapped_err.type == SystemExit


def test_parse_config_file_cmd(
    args_config_file_cmd,
    cazy_dictionary,
    parse_configuration_path,
    monkeypatch,
):
    """Test parse_config when only config file is given."""

    std_classes = list(cazy_dictionary.keys())

    def mock_get_taxonomy_filter(*args, **kwargs):
        return None

    def mock_cazy_dict(*args, **kwargs):
        return cazy_dictionary, std_classes

    monkeypatch.setattr(parse_configuration, "get_cazy_dict_std_names", mock_cazy_dict)
    monkeypatch.setattr(parse_configuration, "get_genera_species_strains", mock_get_taxonomy_filter)

    (
        excluded_classes,
        config_dict,
        cazy_dictionary,
        tax_filter,
        kingdoms,
    ) = parse_configuration.parse_configuration(
        parse_configuration_path,
        args_config_file_cmd["args"],
    )

    expected_excluded_classes_1 = [
        '<strong>Carbohydrate-Binding Modules (CBMs)</strong>',
        '<strong>GlycosylTransferases (GTs)</strong>',
    ]

    expected_excluded_classes_2 = [
        '<strong>GlycosylTransferases (GTs)</strong>',
        '<strong>Carbohydrate-Binding Modules (CBMs)</strong>',
    ]

    expected_config_dict = {
        'classes': [
            'Glycoside Hydrolases (GHs)',
            'Polysaccharide Lyases (PLs)',
            'Carbohydrate Esterases (CEs)',
            'Auxiliary Activities (AAs)',
        ],
        'Glycoside Hydrolases (GHs)': ['GH147', 'GH14'],
        'GlycosylTransferases (GTs)': None,
        'Polysaccharide Lyases (PLs)': ['PL17', 'PL15'],
        'Carbohydrate Esterases (CEs)': None,
        'Auxiliary Activities (AAs)': None,
        'Carbohydrate-Binding Modules (CBMs)': None,
    }

    assert (expected_excluded_classes_1 == excluded_classes) or \
           (expected_excluded_classes_2 == excluded_classes)
    assert expected_config_dict == config_dict


def test_parse_config_file_only(args_config_file, cazy_dictionary, monkeypatch):
    """Test parse_configuration when only configuration file is given."""

    std_classes = list(cazy_dictionary.keys())

    def mock_cazy_dict(*args, **kwargs):
        return cazy_dictionary, std_classes

    def mock_get_tax_filter(*args, **kwargs):
        return

    monkeypatch.setattr(parse_configuration, "get_cazy_dict_std_names", mock_cazy_dict)
    monkeypatch.setattr(parse_configuration, "get_genera_species_strains", mock_get_tax_filter)

    (
        excluded_classes,
        config_dict,
        cazy_dictionary,
        tax_filter,
        kingdom,
    ) = parse_configuration.parse_configuration(
        parse_configuration_path,
        args_config_file["args"],
    )


def test_parse_config_cmd_only(args_config_cmd, cazy_dictionary, monkeypatch):
    """Test parse_configuration when only cmd-line configuration is given."""

    std_classes = list(cazy_dictionary.keys())

    def mock_cazy_dict(*args, **kwargs):
        return cazy_dictionary, std_classes

    monkeypatch.setattr(parse_configuration, "get_cazy_dict_std_names", mock_cazy_dict)

    (
        excluded_classes,
        config_dict,
        cazy_dictionary,
        tax_filter,
        kingdom,
    ) = parse_configuration.parse_configuration(
        parse_configuration_path,
        args_config_cmd["args"],
    )


def test_parse_config_no_config(args_config_none, cazy_dictionary, monkeypatch):
    """Test parse_configuration when no configuration is given."""

    std_classes = list(cazy_dictionary.keys())

    def mock_cazy_dict(*args, **kwargs):
        return cazy_dictionary, std_classes

    monkeypatch.setattr(parse_configuration, "get_cazy_dict_std_names", mock_cazy_dict)

    (
        excluded_classes,
        config_dict,
        cazy_dictionary,
        tax_filter,
        kingdom,
    ) = parse_configuration.parse_configuration(
        parse_configuration_path,
        args_config_none["args"],
    )


# test get_cazy_dict_std_names()


def test_cazy_dict_not_found():
    """Test get_cazy_dict_std_names() when the cazy_dictionary file cannot be found."""
    fake_path = sql.__file__

    with pytest.raises(SystemExit) as pytest_wrapped_err:
        parse_configuration.get_cazy_dict_std_names(fake_path)
    assert pytest_wrapped_err.type == SystemExit


def test_cazy_dict_success():
    """Test get_cazy_dict_std_names() when the retrieval of the data is successful."""
    path = parse_configuration.__file__

    parse_configuration.get_cazy_dict_std_names(path)


# test get_yaml_configuration


def test_get_yaml_not_found(empty_config_dict, cazy_dictionary, args_no_yaml):
    """Test get_yaml_configuration when the configuration cannot be found."""

    std_classes = list(cazy_dictionary.keys())

    with pytest.raises(SystemExit) as pytest_wrapped_err:
        parse_configuration.get_yaml_configuration(
            empty_config_dict,
            cazy_dictionary,
            std_classes,
            args_no_yaml["args"],
        )
    assert pytest_wrapped_err.type == SystemExit


def test_get_yaml_success(empty_config_dict, cazy_dictionary, args_config_file):
    """Test get_yaml_configuration when successfully retrieves and parses data from yaml file."""

    std_classes = list(cazy_dictionary.keys())

    parse_configuration.get_yaml_configuration(
        empty_config_dict,
        cazy_dictionary,
        std_classes,
        args_config_file["args"],
    )


# test get_yaml_cazy_classes


def test_get_yaml_classes_no(cazy_dictionary):
    """Test get_yaml_cazy_classes when classes key cannot be found."""

    std_classes = list(cazy_dictionary.keys())
    yaml_dict = {"Glycoside Hydrolases (GHs)": ["GH1", "GH2"]}

    parse_configuration.get_yaml_cazy_classes(yaml_dict, cazy_dictionary, std_classes)


def test_get_yaml_classes_none(cazy_dictionary):
    """Test get_yaml_cazy_classes when classes key has the value of None."""

    std_classes = list(cazy_dictionary.keys())
    yaml_dict = {"classes": None, "Glycoside Hydrolases (GHs)": ["GH1", "GH2"]}

    parse_configuration.get_yaml_cazy_classes(yaml_dict, cazy_dictionary, std_classes)


def test_get_yaml_classes_success(cazy_dictionary):
    """Test get_yaml_cazy_classes when classes key has the value of CAZy classes."""

    std_classes = list(cazy_dictionary.keys())
    yaml_dict = {"classes": ["PL", "CE"], "Glycoside Hydrolases (GHs)": ["GH1", "GH2"]}

    parse_configuration.get_yaml_cazy_classes(yaml_dict, cazy_dictionary, std_classes)


# test parse_user_cazy_classes()


def test_parse_user_cazy_classes(cazy_dictionary):
    """Test parse_user_cazy_classes to standardise the CAZy classes written by the user."""

    cazy_classes = ['GH', 'pl']
    class_name = [
        'Glycoside Hydrolases (GHs)',
        'GlycosylTransferases (GTs)',
        'Polysaccharide Lyases (PLs)',
        'Carbohydrate Esterases (CEs)',
        'Auxiliary Activities (AAs)',
        'Carbohydrate-Binding Modules (CBMs)',
    ]

    exepected = ['Glycoside Hydrolases (GHs)', 'Polysaccharide Lyases (PLs)']

    assert exepected == parse_configuration.parse_user_cazy_classes(
        cazy_classes,
        cazy_dictionary,
        class_name,
    )


def test_cannot_standardise_user_classes(cazy_dictionary):
    """Test when cannot standardise a name listed in the configuration file."""

    cazy_classes = ['GH', 'pl', 'testtesttest']
    class_name = [
        'Glycoside Hydrolases (GHs)',
        'GlycosylTransferases (GTs)',
        'Polysaccharide Lyases (PLs)',
        'Carbohydrate Esterases (CEs)',
        'Auxiliary Activities (AAs)',
        'Carbohydrate-Binding Modules (CBMs)',
    ]

    exepected = ['Glycoside Hydrolases (GHs)', 'Polysaccharide Lyases (PLs)']

    assert exepected == parse_configuration.parse_user_cazy_classes(
        cazy_classes,
        cazy_dictionary,
        class_name,
    )


# test get_cmd_defined_fams_classes()


def test_get_cmd_configuration(args_config_cmd, cazy_dictionary):
    """Test get_cmd_defined_fams_classes."""

    std_classes = list(cazy_dictionary.keys())
    parse_configuration.get_cmd_defined_fams_classes(
        cazy_dictionary,
        std_classes,
        args_config_cmd["args"],
    )


# test get_excluded_classes


def test_get_no_excluded_classes(cazy_dictionary):
    """Test when no excluded classes should be returned."""
    std_classes = list(cazy_dictionary.keys())
    config_dict = {"classes": [
        'Carbohydrate Esterases (CEs)',
        'Auxiliary Activities (AAs)',
        'GlycosylTransferases (GTs)',
        'Glycoside Hydrolases (GHs)',
        'Polysaccharide Lyases (PLs)',
        'Carbohydrate-Binding Modules (CBMs)',
    ]}

    assert None is parse_configuration.get_excluded_classes(
        std_classes,
        config_dict,
        cazy_dictionary,
    )

# test get_configuration() - retrieves configuraiton for the expand module


def test_get_configuraiton_0(cazy_dictionary, parse_configuration_path, monkeypatch):
    """Tests getting configuration for the expand module when NO tax filters are given."""

    def mock_parse_config(*args, **kwargs):
        return None, None, cazy_dictionary, {}, set()

    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_parse_config)

    parse_configuration.get_configuration(parse_configuration_path, "args")


def test_get_configuraiton_1(cazy_dictionary, parse_configuration_path, monkeypatch):
    """Tests getting configuration for the expand module when ARE tax filters are given."""

    def mock_parse_config(*args, **kwargs):
        return None, None, cazy_dictionary, {"genera": ["Trichoderma", "Aspergillus"]}, set(['viruses'])

    monkeypatch.setattr(parse_configuration, "parse_configuration", mock_parse_config)

    parse_configuration.get_configuration(parse_configuration_path, "args")


# test get_genera_species_strains()


def test_get_tax_filter_no_yaml(tax_no_yaml_args, raw_config_dict):
    """Test the get_genera_species_strains function when the config file can't be found."""
    parse_configuration.get_genera_species_strains(tax_no_yaml_args["args"], raw_config_dict)


def test_get_tax_filter(tax_args, raw_config_dict):
    """Test get_genera_species_strains when cmd-line and config file are parsed."""
    parse_configuration.get_genera_species_strains(tax_args["args"], raw_config_dict)


# test get_kingdoms


def test_get_kingdoms(raw_config_dict, args_config_cmd):
    """Test retrieving taxonomy Kingdoms."""
    raw_config_dict["kingdoms"] = ["archaea", "jackjackjack", "bacteria"]
    parse_configuration.get_kingdoms(args_config_cmd["args"], raw_config_dict)


def test_get_kingdoms_all_wrong(raw_config_dict):
    """Test retrieving taxonomy Kingdoms when none that were parsed are recognisable."""
    raw_config_dict["kingdoms"] = ["jack", "jill", "sam"]
    args = {"args": Namespace(kingdoms="milly")}

    with pytest.raises(SystemExit) as pytest_wrapped_err:
        parse_configuration.get_kingdoms(args["args"], raw_config_dict)
    assert pytest_wrapped_err.type == SystemExit


# test creating logger warning message when enabled streamlined scraping


def test_creating_streamlined_warning_message():
    """Test creating logger warning message when enabled streamlined scraping."""
    args = {"args": Namespace(streamline="genbank,ec,uniprot,pdb")}
    parse_configuration.create_streamline_scraping_warning(args["args"])
