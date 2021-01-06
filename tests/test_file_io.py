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

"""Tests the module file_io which builds the logger and args parser.

These test are intened to be run from the root of the repository using:
pytest -v

"""

import json
import pytest
import sys

import pandas as pd

from argparse import Namespace

from scraper import file_io, parse


@pytest.fixture
def test_output_dir(test_dir):
    path = test_dir / "test_outputs" / "test_outputs_file_io"
    return path


@pytest.fixture
def making_output_dir(test_output_dir):
    path = test_output_dir / "testing_making_dir"
    return path


@pytest.fixture
def args_config_None():
    args_dict = {
        "args": Namespace(
            config=None,
            classes=None,
            families=None,
        )
    }
    return args_dict


@pytest.fixture
def file_io_path():
    file_io_path = file_io.__file__
    return file_io_path


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
    path = test_input_dir / "test_inputs_file_io"
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
        )
    }
    return args_dict


@pytest.fixture
def args_config_cmd(config_file_path):
    args_dict = {
        "args": Namespace(
            config=None,
            classes="CE,AA",
            families="GH14,PL15",
        )
    }
    return args_dict


@pytest.fixture
def args_config_none(config_file_path):
    args_dict = {
        "args": Namespace(
            config=None,
            classes=None,
            families=None,
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
def df_output_file(test_dir):
    df_output = (
        test_dir / "test_targets" / "file_io_test_targets" / "test_writing_df.csv"
    )
    return df_output


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
    path = test_dir / "test_outputs" / "test_outputs_file_io"
    args_dict = {
        "args": Namespace(
            output=path,
        )
    }
    return args_dict


# test make_output_directory()


def test_output_dir_creation_nd_true(making_output_dir, null_logger):
    """Test creation of output dir when nodelete is false"""
    file_io.make_output_directory(making_output_dir, null_logger, True, False)


def test_output_dir_creation_nd_false(making_output_dir, null_logger):
    """Test creation of output dir when nodelete is true"""
    file_io.make_output_directory(making_output_dir, null_logger, True, True)


def test_output_dir_creation_exists(test_dir, null_logger):
    """Test creation of output dir when already exists."""
    file_io.make_output_directory(test_dir, null_logger, False, False)


# test parse_configuration()


def test_parse_config_file_cmd(
    args_config_file_cmd,
    cazy_dictionary,
    file_io_path,
    null_logger,
    monkeypatch,
):
    """Test parse_config when only config file is given."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    std_classes = list(cazy_dict.keys())

    def mock_cazy_dict(*args, **kwargs):
        return cazy_dict, std_classes

    monkeypatch.setattr(file_io, "get_cazy_dict_std_names", mock_cazy_dict)

    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args_config_file_cmd["args"],
        null_logger,
    )

    expected_excluded_classes = [
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

    assert expected_excluded_classes == excluded_classes
    assert expected_config_dict == config_dict


def test_parse_config_file_only(args_config_file, cazy_dictionary, null_logger, monkeypatch):
    """Test parse_configuration when only configuration file is given."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    std_classes = list(cazy_dict.keys())

    def mock_cazy_dict(*args, **kwargs):
        return cazy_dict, std_classes

    monkeypatch.setattr(file_io, "get_cazy_dict_std_names", mock_cazy_dict)

    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args_config_file["args"],
        null_logger,
    )


def test_parse_config_cmd_only(args_config_cmd, cazy_dictionary, null_logger, monkeypatch):
    """Test parse_configuration when only cmd-line configuration is given."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    std_classes = list(cazy_dict.keys())

    def mock_cazy_dict(*args, **kwargs):
        return cazy_dict, std_classes

    monkeypatch.setattr(file_io, "get_cazy_dict_std_names", mock_cazy_dict)

    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args_config_cmd["args"],
        null_logger,
    )


def test_parse_config_no_config(args_config_none, cazy_dictionary, null_logger, monkeypatch):
    """Test parse_configuration when no configuration is given."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    std_classes = list(cazy_dict.keys())

    def mock_cazy_dict(*args, **kwargs):
        return cazy_dict, std_classes

    monkeypatch.setattr(file_io, "get_cazy_dict_std_names", mock_cazy_dict)

    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args_config_none["args"],
        null_logger,
    )


# test get_cazy_dict_std_names()


def test_cazy_dict_not_found(null_logger):
    """Test get_cazy_dict_std_names() when the cazy_dictionary file cannot be found."""
    fake_path = parse.__file__

    with pytest.raises(SystemExit) as pytest_wrapped_err:
        file_io.get_cazy_dict_std_names(fake_path, null_logger)
    assert pytest_wrapped_err.type == SystemExit


def test_cazy_dict_success(null_logger):
    """Test get_cazy_dict_std_names() when the retrieval of the data is successful."""
    path = file_io.__file__

    file_io.get_cazy_dict_std_names(path, null_logger)


# test get_yaml_configuration


def test_get_yaml_not_found(empty_config_dict, cazy_dictionary, args_no_yaml, null_logger):
    """Test get_yaml_configuration when the configuration cannot be found."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    std_classes = list(cazy_dict.keys())

    with pytest.raises(SystemExit) as pytest_wrapped_err:
        file_io.get_yaml_configuration(
            empty_config_dict,
            cazy_dict,
            std_classes,
            args_no_yaml["args"],
            null_logger,
        )
    assert pytest_wrapped_err.type == SystemExit


def test_get_yaml_success(empty_config_dict, cazy_dictionary, args_config_file, null_logger):
    """Test get_yaml_configuration when successfully retrieves and parses data from yaml file."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    std_classes = list(cazy_dict.keys())

    file_io.get_yaml_configuration(
        empty_config_dict,
        cazy_dict,
        std_classes,
        args_config_file["args"],
        null_logger,
    )


# test get_yaml_cazy_classes


def test_get_yaml_classes_no(cazy_dictionary, null_logger):
    """Test get_yaml_cazy_classes when classes key cannot be found."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    std_classes = list(cazy_dict.keys())
    yaml_dict = {"Glycoside Hydrolases (GHs)":["GH1","GH2"]}

    file_io.get_yaml_cazy_classes(yaml_dict, cazy_dict, std_classes, null_logger)


def test_get_yaml_classes_none(cazy_dictionary, null_logger):
    """Test get_yaml_cazy_classes when classes key has the value of None."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    std_classes = list(cazy_dict.keys())
    yaml_dict = {"classes": None, "Glycoside Hydrolases (GHs)":["GH1","GH2"]}

    file_io.get_yaml_cazy_classes(yaml_dict, cazy_dict, std_classes, null_logger)


def test_get_yaml_classes_success(cazy_dictionary, null_logger):
    """Test get_yaml_cazy_classes when classes key has the value of CAZy classes."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    std_classes = list(cazy_dict.keys())
    yaml_dict = {"classes": ["PL","CE"], "Glycoside Hydrolases (GHs)":["GH1","GH2"]}

    file_io.get_yaml_cazy_classes(yaml_dict, cazy_dict, std_classes, null_logger)


# test parse_user_cazy_classes()


def test_parse_user_cazy_classes(cazy_dictionary, null_logger):
    """Test parse_user_cazy_classes to standardise the CAZy classes written by the user."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

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

    assert exepected == file_io.parse_user_cazy_classes(
        cazy_classes,
        cazy_dict,
        class_name,
        null_logger,
    )


def test_cannot_standardise_user_classes(cazy_dictionary, null_logger):
    """Test when cannot standardise a name listed in the configuration file."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

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

    assert exepected == file_io.parse_user_cazy_classes(
        cazy_classes,
        cazy_dict,
        class_name,
        null_logger,
    )


# test get_cmd_defined_fams_classes()


def get_cmd_configuration(args_config_cmd, cazy_dictionary, null_logger):
    """Test get_cmd_defined_fams_classes."""
    with open(cazy_dictionary, "r") as fh:
        cazy_dict = json.load(fh)

    std_classes = list(cazy_dict.keys())
    file_io.get_cmd_defined_fams_classes(
        cazy_dict,
        std_classes,
        args_config_cmd["args"],
        null_logger,
    )


# test write_out_df()


def test_writing_named_df_f_true(testing_df, null_logger, making_output_dir):
    """Tests function for writing out a prenamed dataframe"""
    file_io.write_out_df(
        testing_df, "test_writing_df.csv", making_output_dir, null_logger, True
    )


def test_writing_named_df_f_false(testing_df, null_logger, making_output_dir):
    """Tests function for writing out a prenamed dataframe"""
    file_io.write_out_df(
        testing_df, "test_writing_df.csv", making_output_dir, null_logger, False
    )


def test_writing_df_exists(testing_df, null_logger, making_output_dir):
    """Tests function for writing out a prenamed dataframe"""
    file_io.write_out_df(
        testing_df, "test_writing_existing_df", making_output_dir, null_logger, True
    )


# test write_out_failed_scrapes


def test_write_out_failed_scrapes_stdout(stdout_args, null_logger):
    """Test write_out_failed_scrapes when args.output is sys.STDOUT"""
    file_io.write_out_failed_scrapes(
        ["url_1 - message", "url_2 - fail reason"],
        "time_stamp",
        stdout_args["args"],
        null_logger,
    )


def test_write_out_failed_scrapes(output_args, null_logger):
    """Test write_out_failed_scrapes when args.output is sys.STDOUT"""
    file_io.write_out_failed_scrapes(
        ["url_1 - message", "url_2 - fail reason"],
        "time_stamp",
        output_args["args"],
        null_logger,
    )
