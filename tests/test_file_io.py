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
            config=None
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
def config_no_class(input_dir):
    path = input_dir / "config_no_classes.yaml"
    return path


@pytest.fixture
def args_config_no_class(config_no_class):
    args_dict = {
        "args": Namespace(
            config=config_no_class,
        )
    }
    return args_dict


@pytest.fixture
def config_with_classes(input_dir):
    path = input_dir / "config_with_classes.yaml"
    return path


@pytest.fixture
def args_config_with_classes(config_with_classes):
    args_dict = {
        "args": Namespace(
            config=config_with_classes,
        )
    }
    return args_dict


@pytest.fixture
def config_no_class_tag(input_dir):
    path = input_dir / "config_no_class_tag.yaml"
    return path


@pytest.fixture
def args_config_no_class_tag(config_no_class_tag):
    args_dict = {
        "args": Namespace(
            config=config_no_class_tag,
        )
    }
    return args_dict


@pytest.fixture
def config_all_classes(input_dir):
    path = input_dir / "config_all_classes.yaml"
    return path


@pytest.fixture
def args_config_all_classes(config_all_classes):
    args_dict = {
        "args": Namespace(
            config=config_all_classes,
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


def test_parse_config_sys_exit(args_config_None, null_logger):
    """Tests parse_configuration when cazy_dict cannot be opened and should perform sys.exit(1)."""
    fake_path = parse.__file__

    with pytest.raises(SystemExit) as pytest_wrapped_err:
        file_io.parse_configuration(fake_path, args_config_None["args"], null_logger)
    assert pytest_wrapped_err.type == SystemExit


def test_parse_config_none(cazy_dictionary, file_io_path, args_config_None, null_logger):
    """Test parse_configuration when no configuration file provided."""
    with open(cazy_dictionary, "r") as fh:
        expected_cazy_dict = json.load(fh)

    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args_config_None["args"],
        null_logger,
    )

    assert excluded_classes is None
    assert config_dict is None
    assert cazy_dict == expected_cazy_dict


def test_parse_config_cant_find(cazy_dictionary, file_io_path, args_config_fake, null_logger):
    """Test parse_configuration when the configuration file cannot be found."""
    with open(cazy_dictionary, "r") as fh:
        expected_cazy_dict = json.load(fh)

    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args_config_fake["args"],
        null_logger,
    )

    assert excluded_classes is None
    assert config_dict is None
    assert cazy_dict == expected_cazy_dict


def test_parse_config_no_classes(
    cazy_dictionary,
    file_io_path,
    args_config_no_class,
    null_logger,
):
    """Test parse_configuration when no classes are included."""
    with open(cazy_dictionary, "r") as fh:
        expected_cazy_dict = json.load(fh)

    expected_excluded_classes = [
        '<strong>Auxiliary Activities (AAs)</strong>',
        '<strong>Carbohydrate Esterases (CEs)</strong>',
        '<strong>Polysaccharide Lyases (PLs)</strong>',
        '<strong>Carbohydrate-Binding Modules (CBMs)</strong>',
        '<strong>GlycosylTransferases (GTs)</strong>',
    ]

    expected_config_dict = {
        'classes': None,
        'Glycoside Hydrolases (GHs)': ['GH147'],
        'GlycosylTransferases (GTs)': None,
        'Polysaccharide Lyases (PLs)': None,
        'Carbohydrate Esterases (CEs)': None,
        'Auxiliary Activities (AAs)': None,
        'Carbohydrate-Binding Modules (CBMs)': None,
    }

    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args_config_no_class["args"],
        null_logger,
    )

    result = None

    for item in excluded_classes:
        if item not in expected_excluded_classes:
            result = "fail"
    for item in expected_excluded_classes:
        if item not in excluded_classes:
            result = "fail"

    assert result is None
    assert config_dict == expected_config_dict
    assert cazy_dict == expected_cazy_dict


def test_parse_config_with_classes(
    cazy_dictionary,
    file_io_path,
    args_config_with_classes,
    null_logger,
):
    """Test parse_configuration when no classes are included."""
    with open(cazy_dictionary, "r") as fh:
        expected_cazy_dict = json.load(fh)

    expected_excluded_classes = [
        '<strong>Auxiliary Activities (AAs)</strong>',
        '<strong>Carbohydrate-Binding Modules (CBMs)</strong>',
        '<strong>GlycosylTransferases (GTs)</strong>',
        '<strong>Carbohydrate Esterases (CEs)</strong>',
    ]

    expected_config_dict = {
        'classes': ['Glycoside Hydrolases (GHs)', 'Polysaccharide Lyases (PLs)'],
        'Glycoside Hydrolases (GHs)': ['GH147'],
        'GlycosylTransferases (GTs)': None,
        'Polysaccharide Lyases (PLs)': None,
        'Carbohydrate Esterases (CEs)': None,
        'Auxiliary Activities (AAs)': None,
        'Carbohydrate-Binding Modules (CBMs)': None,
    }

    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args_config_with_classes["args"],
        null_logger,
    )

    result = None

    for item in excluded_classes:
        if item not in expected_excluded_classes:
            result = "fail"
    for item in expected_excluded_classes:
        if item not in excluded_classes:
            result = "fail"

    assert result is None
    assert config_dict == expected_config_dict
    assert cazy_dict == expected_cazy_dict


def test_parse_config_no_class_tag(
    cazy_dictionary,
    file_io_path,
    args_config_no_class_tag,
    null_logger,
):
    """Test parse_configuration when no classes are included."""
    with open(cazy_dictionary, "r") as fh:
        expected_cazy_dict = json.load(fh)

    expected_excluded_classes = [
        '<strong>Auxiliary Activities (AAs)</strong>',
        '<strong>Carbohydrate Esterases (CEs)</strong>',
        '<strong>Polysaccharide Lyases (PLs)</strong>',
        '<strong>Carbohydrate-Binding Modules (CBMs)</strong>',
        '<strong>GlycosylTransferases (GTs)</strong>',
    ]

    expected_config_dict = {
        'Glycoside Hydrolases (GHs)': ['GH147'],
        'GlycosylTransferases (GTs)': None,
        'Polysaccharide Lyases (PLs)': None,
        'Carbohydrate Esterases (CEs)': None,
        'Auxiliary Activities (AAs)': None,
        'Carbohydrate-Binding Modules (CBMs)': None,
    }

    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args_config_no_class_tag["args"],
        null_logger,
    )

    result = None

    for item in excluded_classes:
        if item not in expected_excluded_classes:
            result = "fail"
    for item in expected_excluded_classes:
        if item not in excluded_classes:
            result = "fail"

    assert result is None
    assert config_dict == expected_config_dict
    assert cazy_dict == expected_cazy_dict


def test_config_all_classes(
    cazy_dictionary,
    file_io_path,
    args_config_all_classes,
    null_logger,
):
    """Test parse_configuration when there are no excluded classes."""
    with open(cazy_dictionary, "r") as fh:
        expected_cazy_dict = json.load(fh)

    expected_config_dict = {
        'classes': [
            'Glycoside Hydrolases (GHs)',
            'Polysaccharide Lyases (PLs)',
            'GlycosylTransferases (GTs)',
            'Carbohydrate Esterases (CEs)',
            'Auxiliary Activities (AAs)',
            'Carbohydrate-Binding Modules (CBMs)'
        ],
        'Glycoside Hydrolases (GHs)': ['GH147'],
        'GlycosylTransferases (GTs)': ['GT2'],
        'Polysaccharide Lyases (PLs)': ['PL1'],
        'Carbohydrate Esterases (CEs)': ['CE1'],
        'Auxiliary Activities (AAs)': ['AA10'],
        'Carbohydrate-Binding Modules (CBMs)': ['CBM5'],
    }

    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args_config_all_classes["args"],
        null_logger,
    )

    assert excluded_classes is None
    assert config_dict == expected_config_dict
    assert cazy_dict == expected_cazy_dict


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


def test_cannot_standardise(cazy_dictionary, null_logger):
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
