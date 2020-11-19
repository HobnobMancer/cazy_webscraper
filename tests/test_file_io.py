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

import pandas as pd

from scraper import file_io


@pytest.fixture
def test_output_dir(test_dir):
    path = test_dir / "test_outputs" / "test_outputs_file_io"
    return path


@pytest.fixture
def making_output_dir(test_output_dir):
    path = test_output_dir / "testing_making_dir"
    return path


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


# test make_output_directory()


def test_output_dir_creation_nd_true(making_output_dir, null_logger):
    """Test creation of output dir when nodelete is false"""
    file_io.make_output_directory(making_output_dir, null_logger, True, False)


def test_output_dir_creation_nd_false(making_output_dir, null_logger):
    """Test creation of output dir when nodelete is true"""

    file_io.make_output_directory(making_output_dir, null_logger, True, True)


# test parse_configuration()


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
    
# only one test and assert get the correct result back

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
