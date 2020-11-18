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

import pytest

from scraper import file_io


@pytest.fixture
def test_output_dir(test_dir):
    path = test_dir / "test_outputs" / "test_outputs_file_io"
    return path


@pytest.fixture
def making_output_dir(test_output_dir):
    path = test_output_dir / "testing_making_dir"
    return path


# test make_output_directory()


def test_output_dir_creation_nd_true(making_output_dir, null_logger):
    """Test creation of output dir when nodelete is false"""
    file_io.make_output_directory(making_output_dir, null_logger, True, False)


def test_output_dir_creation_nd_false(making_output_dir, null_logger):
    """Test creation of output dir when nodelete is true"""

    file_io.make_output_directory(making_output_dir, null_logger, True, True)


# test parse_configuration()


# test parse_user_cazy_classes()


# test write_out_df()
