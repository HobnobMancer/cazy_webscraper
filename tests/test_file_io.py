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

"""Tests the submodule file_io from utilities.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest
import shutil
import sys

import pandas as pd

from argparse import Namespace

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scraper.utilities import file_io


@pytest.fixture
def test_output_dir(test_dir):
    path = test_dir / "test_outputs" / "test_outputs_file_io"
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
def args_config_cmd():
    args_dict = {
        "args": Namespace(
            config=None,
            classes="CE,AA",
            families="GH14,GT1,PL15,CE6,AA10,CBM50,DD21,CBMAA1",
            genera=None,
            species=None,
            strains=None,
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
    path = test_dir / "test_outputs" / "test_outputs_file_io"
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


# test make_output_directory()


def test_making_new_dir_1(making_output_dir):
    """Test making a new directory."""
    path_ = making_output_dir / "test_build"
    file_io.make_output_directory(path_, True, True)
    shutil.rmtree(path_)


def test_making_new_dir_2(making_output_dir):
    """Test making a new directory."""
    path_ = making_output_dir / "test_build"
    file_io.make_output_directory(path_, False, True)
    shutil.rmtree(path_)


def test_output_dir_creation_nd_true(making_output_dir):
    """Test creation of output dir when nodelete is false"""
    file_io.make_output_directory(making_output_dir, True, False)


def test_output_dir_creation_nd_false(making_output_dir):
    """Test creation of output dir when nodelete is true"""
    file_io.make_output_directory(making_output_dir, True, True)


def test_output_dir_creation_exists(test_dir):
    """Test creation of output dir when already exists."""
    with pytest.raises(SystemExit) as pytest_wrapped_err:
        file_io.make_output_directory(test_dir, False, False)
    assert pytest_wrapped_err.type == SystemExit


# test write_out_failed_scrapes


def test_write_out_failed_scrapes_stdout(stdout_args):
    """Test write_out_failed_scrapes when args.output is sys.STDOUT"""
    file_io.write_out_failed_scrapes(
        ["url_1 - message", "url_2 - fail reason"],
        "time_stamp",
        stdout_args["args"],
    )


def test_write_out_failed_scrapes(output_args):
    """Test write_out_failed_scrapes when args.output is sys.STDOUT"""
    file_io.write_out_failed_scrapes(
        ["url_1 - message", "url_2 - fail reason"],
        "time_stamp",
        output_args["args"],
    )


# test write_out_failed_proteins()


def test_write_out_failed_proteins_stdout(stdout_args):
    """Test write_out_failed_proteins when args.output is sys.STDOUT"""
    file_io.write_out_failed_proteins(
        ["sql_1 - message", "sql_2 - fail reason"],
        "time_stamp",
        stdout_args["args"],
    )


def test_write_out_failed_proteins(output_args):
    """Test write_out_failed_proteins when args.output is sys.STDOUT"""
    file_io.write_out_failed_proteins(
        ["sql_1 - message", "sql_2 - fail reason"],
        "time_stamp",
        output_args["args"],
    )


# test writing out FASTA files


def test_write_out_fasta(making_output_dir):
    """Test writing out FASTA file."""
    path_ = making_output_dir / "fasta_test"
    file_io.make_output_directory(path_, True, True)
    genbank_accession = "test_accession"
    record = SeqRecord(
        Seq("MKQHKAMIVALIVTAVVAALVTRKDLCEHIRTGQTEVAVAVF"),
        id="fake_protein.1",
        name="fake",
        description="test protein record",
    )
    args = {"args": Namespace(write=path_)}
    file_io.write_out_fasta(record, genbank_accession, args["args"])
    shutil.rmtree(path_)
