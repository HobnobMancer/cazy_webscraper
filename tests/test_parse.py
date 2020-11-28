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

"""Tests the module parse which parses the protein data to produce a dataframe.

These test are intened to be run from the root of the repository using:
pytest -v

"""

import pytest

import numpy as np
import pandas as pd

from argparse import Namespace

from Bio.PDB import PDBList

from scraper import cazy_webscraper, parse


@pytest.fixture
def out_dir(test_dir):
    path = test_dir / "test_outputs" / "test_outputs_file_io"
    return path


@pytest.fixture
def family():
    family = cazy_webscraper.Family("GH1_test", "GHtest")

    protein = cazy_webscraper.Protein(
        "protein_name",
        "GH1",
        "1.2.3.4",
        "organism",
        {"GenBank": ["link1"], "UniProt": ["link2"], "PDB/3D": ["link3"]},
        {"link1": ["A", "B", "C"]},
    )
    family.members.add(protein)

    protein = None
    family.members.add(protein)

    protein = cazy_webscraper.Protein(
        "protein",
        "GH1",
        "",
        "organism",
        {"GenBank": ["link1"], "UniProt": ["link2"], "PDB": ["link3"]},
    )
    family.members.add(protein)

    return family


@pytest.fixture
def args_ds_fam_no_subfam(test_dir):
    args_dict = {
        "args": Namespace(
            subfamilies=False,
            data_split="family",
            output=test_dir,
            force=False,
            genbank="dummy_email",
            genbank_output=test_dir,
            pdb="pdb",
            pdb_output=test_dir,
        )
    }
    return args_dict


@pytest.fixture
def args_ds_class_subfam(test_dir):
    args_dict = {
        "args": Namespace(
            subfamilies=True,
            data_split="class",
            output=test_dir,
            force=False,
            genbank="dummy_email",
            genbank_output=test_dir,
            pdb=None,
            pdb_output=None,
        )
    }
    return args_dict


@pytest.fixture
def args_ds_none(test_dir):
    args_dict = {
        "args": Namespace(
            subfamilies=False,
            data_split=None,
            output=test_dir,
            force=False,
            genbank=None,
            genbank_output=None,
            pdb="pdb",
            pdb_output=test_dir,
        )
    }
    return args_dict


@pytest.fixture
def protein_df():
    column_names = [
        "Protein_name",
        "CAZy_family",
        "EC#",
        "Source_organism",
        "GenBank",
        "UniProt",
        "PDB/3D",
    ]
    
    df_data = [
        ["protein_1", "GH1", "1.2.3.4", "bact", "GB1", "U1,\nU2", "P1[A]"],
        ["protein_2", "PL2", "1.2.3.4\n2.4.5.6", "euk", "GB2", "U1", np.nan],
        ["protein_3", "CE3", np.nan, "bact", "GB2", np.nan, "P1,\nP2[B]"],
    ]

    df = pd.DataFrame(df_data, columns=column_names)
    
    return df


# test proteins_to_dataframe() (dataframe building function)


def test_prt_to_df_ds_fam_no_subfams(args_ds_fam_no_subfam, family, null_logger, monkeypatch):
    """Test proteins_to_dataframe when data split is family and subfamilies is False."""

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "write_out_df", mock_no_return)
    monkeypatch.setattr(parse, "get_structures_and_sequences", mock_no_return)

    parse.proteins_to_dataframe([family], args_ds_fam_no_subfam["args"], null_logger)


def test_prt_to_df_ds_class_subfams(args_ds_class_subfam, family, null_logger, monkeypatch):
    """Test proteins_to_dataframe when data split is class and subfamilies is True."""

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "write_out_df", mock_no_return)
    monkeypatch.setattr(parse, "get_genbank_fasta", mock_no_return)

    parse.proteins_to_dataframe([family], args_ds_class_subfam["args"], null_logger)


def test_prt_to_df_ds_none(args_ds_none, family, null_logger, monkeypatch):
    """Test proteins_to_dataframe when data split is None."""

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "write_out_df", mock_no_return)
    monkeypatch.setattr(parse, "get_pdb_structures", mock_no_return)

    parse.proteins_to_dataframe([family], args_ds_none["args"], null_logger)


# test get_structures_and_sequences()


def test_gt_s_p_pdb_out_given():
    # protein_df
    # df name
    # args
    # logger
def test_gt_s_p_pdb_to_output():
def test_gt_s_p_pdb_to_cwd():