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

from argparse import Namespace

from scraper import cazy_webscraper, parse


@pytest.fixture
def family():
    family = cazy_webscraper.Family("GH1_test", "GHtest")

    protein = cazy_webscraper.Protein("protein_name", "GH1", "1.2.3.4", "organism")
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
        )
    }
    return args_dict


def test_parse_ds_fam_no_subfams(args_ds_fam_no_subfam, family, null_logger, monkeypatch):
    """Test proteins_to_dataframe when data split is family and subfamilies is False."""

    def mock_write_df(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "write_out_df", mock_write_df)

    parse.proteins_to_dataframe([family], args_ds_fam_no_subfam["args"], null_logger)



def test_parse_ds_class_subfams(args_ds_class_subfam, family, null_logger, monkeypatch):
    """Test proteins_to_dataframe when data split is class and subfamilies is True."""

    def mock_write_df(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "write_out_df", mock_write_df)

    parse.proteins_to_dataframe([family], args_ds_class_subfam["args"], null_logger)


def test_parse_ds_none(args_ds_none, family, null_logger, monkeypatch):
    """Test proteins_to_dataframe when data split is None."""

    def mock_write_df(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "write_out_df", mock_write_df)

    parse.proteins_to_dataframe([family], args_ds_none["args"], null_logger)
