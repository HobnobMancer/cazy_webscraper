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


import json
import pytest
import sys

import numpy as np
import pandas as pd

from argparse import Namespace
from unittest.mock import patch, mock_open

from Bio.PDB import PDBList

from scraper import cazy_webscraper, parse


@pytest.fixture
def out_dir(test_dir):
    path = test_dir / "test_outputs" / "test_outputs_parse"
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
def args_ds_fam_no_subfam(out_dir):
    args_dict = {
        "args": Namespace(
            subfamilies=False,
            data_split="family",
            output=out_dir,
            force=False,
            genbank="dummy_email",
            genbank_output=out_dir,
            pdb="pdb",
            pdb_output=out_dir,
        )
    }
    return args_dict


@pytest.fixture
def args_ds_class_subfam(test_dir, out_dir):
    args_dict = {
        "args": Namespace(
            subfamilies=True,
            data_split="class",
            output=out_dir,
            force=False,
            genbank="dummy_email",
            genbank_output=test_dir,
            pdb=None,
            pdb_output=None,
        )
    }
    return args_dict


@pytest.fixture
def args_ds_none(test_dir, out_dir):
    args_dict = {
        "args": Namespace(
            subfamilies=False,
            data_split=None,
            output=out_dir,
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
        ["protein_1", "GH1", "1.2.3.4", "bact", "ABC12345.1", "U1", "P1,\nP2[B]"],
        ["protein_2", "PL2", "1.2.3.4\n2.4.5.6", "euk", "1235G", "U1", ""],
        ["protein_3", "CE3", np.nan, "bact", np.nan, np.nan, np.nan],
    ]

    df = pd.DataFrame(df_data, columns=column_names)

    return df


@pytest.fixture
def df_name():
    name = "test_df_name"
    return name


@pytest.fixture
def args_output_no_pdb_output(out_dir):
    args_dict = {
        "args": Namespace(
            subfamilies=False,
            data_split=None,
            output=out_dir,
            force=False,
            genbank="dummy_email",
            genbank_output=out_dir,
            pdb="pdb",
            pdb_output=None,
        )
    }
    return args_dict


@pytest.fixture
def args_no_output_no_pdb_output(out_dir):
    args_dict = {
        "args": Namespace(
            subfamilies=False,
            data_split=None,
            output=sys.stdout,
            force=False,
            genbank="dummy_email",
            genbank_output=out_dir,
            pdb="pdb",
            pdb_output=None,
        )
    }
    return args_dict


@pytest.fixture
def df_row_float(protein_df):
    row = protein_df.iloc[2]
    return row


@pytest.fixture
def df_row_gb_format_wrong(protein_df):
    row = protein_df.iloc[1]
    return row


@pytest.fixture
def df_row_success(protein_df):
    row = protein_df.iloc[0]
    return row


@pytest.fixture
def args_args_gen_out_stdout(test_dir, out_dir):
    args_dict = {
        "args": Namespace(
            subfamilies=False,
            data_split="family",
            output=test_dir,
            force=False,
            genbank="dummy_email",
            genbank_output=sys.stdout,
            pdb="pdb",
            pdb_output=out_dir,
        )
    }
    return args_dict


@pytest.fixture
def sequence_fasta(test_dir):
    path = test_dir / "test_inputs" / "test_inputs_parse" / "sequence.fasta"
    return path


# test proteins_to_dataframe() (dataframe building function)


def test_prt_to_df_ds_fam_no_subfams(args_ds_fam_no_subfam, family, null_logger, monkeypatch):
    """Test proteins_to_dataframe when data split is family and subfamilies is False."""

    with patch("builtins.open", mock_open(read_data="data")) as mock_file:

        def mock_no_return(*args, **kwargs):
            return

        monkeypatch.setattr(parse, "write_out_df", mock_no_return)
        monkeypatch.setattr(parse, "get_structures_and_sequences", mock_no_return)
        monkeypatch.setattr(json, "dump", mock_no_return)

        parse.proteins_to_dataframe([family], args_ds_fam_no_subfam["args"], null_logger)


def test_prt_to_df_ds_class_subfams(args_ds_class_subfam, family, null_logger, monkeypatch):
    """Test proteins_to_dataframe when data split is class and subfamilies is True."""

    with patch("builtins.open", mock_open(read_data="data")) as mock_file:

        def mock_no_return(*args, **kwargs):
            return

        monkeypatch.setattr(parse, "write_out_df", mock_no_return)
        monkeypatch.setattr(parse, "get_genbank_fasta", mock_no_return)

        parse.proteins_to_dataframe([family], args_ds_class_subfam["args"], null_logger)


def test_prt_to_df_ds_none(args_ds_none, family, null_logger, monkeypatch):
    """Test proteins_to_dataframe when data split is None."""

    with patch("builtins.open", mock_open(read_data="data")) as mock_file:

        def mock_no_return(*args, **kwargs):
            return

        monkeypatch.setattr(parse, "write_out_df", mock_no_return)
        monkeypatch.setattr(parse, "get_pdb_structures", mock_no_return)

        parse.proteins_to_dataframe([family], args_ds_none["args"], null_logger)


# test get_structures_and_sequences()


def test_gt_s_p_pdb_out_given(protein_df, df_name, args_ds_fam_no_subfam, null_logger, monkeypatch):
    """Test get_structures_and_sequences when no pdb_output given."""

    def mock_genbank_accessions(*args, **kwargs):
        return("genbank_accession", "GH1")

    def mock_pdb_accessions(*args, **kwargs):
        return["pdb1", "pdb2", "pdb3"]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "get_genbank_accession", mock_genbank_accessions)
    monkeypatch.setattr(parse, "get_pdb_accessions", mock_pdb_accessions)
    monkeypatch.setattr(parse, "download_fasta", mock_no_return)
    monkeypatch.setattr(PDBList, "retrieve_pdb_file", mock_no_return)

    parse.get_structures_and_sequences(
        protein_df,
        df_name,
        args_ds_fam_no_subfam["args"],
        null_logger,
    )


def test_gt_s_p_pdb_to_output(
    protein_df,
    df_name,
    args_output_no_pdb_output,
    null_logger,
    monkeypatch,
):
    """Test get_structures_and_sequences when output but not pdb_output given."""

    def mock_genbank_accessions(*args, **kwargs):
        return("genbank_accession", "GH1")

    def mock_pdb_accessions(*args, **kwargs):
        return["pdb1", "pdb2", "pdb3"]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "get_genbank_accession", mock_genbank_accessions)
    monkeypatch.setattr(parse, "get_pdb_accessions", mock_pdb_accessions)
    monkeypatch.setattr(parse, "download_fasta", mock_no_return)
    monkeypatch.setattr(PDBList, "retrieve_pdb_file", mock_no_return)

    parse.get_structures_and_sequences(
        protein_df,
        df_name,
        args_output_no_pdb_output["args"],
        null_logger,
    )


def test_gt_s_p_pdb_to_cwd(
    protein_df,
    df_name,
    args_no_output_no_pdb_output,
    null_logger,
    monkeypatch,
):
    """Test get_structures_and_sequences when no pdb_output given and output is stdout."""

    def mock_genbank_accessions(*args, **kwargs):
        return("genbank_accession", "GH1")

    def mock_pdb_accessions(*args, **kwargs):
        return["pdb1", "pdb2", "pdb3"]

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "get_genbank_accession", mock_genbank_accessions)
    monkeypatch.setattr(parse, "get_pdb_accessions", mock_pdb_accessions)
    monkeypatch.setattr(parse, "download_fasta", mock_no_return)
    monkeypatch.setattr(PDBList, "retrieve_pdb_file", mock_no_return)

    parse.get_structures_and_sequences(
        protein_df,
        df_name,
        args_no_output_no_pdb_output["args"],
        null_logger,
    )


# test get_genbank_accession()


def test_gb_acc_float(df_row_float, df_name, null_logger):
    """Test get_genbank_accession when the GenBank cell contains a float data type."""
    row_index = 1

    return_acc, returned_fam = parse.get_genbank_accession(
        df_row_float,
        df_name,
        row_index,
        null_logger,
    )

    expected_acc = None
    expected_fam = None

    assert expected_acc == return_acc
    assert expected_fam == returned_fam


def test_gb_acc_re_fail(df_row_gb_format_wrong, df_name, null_logger):
    """Test get_genbank_accession when the retrieve item does not match expected GenBank format."""
    row_index = 1

    return_acc, returned_fam = parse.get_genbank_accession(
        df_row_gb_format_wrong,
        df_name,
        row_index,
        null_logger,
    )

    expected_acc = None
    expected_fam = None

    assert expected_acc == return_acc
    assert expected_fam == returned_fam


def test_gb_acc_successful(df_row_success, df_name, null_logger):
    """Test get_genbank_accession when fully successful."""
    row_index = 1

    return_acc, returned_fam = parse.get_genbank_accession(
        df_row_success,
        df_name,
        row_index,
        null_logger,
    )

    expected_acc = "ABC12345.1"
    expected_fam = "GH1"

    assert expected_acc == return_acc
    assert expected_fam == returned_fam


# test get_pdb_accession()


def test_pdb_acc_float(df_row_float, df_name, null_logger):
    """Test get_pdb_accessions when value in PDB cell is a float."""
    assert None is parse.get_pdb_accessions(df_row_float, df_name, null_logger)


def test_pdb_acc_suss(df_row_success, df_name, null_logger):
    """Test get_pdb_accessions when successfull."""
    expected = ["P1", "P2"]
    assert expected == parse.get_pdb_accessions(df_row_success, df_name, null_logger)


# test get_pdb_structures()


def test_pdb_struc_pdb_out_given(
    protein_df,
    df_name,
    args_ds_fam_no_subfam,
    null_logger,
    monkeypatch,
):
    """Test get_pdb_structures when pdb_output is given."""

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(PDBList, "retrieve_pdb_file", mock_no_return)

    parse.get_pdb_structures(protein_df, df_name, args_ds_fam_no_subfam["args"], null_logger)


def test_pdb_struc_pdb_to_output(
    protein_df,
    df_name,
    args_output_no_pdb_output,
    null_logger,
    monkeypatch,
):
    """Test get_pdb_structures when output not pdb_output is given."""

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(PDBList, "retrieve_pdb_file", mock_no_return)

    parse.get_pdb_structures(protein_df, df_name, args_output_no_pdb_output["args"], null_logger)


def test_pdb_struc_pdb_to_cwd(
    protein_df,
    df_name,
    args_no_output_no_pdb_output,
    null_logger,
    monkeypatch,
):
    """Test get_pdb_structures when output is sys.stdout and no pdb_output is given."""

    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(PDBList, "retrieve_pdb_file", mock_no_return)

    parse.get_pdb_structures(protein_df, df_name, args_no_output_no_pdb_output["args"], null_logger)


# Test get_genbank_fasta()

def test_get_genbank_fasta(protein_df, df_name, args_ds_fam_no_subfam, null_logger, monkeypatch):
    """Test get_genbank_fasta."""

    def mock_download_fasta(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "download_fasta", mock_download_fasta)

    parse.get_genbank_fasta(protein_df, df_name, args_ds_fam_no_subfam["args"], null_logger)


# Test download_fasta()


def test_download_fasta_file_exists(args_ds_fam_no_subfam, out_dir, null_logger, monkeypatch):
    """Test download_fasta() when file exists."""
    file_name = out_dir / "existing_file.txt"
    accession = "test_accession"

    def mock_efetch(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "entrez_retry", mock_efetch)

    parse.download_fasta(accession, file_name, args_ds_fam_no_subfam["args"], null_logger)


def test_download_fasta_handle_none(test_dir, args_ds_fam_no_subfam, null_logger, monkeypatch):
    """Test download_fasta() when returned handle is None."""
    file_name = test_dir / "test_outputs" / "test_output_parse" / "novel_fasta.fasta"
    accession = "test_accession"

    def mock_efetch(*args, **kwargs):
        return

    monkeypatch.setattr(parse, "entrez_retry", mock_efetch)

    parse.download_fasta(accession, file_name, args_ds_fam_no_subfam["args"], null_logger)


def test_download_fasta_when_out_dir_given(
    sequence_fasta,
    out_dir,
    args_ds_fam_no_subfam,
    null_logger,
    monkeypatch,
):
    file_name = out_dir / "novel_fasta.fasta"
    accession = "test_accession"

    def mock_efetch(*args, **kwargs):
        return sequence_fasta

    monkeypatch.setattr(parse, "entrez_retry", mock_efetch)

    parse.download_fasta(accession, file_name, args_ds_fam_no_subfam["args"], null_logger)


def test_download_fasta_when_out_stdout(
    sequence_fasta,
    out_dir,
    args_args_gen_out_stdout,
    null_logger,
    monkeypatch,
):
    file_name = out_dir / "novel_fasta.fasta"
    accession = "test_accession"

    def mock_efetch(*args, **kwargs):
        return sequence_fasta

    monkeypatch.setattr(parse, "entrez_retry", mock_efetch)

    parse.download_fasta(accession, file_name, args_args_gen_out_stdout["args"], null_logger)


# test entrez_retry()


def test_entry_retry(null_logger):
    """Test entrez_retry."""

    def mock_record(*args, **kwargs):
        return "test_record"

    assert "test_record" == parse.entrez_retry(null_logger, mock_record)


# def test_entrez_retry_none(null_logger):
#     """Test entrez_retry when nothing is returned."""

#     def mock_record(*args, **kwargs):
#         return

#     assert parse.entrez_retry(null_logger, mock_record) is None
