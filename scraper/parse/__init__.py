#!/usr/bin/env python3
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
"""Module for parsing the webpages from CAZy."""


import re
import sys

import pandas as pd

from datetime import datetime

from Bio import Entrez
from tqdm import tqdm

from scraper.file_io import make_output_directory, write_out_df
from scraper.genbank import download_fasta, get_genbank_fasta
from scraper.pdb import download_pdb


def proteins_to_dataframe(families, args, logger):
    """Write protein members of CAZy families to a pandas dataframe.

    Duplicate rows are removed. These rows must be identical across the entire row, not only share
    the same protein name. Duplicates are mostly likely to arise when scraping families and
    subfamilies becuase proteins under the subfamily are also catalogued under the parent family.

    :param families: list, list of CAZy families containing protein members
    :param args: args parser
    :param logger: logger object

    Return nothing.
    """
    # Build empty dataframe to add proteins to
    protein_dataframe = pd.DataFrame(
        data={},
        columns=[
            "Protein_name",
            "CAZy_family",
            "EC#",
            "Source_organism",
            "GenBank",
            "UniProt",
            "PDB/3D",
        ]
    )

    # add proteins to dataframe
    for family in families:
        proteins = family.get_proteins()
        for protein in tqdm(proteins, desc=f"Writing {family.name} proteins to df"):
            if protein is not None:
                protein_dict = protein.get_protein_dict()
                df = pd.DataFrame(protein_dict)
                protein_dataframe = protein_dataframe.append(df, ignore_index=True)

    time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

    if args.subfamilies is True:
        subfam = "_incld_subfams"
    else:
        subfam = ""
    # compile dataframe name
    if args.data_split == "family":
        df_name = f"cazy_{families[0].name}{subfam}_{time_stamp}"
    elif args.data_split == "class":
        df_name = f"cazy_{families[0].cazy_class}{subfam}_{time_stamp}"
    else:
        df_name = f"cazy_scrape{subfam}_{time_stamp}"

    # Remove duplicate
    # This can arise when scraping subfamilies and families within a class
    # The proteins within the subfamilies will also be listed under their parent family
    protein_dataframe = protein_dataframe.drop_duplicates()

    # Additional parsing and retrieval of data: protein sequences and/or structures
    if args.genbank and args.pdb:
        get_structures_and_sequences(protein_dataframe, df_name, args, logger)

    elif (args.genbank is not None) and (args.pdb is None):
        get_genbank_fasta(protein_dataframe, df_name, args, logger)

    elif (args.genbank is None) and (args.pdb is not None):
        get_pdb_structures(protein_dataframe, df_name, args, logger)

    # write out dataframe
    write_out_df(protein_dataframe, df_name, args.output, logger, args.force)
    return


def get_structures_and_sequences(protein_dataframe, df_name, args, logger):
    """Coordinate the retrieval of CAZyme structure and sequences from PDB and GenBank.

    :param protein_dataframe: Pandas dataframe, dataframe containing protein data from CAZy
    :param df_name: str, name of the CAZy dataframe
    :param args: cmd args parser
    :param logger: logger object


    Return nothing.
    """
    logger.info(f"Retrieving structures and sequences, for proteins in {df_name}")

    # set user email address for Entrez
    Entrez.email = args.genbank

    # create directory to write FASTA files to
    if (args.genbank_output is not sys.stdout) and (args.genbank_output != args.output):
        make_output_directory(args.genbank_output, logger, args.force, args.nodelete)
    # create directory to write strucure files to
    if (args.pdb_output is not sys.stdout) and (args.pdb_output != args.output):
        make_output_directory(args.pdb_output, logger, args.force, args.nodelete)

    index = 0
    for index in tqdm(range(len(protein_dataframe["Protein_name"])), desc="Downloading PDBs and FASTAs"):
        # Retrieve accession from GenBank cell
        df_row = protein_dataframe.iloc[index]
        accession, cazy_family = get_accession(df_row, df_name, index, logger)

        if accession is None:
            continue

        # get FASTA file from GenBank
        fasta_name = f"{accession}_{cazy_family}.fasta"

        if args.genbank_output is not sys.stdout:
            fasta_name = args.genbank_output / fasta_name

        download_fasta(accession, fasta_name, args, logger)

        # Get .pdb from PDB
        pdb_name = f"{accession}_{cazy_family}.pdb"

        download_pdb(accession, pdb_name, args, logger)

    return


def get_accession(df_row, df_name, row_index, logger):
    """Retrieve GenBank accession for protein in the dataframe (df) row.

    Retrieve the first GenBank accession becuase if multiple are given, CAZy only links
    to the first GenBank accession.

    :param df_row: Pandas series, from protein protein dataframe
    :param df_name: str, name of the Pandas dataframe from which row was retrieved
    :param row_index: int, index of the row in the protein dataframe
    :param logger: logger object

    Return two strings, GenBank accession of the protein and the protein's CAZy family.
    """
    family = df_row[1]
    genbank_cell = df_row[4]

    # Retrieve the first accession number, and separate from assoicated HTML link
    find_result = genbank_cell.find(" ")  # finds space separating the first accession and html link
    first_accession = genbank_cell[:find_result]

    # check result is in expected genbank format
    try:
        re.match(r"\D{3}\d+.\d+", first_accession)
    except AttributeError:
        logger.warning(
            f"Could not return accession for protein in row {index} in\n"
            "{df_name}.\n"
            "Not retrieving FASTA file for this protein."
        )
        return None, None

    return first_accession, family
