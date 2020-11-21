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
"""Module for retrieving FASTA files from GenBank."""

import re
import sys
import time

from socket import timeout
from typing import List, Optional
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import pandas as pd

from Bio import Entrez, SeqIO
from tqdm import tqdm

from scraper.file_io import make_output_directory


def get_genbank_fasta(dataframe, df_name, args, logger):
    """Retrieve GenBank accessions, coordinate downloading associated GenBank FASTA file.

    :param:

    Return  nothing.
    """
    logger.info(f"Retrieving FASTA files from GenBank, for proteins in {df_name}")

    # set user email address for Entrez
    Entrez.email = args.genbank

    # create directory to write FASTA files to
    if (args.genbank_output is not sys.stdout) and (args.genbank_output != args.output):
        make_output_directory(args.genbank_output, logger, args.force, args.nodelete)

    index = 0
    for index in tqdm(range(len(dataframe["Protein_name"])), desc="Downloading GenBank FASTAs"):
        # Retrieve accession from GenBank cell
        df_row = dataframe.iloc[index]
        accession, cazy_family = get_accession(df_row, df_name, index, logger)

        if accession is None:
            continue

        # create file name
        file_name = f"{accession}_{cazy_family}.fasta"
        if args.genbank_output is not sys.stdout:
            file_name = args.genbank_output / file_name

        download_fasta(accession, file_name, args, logger)

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


def download_fasta(accession, file_name, args, logger):
    """Download FASTA file from GenBank FTP server.

    :param accession: str, accession of protein
    :param args: cmd args parser
    :param logger: logger object

    Returning nothing.
    """
    if args.genbank_output is not sys.stdout:
        if file_name.exists():
            logger.warning(f"FASTA file {file_name} aleady exists, not downloading again")
            return

    record, retries, tries = None, args.retries, 0

    while record is None and tries < retries:
        try:
            record = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")

        except IOError:
            if tries < retries:
                logger.warning(
                    "Network error encoutered when trying to download FASTA file from GenBank \n"
                    f"for {accession}, during try no. {tries}.\n"
                    "Retrying in 10s"
                )
                time.sleep(10)

            tries += 1

    if record is None:
        logger.error(
            "Network error encountered too many times.\n"
            f"Not retrieving FASTA file for {accession}"
        )
        return

    if args.genbank_output is not sys.stdout:
        with open(file_name, "w") as fh:
            fh.write(record)

    else:
        SeqIO.write(record, sys.stdout, "fasta")

    return
