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

from tqdm import tqdm


def get_genbank_fasta(dataframe, args, logger):
    """Retrieve GenBank accessions, coordinate downloading associated GenBank FASTA file.

    :param:

    Return  nothing.
    """
    index = 0
    for index in tqdm(range(len(dataframe["Protein_name"])), desc="Downloading GenBank FASTAs"):
        # Retrieve accession from GenBank cell
        df_row = dataframe.iloc[index]
        accession = get_accession(df_row, logger)
    
    if accession is None:
        return

    
    

def get_accession(df_row, logger):
    """Retrieve GenBank accession for protein in the dataframe (df) row.

    Retrieve the first GenBank accession becuase if multiple are given, CAZy only links
    to the first GenBank accession.

    :param df_row: Pandas series, from protein protein dataframe
    :param logger: logger object

    Return str, GenBank accession of the protein.
    """
    genbank_cell = dr_row[4]

    # Retrieve the first accession number, and separate from assoicated HTML link
    find_result = genbank_cell.find(" ")  # finds space separating the first accession and html link
    first_accession = genbank_cell[:find_result]

    # check result is in expected genbank format
    try:
        re.match(r"\D{3}\d+.\d+", first_accession)
    except AttributeError:
        return

    return first_accession