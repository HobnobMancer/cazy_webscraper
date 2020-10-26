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
"""Module for parsing the webpages from CAZy."""

import pandas as pd

def parse_protein_table(url):
    """Parses protein table from CAZy family page.

    Writes out dataframe to .csv file.

    :param url: str, url to web page containing protein table

    Return parsed dataframe.
    """
    # scrape the protien table page and retrieve all dataframes
    list_of_dfs = pd.read_html(url, header=2, index_col=0, encoding="utf-8")

    # 2 dataframes are returned
    # The first is the summary table at the top of page when viewed in a browser
    # The second df is the protein table of interest
    protein_df= list_of_dfs[1]

    # final row contains links to the top of the webpage
    # remove final row
    protein_df = protein_df.drop((len(protein_df["Organism"]) - 1))

    # column 'Unnamed: 6' only contains content at the end of the dataframe when
    # changing kingdom of documentary species, e.g. the start of the Bacteria species
    # idenfity if this 'kingdom' row is present in the scapped table
    row_index = 0
    for row_index in range(len(protein_df["Organism"])):
        if pd.isnull(protein_df["Unnamed: 6"].iloc[row_index]) is False:
            # 'Unnamed: 6' cell contains data
            # remove this row nad row below containing 
            protein_df.drop(row_index)
            protein_df.drop((row_index + 1))
    
    # remove column 'unnamed: 6'
    protein_df.drop("Unnamed: 6", axis=1)

    return protein_df
