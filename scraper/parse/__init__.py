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

import pandas as pd

from datetime import datetime

from tqdm import tqdm

from scraper.file_io import write_out_df


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
            "EC#",
            "Source_organism",
            "GenBank",
            "UniProt",
            "PDB/3D"
        ]
    )

    # add proteins to dataframe
    for family in families:
        proteins = family.get_proteins()
        for protein in tqdm(proteins, desc=f"Writing {family.name} proteins to df"):
            if protein is not None:
                protein_dict = protein.get_protein_dict()
                df = pd.DataFrame(protein_dict)
                protein_dataframe = protein_dataframe.append(df)

    time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

    # compile dataframe name
    if args.data_split == "family":
        df_name = f"cazy_{families[0].name}_{time_stamp}.csv"
    elif args.data_split == "class":
        df_name = f"cazy_{families[0].cazy_class}_{time_stamp}.csv"
    else:
        df_name = f"cazy_scrape_{time_stamp}.csv"

    # Remove duplicate
    # This can arise when scraping subfamilies and families within a class
    # The proteins within the subfamilies will also be listed under their parent family
    protein_dataframe = protein_dataframe.drop_duplicates()

    # write out dataframe
    write_out_df(protein_dataframe, df_name, args.output, logger, args.force)

    return
