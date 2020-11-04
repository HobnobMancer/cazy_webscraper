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

from scraper.file_io import write_out_df


def proteins_to_dataframe(families, args, logger):
    """Write protein members of CAZy families to a pandas dataframe.

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
        for protein in proteins:
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

    # write out dataframe
    write_out_df(protein_dataframe, df_name, args.output, logger, args.force)

    return
