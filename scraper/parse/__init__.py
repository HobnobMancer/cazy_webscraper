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


def parse_cazy_protein_data(site, cazy_classes, args, logger):
    """Coordinate parsing protein tables from CAZy, and spliting data.

    Data can no be split (option 'all') so a single dataframe containing
    all protein data is created.

    Data can be split by class (option 'class') creating a single dataframe
    per CAZy class, or be split by family (option 'family') creating a single
    dataframe per CAZy family.

    :param site: site class object, contains links to all CAZy pages
    :param cazy_classes: tuple, each list containing the full length and abbreviated class name
    :param args: parser object
    :param logger: logger object

    Return nothing.
    """
    if args.data_split == "all":  # create a single dataframe containing ALL data
        # retrieve the urls to all protein tables pages in site
        all_protein_table_page_urls = 

        # build empty df to store all protein data from CAZy
        all_protein_df = pd.DataFrame(columns=["Protein Name", "EC#", "Organism", "GenBank", "Uniprot", "PDB/3D", "Unnamed: 6"])

        for url in all_protein_table_pages:
            new_protein_df = parse_protein_table(url)
            all_protein_df = all_protein_df.append(new_protein_df, ignore_index=True)
        
        # create time stamp
        time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")
        
        # write out protein dataframe to user specified output directory
        write_out_df(all_protein_df, f"cazy_scrape_{time_stamp}.csv", args.output)
    
    elif args.data_split == "class":  # create a dataframe PER CAZy CLASS
        for class_list in cazy_classes:
            cazy_class = class_list[0]
            class_abbrev = class_list[1]

            # retrieve urls to all protein tables pages for given CAZy class
            class_urls = 
        
            # build empty dataframe to store protein data for CAZy class
            class_protein_df = pd.DataFrame(columns=["Protein Name", "EC#", "Organism", "GenBank", "Uniprot", "PDB/3D", "Unnamed: 6"])
            for url in class_urls:
                new_protein_df = parse_protein_table(url)
                class_protein_df = all_protein_df.append(new_protein_df, ignore_index=True)

            # create time stamp
            time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

            # write out dataframe for given CAZy class to user output directory
            write_out_df(class_protein_df, f"cazy_scrape_{cazy_class}_{time_stamp}.csv", args.output)
        
    elif args.data_split == "family":  # create a dataframe PER CAZy FAMILY
        for class_list in cazy_classes:
            cazy_class = class_list[0]
            class_abbrev = class_list[1]

            # retrieve all Family pages for given class
            families = 

            for family in families:
                # retrieve urls to all protein tables for given CAZy family
                family_urls = 

                # create empty df to store protein data from family
                family_protein_df = pd.DataFrame(columns=["Protein Name", "EC#", "Organism", "GenBank", "Uniprot", "PDB/3D", "Unnamed: 6"])
                
                for url in family_urls:
                    new_protein_df = parse_protein_table(url)
                    family_protein_df = all_protein_df.append(new_protein_df, ignore_index=True)

                # create time stamp
                time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")
                # write out dataframe for CAZy Family
                write_out_df(class_protein_df, f"cazy_scrape_{cazy_family}_{time_stamp}.csv", args.output)
        

def parse_protein_table(url):
    """Parses protein table from CAZy family page.

    Writes out dataframe to .csv file.

    :param url: str, url to web page containing protein table

    Return parsed dataframe.
    """
    # scrape the protien table page and retrieve all dataframes
    list_of_dfs = pd.read_html(url, header=2, encoding="utf-8")

    # 2 dataframes are returned
    # The first is the summary table at the top of page when viewed in a browser
    # The second df is the protein table of interest
    protein_df = list_of_dfs[1]

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
            # remove this row nad row below containing repeated column headings
            protein_df.drop(row_index)
            protein_df.drop((row_index + 1))

    # remove column 'unnamed: 6'
    protein_df.drop("Unnamed: 6", axis=1)

    return protein_df
