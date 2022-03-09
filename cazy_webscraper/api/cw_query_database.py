#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2021
# (c) University of Strathclyde 2021
# (c) James Hutton Institute 20201
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Interrogate the local CAZyme database"""


import logging
import json
import re

import pandas as pd

from datetime import datetime
from typing import List, Optional

from saintBioutils.utilities.logger import config_logger
from saintBioutils.utilities import file_io
from tqdm import tqdm

from cazy_webscraper import cazy_scraper
from cazy_webscraper.sql.sql_interface import get_selected_gbks, get_api_data
from cazy_webscraper.utilities.parsers import api_parser
from cazy_webscraper.utilities.parse_configuration import get_expansion_configuration


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy."""
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # Program preparation
    if argv is None:
        parser = api_parser.build_parser()
        args = parser.parse_args()
    else:
        parser = api_parser.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)

    connection, logger_name, cache_dir = cazy_scraper.connect_existing_db(args, time_stamp, start_time)

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "uniprot_data_retrieval"
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)

    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = get_expansion_configuration(args)

    # get the records of GenBank accessions matching the criteria of interest
    # {gbk_acc: gbk_id}
    gbk_dict = get_selected_gbks.get_genbank_accessions(
        class_filters,
        family_filters,
        taxonomy_filter_dict,
        kingdom_filters,
        ec_filters,
        connection,
    )

    query_data = get_query_data(gbk_dict)

    output_path = compile_output_name(args, time_stamp)
    
    if 'json' in args.file_types:
        json_output_path = output_path + ".json"
        with open(json_output_path, 'w') as fh:
            json.dump(query_data, fh)

    if 'csv' in args.file_types:
        write_csv_output(query_data, args)


def get_query_data(gbk_dict, connection, args):
    """Retrieve additional data as requested per retrieved GenBank accession
    
    :param gbk_dict: dict, {gbk_acc: gbk_id} GenBank accessions matching user criteria
    :param connection: open sqlaclchemy connection for an SQLite db
    :param args: cmd-line args parser

    Structure of dict if all data is retrievied:
    query_data = {
        genbank_accession: 
            'class': {CAZy classes},
            'family': {CAZy families},
            'subfamily': {CAZy subfamilies},
            'kingdom': kingdom,
            'genus': genus,
            'organism': genus species strain,
            'ec_numbers': {EC number annotations},
            'pdb_accessions': {PDB accessions},
            'uniprot_accession': UniProt protein accession,
            'uniprot_name': Name of the protein from UniProt,
            'uniprot_sequence': Protein Aa seq from UniProt,
            'uniprot_sequence_date': Date the seq was last modified in UniProt,
            'gbk_sequence': Protein Aa seq from GenBank
            'gbk_sequence_date': Date the seq was last modified in Gbk,
    }
    
    Retrun dict of retrieved data, keyed by GenBank accession and valued by dict containing
        requested data.
    """
    # create dict to store all data retrieved from the databsae query
    query_data = {}
    
    # add the GenBank accessions of interest to the query dict
    for gbk_acc in gbk_dict:
        query_data[gbk_acc] = {}

    if ('class' in args.include) or ('family' in args.include) or ('subfamily' in args.include):
        # retrieve the CAZy family annotations from the local CAZyme database
        query_data = get_api_data.get_class_fam_annotations(gbk_dict, query_data, connection, args)

    if ('kingdom' in args.include) or ('genus' in args.include) or ('organism' in args.include):
        # retrieve the taxonomy data from the local CAZyme database
        query_data = get_api_data.get_tax_annotations(gbk_dict, query_data, connection, args)

    if 'ec' in args.include:
       # retrieve the ec numbers from the local CAZyme database
       query_data = get_api_data.get_ec_annotations(gbk_dict, query_data, connection)

    if 'pdb' in args.include:
        # retrieve the PDB accessions from the local CAZyme database
        query_data = get_api_data.get_pdb_accessions(gbk_dict, query_data, connection)

    if ('uniprot_acc' in args.include) or ('uniprot_name' in args.include) or ('uniprot_seq' in args.include):
        # retrieve the UniProt data from the local CAZyme database
        query_data = get_api_data.get_uniprot_data(gbk_dict, query_data, connection, args)

    if 'genbank_seq' in args.include:
        # retrieve GenBank protein sequences from the local CAZyme database
        query_data = get_api_data.get_gbk_seq(gbk_dict, query_data, connection)


def write_csv_output(query_data, args, output_path, time_stamp):
    """Parse the output into a df structure and write out to a csv file.
    
    :param query_data: dict containing db query data.
    :param args: cmd-line args parser
    :param output_path: Path to write out output
    :param time_stamp: date and time script was invoked
    
    Return nothing.
    """
    output_path += ".csv"

    column_names = get_column_names(args)

    df_data = []  # list of nested lists, one nested list per df row

    for gbk_acc in tqdm(query_data, desc="Compiling output dataframe"):
        new_rows = []

        if (('class' not in args.include) and ('family' not in args.include) and ('subfamily' not in args.include)) or \
            (('class' in args.include) and ('family' not in args.include) and ('subfamily' not in args.include)) or \
                (('class' not in args.include) and ('family' in args.include) and ('subfamily' not in args.include)) \
                    (('class' not in args.include) and ('family' not in args.include) and ('subfamily' in args.include)):
                    new_rows.append([gbk_acc])  # don't need to create multiple rows to separate the class/fam/subfam annotations
        
        else:
            # create one row for each CAZy class, CAZy family and CAZy subfamily annotation
            # link the parent-child relationships between CAZy class, family and subfamily
            new_rows = get_class_fam_relationships(gbk_acc, query_data, args)

        # data from CAZy
        if 'kingdom' in args.include:
            for row in new_rows:
                row.append(query_data[gbk_acc]["kingdom"])
        if 'genus' in args.include:
            for row in new_rows:
                row.append(query_data[gbk_acc]["genus"])
        if 'organism' in args.include:
            for row in new_rows:
                row.append(query_data[gbk_acc]["organism"])
        
        # data from GenBank
        if "genbank_seq" in args.include:
            for row in new_rows:
                row.append(query_data[gbk_acc]["gbk_sequence"])
                row.append(query_data[gbk_acc]["gbk_sequence_date"])
        
        # data from UniProt
        if 'uniprot_acc' in args.include:
            for row in new_rows:
                row.append(query_data[gbk_acc]["uniprot_accession"])
        if 'uniprot_name' in args.include:
            for row in new_rows:
                row.append(query_data[gbk_acc]["uniprot_name"])
        if "ec" in args.include:
            for row in new_rows:
                row.append(" ".join(list(query_data[gbk_acc]["ec_numbers"])))
        if "pdb" in args.include:
            for row in new_rows:
                row.append(" ".join(list(query_data[gbk_acc]["pdb_accessions"])))
        if "uniprot_seq" in args.include:
            for row in new_rows:
                row.append(query_data[gbk_acc]["uniprot_sequence"])
                row.append(query_data[gbk_acc]["uniprot_sequence_date"])

        for row in new_rows:
            df_data.append(row)

    query_df = pd.DataFrame(df_data, columns=column_names)

    query_df.to_csv(output_path)


def compile_output_name(args, time_stamp):
    db_name = args.database.name.replace(".db","")
    output_path = args.output_dir / f"{db_name}_{time_stamp}"

    if args.cazy_class:
        output_path += "_classes"
    if args.cazy_family:
        output_path += "_fams"
    if args.cazy_subfamily:
        output_path += "_subfams"
    if args.kingdom:
        output_path += "_kngdm"
    if args.genus:
        output_path += "_genus"
    if args.organism:
        output_path += "_orgnsm"
    if args.ec:
        output_path += "_ec"
    if args.pdb:
        output_path += "_pdb"
    if args.uniprot:
        output_path += "_uniprot"
    if args.seq_uniprot:
        output_path += "_uniprotSeq"
    if args.seq_genbank:
        output_path += "_gbkSeq"

    return output_path


def get_column_names(args):
    """Compile the column names for the output df.
    
    :param args: cmd-line args parser
    
    Return list of column names"""
    column_names = ["genbank_accession"]

    # data from CAZy
    if 'class' in args.include:
        column_names.append("class")
    if 'family' in args.include:
        column_names.append("family")
    if 'subfamily' in args.include:
        column_names.append("subfamily")
    if 'kingdom' in args.include:
        column_names.append("kingdom")
    if 'genus' in args.include:
        column_names.append("genus")
    if 'organism' in args.include:
        column_names.append("source_organism")

    # data from GenBank
    if 'genbank_seq' in args.include:
        column_names.append("genbank_sequence")
        column_names.append("genbank_sequence_date")

    # data from UniProt
    if 'uniprot_acc' in args.include:
        column_names.append("uniprot_accession")
    if 'uniprot_name' in args.include:
        column_names.append("uniprot_name")
    if 'ec' in args.include:
        column_names.append("ec_number")
    if 'pdb' in args.include:
        column_names.append("pdb_accession")
    if 'uniprot_seq' in args.include:
        column_names.append("uniprot_sequence")
        column_names.append("uniprot_sequence_date")

    return column_names


def get_class_fam_relationships(gbk_acc, protein_query_data, args):
    """Define CAZy class-family-subfamily relationships in the query data.
    
    :param gbk_acc: str, GenBank protein accession
    :param query_data: dict, data retrieved from the local CAZyme database for the specific protein
        {
            'class': {CAZy classes},
            'family': {CAZy families},
            'subfamily': {CAZy subfamilies},
        }
    :param args: cmd-line args parser
    
    Return list of nested lists, one nested list per relationship.
    """
    logger =logging.getLogger(__name__)

    new_rows = []

    # class and family, NO subfamily
    if ('class' in args.include) and ('family' in args.include) and ('subfamily' not in args.include):
        families = protein_query_data['family']
        for family in families:
            try:
                parent_class = re.match(r'\D{2,3}', family).group()
            except AttributeError:
                logger.warning(f"Could not retrieve CAZy class from {family}, setting CAZy class as 'NA'")
                parent_class = 'NA'
            new_rows.append([gbk_acc, parent_class, family])

    # class and subfamily, NO family
    elif  'class' in args.include and 'family' not in args.include and 'subfamily' in args.include:
        subfamilies = protein_query_data['subfamily']
        for subfamily in subfamilies:
            try:
                parent_class = re.match(r'\D{2,3}', subfamily).group()
            except AttributeError:
                logger.warning(f"Could not retrieve CAZy class from {subfamily}, setting CAZy class as 'NA'")
                parent_class = 'NA'
            new_rows.append([gbk_acc, parent_class, subfamily])

    # family and subfamily, NO class
    elif  'class' not in args.include and 'family' in args.include and 'subfamily' in args.include:
        subfamilies = protein_query_data['subfamily']
        added_families = set()
        for subfamily in subfamilies:
            try:
                parent_fam = re.match(r'\D{2,3}\d+_', subfamily).group()[:-1]
                added_families.add(parent_fam)
            except AttributeError:
                logger.warning(f"Could not retrieve CAZy family from {subfamily}, setting CAZy family as 'NA'")
                parent_fam = 'NA'
            new_rows.append([gbk_acc, parent_fam, subfamily])

        # add remaining families for which there is no subfamily annotation
        families = protein_query_data['family']
        for family in families:
            if family not in added_families:
                new_rows.append([gbk_acc, family, 'NA'])

    # class, family and subfamily
    else:
        subfamilies = protein_query_data['subfamily']
        added_families = set()
        for subfamily in subfamilies:

            # get the parent family
            try:
                parent_fam = re.match(r'\D{2,3}', subfamily).group()
                added_families.add(parent_fam)
            except AttributeError:
                logger.warning(f"Could not retrieve CAZy family from {subfamily}, setting CAZy family as 'NA'")
                parent_fam = 'NA'

            # get the parent class
            try:
                parent_class = re.match(r'\D{2,3}', subfamily).group()
            except AttributeError:
                logger.warning(f"Could not retrieve CAZy class from {subfamily}, setting CAZy class as 'NA'")
                parent_class = 'NA'

            new_rows.append([gbk_acc, parent_class, parent_fam, subfamily])
        
        # add remaining families for which there is no subfamily annotation
        families = protein_query_data['family']

        for family in families:
            if family not in added_families:
                try:
                    parent_class = re.match(r'\D{2,3}', family).group()
                except AttributeError:
                    logger.warning(f"Could not retrieve CAZy class from {family}, setting CAZy class as 'NA'")
                    parent_class = 'NA'
                new_rows.append([gbk_acc, parent_class, family, 'NA'])

    return new_rows


if __name__ == "__main__":
    main()
