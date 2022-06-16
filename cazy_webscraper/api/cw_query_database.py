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
import sys

import pandas as pd

from datetime import datetime
from pathlib import Path
from typing import List, Optional

from saintBioutils.utilities.logger import config_logger
from saintBioutils.utilities import file_io
from tqdm import tqdm

from cazy_webscraper import cazy_scraper, closing_message
from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks, get_api_data
from cazy_webscraper.utilities.parsers.api_parser import build_parser
from cazy_webscraper.utilities.parse_configuration import get_expansion_configuration


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy."""
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # Program preparation
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)

    connection, logger_name, cache_dir = cazy_scraper.connect_existing_db(args, time_stamp, start_time)
    logger.info("Open connection to local cazyme database:", args.database)

    if args.output_dir is not None:
        file_io.make_output_directory(args.output_dir, args.force, args.nodelete)

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        logger.info("Building cache dir")
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        logger.info("Building cache dir")
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

    output_path = compile_output_name(args)

    existing_files = ""
    if 'json' in args.file_types:
        json_output_path = str(output_path) + ".json"
        logger.warning(f"JSON output path: {json_output_path}")
        if Path(json_output_path).exists():
            existing_files = existing_files + " " + f"{json_output_path}\n"
    
    if 'csv' in args.file_types:
        csv_output_path = str(output_path) + ".csv"
        logger.warning(f"CSV output path: {csv_output_path}")
        if Path(csv_output_path).exists():
            existing_files = existing_files + " " + f"{csv_output_path}\n"
    
    existing_files = existing_files.strip()
    if len(existing_files) != 0:
        if args.overwrite:
            logger.warning(
                "The output files\n"
                f"{existing_files}"
                "Exist already. Overwrite is True. Overwriting output files"
            )
        else:
            logger.warning(
                "The output files\n"
                f"{existing_files}"
                "Exist already. Overwrite is False\n"
                "To overwrite the files use the --overwrite flag, or "
                "change the output file prefix using the --prefix flag\n"
                "Terminating program"
            )
            closing_message("cw_query_database", start_time, args)
            sys.exit(1)

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

    query_data = get_query_data(gbk_dict, connection, args)
    logger.warning(f"Retrieved {len(list(query_data.keys()))} proteins from the local db")
    
    if 'json' in args.file_types:
        write_json_output(json_output_path, query_data, args)

    if 'csv' in args.file_types:
        write_csv_output(query_data, args, csv_output_path)
    
    closing_message("cw_query_database", start_time, args)


def compile_output_name(args):
    """Compile the file name for output files.
    
    :param args: cmd-line args parser
    
    Return str"""
    logger = logging.getLogger(__name__)
    logger.info("Compiling names of output files")

    file_prefix = f"{args.database.name.replace('.db', '')}"
    if args.prefix is not None:
        file_prefix = f"{args.prefix}_{file_prefix}"

    file_prefix += "_gbkAcc"
    
    # CAZy classification info
    if 'class' in args.include:
        file_prefix += "_classes"
    if 'family' in args.include:
        file_prefix += "_fams"
    if 'subfamily' in args.include:
        file_prefix += "_subfams"
    
    # tax infor
    if 'kingdom' in args.include:
        file_prefix += "_kngdm"
    if 'genus' in args.include:
        file_prefix += "_genus"
    if 'organism' in args.include:
        file_prefix += "_orgnsm"

    # data from genbank
    if "genbank_seq" in args.include:
        file_prefix += "_gbkSeq"

    # data from uniprot
    if 'uniprot_acc' in args.include:
        file_prefix += "_uni_acc"
    if 'uniprot_name' in args.include:
        file_prefix += "_uni_name"
    if "ec" in args.include:
        file_prefix += "_ec"
    if "pdb" in args.include:
        file_prefix += "_pdb"
    if "uniprot_seq" in args.include:
        file_prefix += "_uniprotSeq"

    if args.output_dir is not None:
        output_path = args.output_dir / file_prefix
    else:
        output_path = file_prefix
    
    logger.info("Compiled output path: ", output_path)
    
    return output_path


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
    logger = logging.getLogger(__name__)

    # create dict to store all data retrieved from the databsae query
    query_data = {}
    
    # add the GenBank accessions of interest to the query dict
    for gbk_acc in gbk_dict:
        query_data[gbk_acc] = {}

    if ('class' in args.include) or ('family' in args.include) or ('subfamily' in args.include):
        logger.info("Retrieving CAZy family annotations")
        # retrieve the CAZy family annotations from the local CAZyme database
        query_data = get_api_data.get_class_fam_annotations(gbk_dict, query_data, connection, args)

    if ('kingdom' in args.include) or ('genus' in args.include) or ('organism' in args.include):
        logger.info("Retrieving taxonomy data")
        # retrieve the taxonomy data from the local CAZyme database
        query_data = get_api_data.get_tax_annotations(gbk_dict, query_data, connection, args)

    if 'ec' in args.include:
        logger.info("Retrieving EC numbers")
        # retrieve the ec numbers from the local CAZyme database
        query_data = get_api_data.get_ec_annotations(gbk_dict, query_data, connection)

    if 'pdb' in args.include:
        logger.info("Retrieving PDB accessions")
        # retrieve the PDB accessions from the local CAZyme database
        query_data = get_api_data.get_pdb_accessions(gbk_dict, query_data, connection)

    if ('uniprot_acc' in args.include) or ('uniprot_name' in args.include) or ('uniprot_seq' in args.include):
        logger.info("Retrieving UniProt data")
        # retrieve the UniProt data from the local CAZyme database
        query_data = get_api_data.get_uniprot_data(gbk_dict, query_data, connection, args)

    if 'genbank_seq' in args.include:
        logger.info("Retrieving GenBank protein sequences")
        # retrieve GenBank protein sequences from the local CAZyme database
        query_data = get_api_data.get_gbk_seq(gbk_dict, query_data, connection)

    return query_data


def write_json_output(json_output_path, query_data, args):
    """Parse dict to be suitable for JSON serialisation and write out output.
    
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

    Return nothing
    """
    logger = logging.getLogger(__name__)
    logger.info("Writing JSON output")

    set_keys = [
        'class',
        'family',
        'subfamily',
        'ec_numbers',
        'pdb_accessions',
    ]
    for genbank_accession in query_data:
        for key in set_keys:
            try:
                query_data[genbank_accession][key] = list(query_data[genbank_accession][key])
            
            except KeyError:  # raised when data was not retrieved from the db
                pass

    with open(json_output_path, 'w') as fh:
        json.dump(query_data, fh)


def write_csv_output(query_data, args, output_path):
    """Parse the output into a df structure and write out to a csv file.
    
    :param query_data: dict containing db query data.
    :param args: cmd-line args parser
    :param output_path: Path to write out output
    
    Return nothing.
    """
    logger = logging.getLogger(__name__)
    logger.info("Writing CSV output")

    column_names = get_column_names(args)

    df_data = []  # list of nested lists, one nested list per df row

    for gbk_acc in tqdm(query_data, desc="Compiling output dataframe"):
        new_rows = []

        if (('class' not in args.include) and ('family' not in args.include) and ('subfamily' not in args.include)):
            new_rows.append([gbk_acc])  # don't need to create multiple rows to separate the class/fam/subfam annotations
        
        else:
            # create one row for each CAZy class, CAZy family and CAZy subfamily annotation
            # link the parent-child relationships between CAZy class, family and subfamily
            new_rows = get_class_fam_relationships(gbk_acc, query_data, args)

        # data from CAZy
        if 'kingdom' in args.include:
            new_rows = add_single_value_to_rows(query_data, gbk_acc, 'kingdom', new_rows)
        if 'genus' in args.include:
            new_rows = add_single_value_to_rows(query_data, gbk_acc, 'genus', new_rows)
        if 'organism' in args.include:
            new_rows = add_single_value_to_rows(query_data, gbk_acc, 'organism', new_rows)
        
        # data from GenBank
        if "genbank_seq" in args.include:
            try:
                seq = query_data[gbk_acc]["gbk_sequence"]
            except KeyError:
                seq = None
            
            try:
                date = query_data[gbk_acc]["gbk_sequence_date"]
            except KeyError:
                date = None
    
            for row in new_rows:
                row.append(seq)
                row.append(date)
        
        # data from UniProt
        if 'uniprot_acc' in args.include:
            new_rows = add_single_value_to_rows(query_data, gbk_acc, 'uniprot_accession', new_rows)

        if 'uniprot_name' in args.include:
            new_rows = add_single_value_to_rows(query_data, gbk_acc, 'uniprot_name', new_rows)

        if "ec" in args.include:
            try:
                ec_numbers = " ".join(list(query_data[gbk_acc]["ec_numbers"]))
            except KeyError:
                ec_numbers = None
            for row in new_rows:
                row.append(ec_numbers)
            
        if "pdb" in args.include:
            try:
                pdb_accessions = " ".join(list(query_data[gbk_acc]["pdb_accessions"]))
            except KeyError:
                pdb_accessions = None
            for row in new_rows:
                row.append(pdb_accessions)
            
        if "uniprot_seq" in args.include:
            try:
                seq = query_data[gbk_acc]["uniprot_sequence"]
            except KeyError:
                seq = None
            
            try:
                date = query_data[gbk_acc]["uniprot_sequence_date"]
            except KeyError:
                date = None
    
            for row in new_rows:
                row.append(seq)
                row.append(date)

        for row in new_rows:
            df_data.append(row)

    query_df = pd.DataFrame(df_data, columns=column_names)

    query_df.to_csv(output_path)


def add_single_value_to_rows(query_data, gbk_acc, key, new_rows):
    """Retrieve value for key from query data for the given protein, add value to all new_rows.
    
    :param query_data: dict
    :param gbk_acc: str
    :param key: str
    :param new_rows: list of lists
    
    Return new_rows
    """
    try:
        value = query_data[gbk_acc][key]
    except KeyError:
        value = None
    for row in new_rows:
        row.append(value)
    
    return new_rows


def get_column_names(args):
    """Compile the column names for the output df.
    
    :param args: cmd-line args parser
    
    Return list of column names"""
    logger = logging.getLogger(__name__)
    logger.info("Compiling column names for the CSV output")

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

    # class only
    if ('class' in args.include) and ('family' not in args.include) and ('subfamily' not in args.include):
        logging.info("Only retrieving CAZy classes annotations")
        classes = protein_query_data[gbk_acc]['class']
        for cazy_class in classes:
            new_rows.append([gbk_acc, cazy_class])

    # family only
    elif ('class' not in args.include) and ('family' in args.include) and ('subfamily' not in args.include):
        logger.info("Only retrieving CAZy family annotations")
        families = protein_query_data[gbk_acc]['family']
        for family in families:
            new_rows.append([gbk_acc, family])

    # subfamily only
    elif ('class' not in args.include) and ('family' not in args.include) and ('subfamily' in args.include):
        logger.info("Retrieving only CAZy subfamily annotations")
        subfamilies = protein_query_data[gbk_acc]['subfamily']
        for subfamily in subfamilies:
            new_rows.append([gbk_acc, subfamily])

    # class and family, NO subfamily
    elif ('class' in args.include) and ('family' in args.include) and ('subfamily' not in args.include):
        logger.info("Retrieving CAZy class and family annotations")
        families = protein_query_data[gbk_acc]['family']
        for family in families:
            try:
                parent_class = re.match(r'\D{2,3}', family).group()
            except AttributeError:
                logger.warning(f"Could not retrieve CAZy class from {family}, setting CAZy class as 'NA'")
                parent_class = 'NA'
            new_rows.append([gbk_acc, parent_class, family])

    # class and subfamily, NO family
    elif  ('class' in args.include) and ('family' not in args.include) and ('subfamily' in args.include):
        logger.info("Retrieving CAZy class and subfamily annotations")
        subfamilies = protein_query_data[gbk_acc]['subfamily']
        for subfamily in subfamilies:
            try:
                parent_class = re.match(r'\D{2,3}', subfamily).group()
            except AttributeError:
                logger.warning(f"Could not retrieve CAZy class from {subfamily}, setting CAZy class as 'NA'")
                parent_class = 'NA'
            new_rows.append([gbk_acc, parent_class, subfamily])

    # family and subfamily, NO class
    elif  ('class' not in args.include) and ('family' in args.include) and ('subfamily' in args.include):
        logger.info("Retrieving CAZy family and subfamily annotations")
        subfamilies = protein_query_data[gbk_acc]['subfamily']
        added_families = set()
        for subfamily in subfamilies:
            try:
                parent_fam = subfamily.split("_")[0]
                added_families.add(parent_fam)
            except AttributeError:
                logger.warning(f"Could not retrieve CAZy family from {subfamily}, setting CAZy family as 'NA'")
                parent_fam = 'NA'
            new_rows.append([gbk_acc, parent_fam, subfamily])

        # add remaining families for which there is no subfamily annotation
        families = protein_query_data[gbk_acc]['family']
        for family in families:
            if family not in added_families:
                new_rows.append([gbk_acc, family, None])

    # class, family and subfamily
    else:
        logger.info("Retrieving CAZy class, family and subfamily annotations")
        subfamilies = protein_query_data[gbk_acc]['subfamily']
        added_families = set()
        for subfamily in subfamilies:

            # get the parent family
            try:
                parent_fam = subfamily.split("_")[0]
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
        families = protein_query_data[gbk_acc]['family']

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
