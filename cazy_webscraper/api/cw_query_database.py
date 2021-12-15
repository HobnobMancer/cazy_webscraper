#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2021
# (c) University of Strathclyde 2021
# (c) James Hutton Institute 20201
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
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

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

import pandas as pd

from datetime import datetime
from typing import List, Optional

from saintBioutils.utilities.logger import config_logger
from saintBioutils.utilities import file_io
from tqdm import tqdm

from cazy_webscraper import cazy_scraper
from cazy_webscraper.sql.sql_interface import get_selected_gbks, get_api_data
from cazy_webscraper.utilities.parsers import query_db_parser
from cazy_webscraper.utilities.parse_configuration import get_expansion_configuration


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy."""
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # Program preparation
    if argv is None:
        parser = query_db_parser.build_parser()
        args = parser.parse_args()
    else:
        parser = query_db_parser.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)

    connection, logger_name, cache_dir = cazy_scraper.connect_existing_db(args, time_stamp)

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

    write_output(query_data, args)


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

    if args.cazy_class or args.cazy_family or args.cazy_subfamily:
        # retrieve the CAZy family annotations from the local CAZyme database
        query_data = get_api_data.get_class_fam_annotations(gbk_dict, query_data, connection, args)

    if args.kingdom or args.genus or args.organism:
        # retrieve the taxonomy data from the local CAZyme database
        query_data = get_api_data.get_tax_annotations(gbk_dict, query_data, connection, args)

    if args.ec:
       # retrieve the ec numbers from the local CAZyme database
       query_data = get_api_data.get_ec_annotations(gbk_dict, query_data, connection)

    if args.pdb:
        # retrieve the PDB accessions from the local CAZyme database
        query_data = get_api_data.get_pdb_accessions(gbk_dict, query_data, connection)

    if args.uniprot or args.seq_uniprot:
        # retrieve the UniProt data from the local CAZyme database
        query_data = get_api_data.get_uniprot_data(gbk_dict, query_data, connection, args)

    if args.seq_genbank:
        # retrieve GenBank protein sequences from the local CAZyme database
        query_data = get_api_data.get_gbk_seq(gbk_dict, query_data, connection)


def write_output(query_data, args, time_stamp):
    """Write out query data to disk.
    
    :param query_data: dict containing db query data.
    :param args: cmd-line args parser
    
    Return nothing"""
    output_path = compile_output_name(args, time_stamp)

    if 'json' in args.file_types:
        output_path += ".json"
        with open(output_path, 'w') as fh:
            json.dump(query_data, fh)

    if 'csv' in args.file_types:
        output_path += ".csv"

        # compile pandas df of the data
        column_names = get_column_names(args)
        
        df_data = []  # list of nested lists, one nested list per df row

        for gbk_acc in tqdm(query_data, desc="Compiling output dataframe"):
            new_rows = []

            # create one row for each CAZy class, CAZy family and CAZy subfamily annotation
            # link the parent-child relationships between CAZy class, family and subfamily
            if args.cazy_class is False and args.cazy_family is False and args.cazy_subfamily is False:
                new_rows.append([gbk_acc])  # don't need to create multiple rows to separate the class/fam/subfam annotations

            if args.kingdom:
                for row in new_rows:
                    row.append(query_data[gbk_acc]["kingdom"])
            if args.genus:
                for row in new_rows:
                    row.append(query_data[gbk_acc]["genus"])
            if args.organism:
                for row in new_rows:
                    row.append(query_data[gbk_acc]["organism"])
            if args.ec:
                for row in new_rows:
                    row.append(" ".join(list(query_data[gbk_acc]["ec_numbers"])))
            if args.pdb:
                for row in new_rows:
                    row.append(" ".join(list(query_data[gbk_acc]["pdb_accessions"])))
            if args.uniprot:
                for row in new_rows:
                    row.append(query_data[gbk_acc]["uniprot_accession"])
                    row.append(query_data[gbk_acc]["uniprot_name"])
            if args.seq_uniprot:
                for row in new_rows:
                    row.append(query_data[gbk_acc]["uniprot_sequence"])
                    row.append(query_data[gbk_acc]["uniprot_sequence_date"])
            if args.seq_genbank:
                for row in new_rows:
                    row.append(query_data[gbk_acc]["gbk_sequence"])
                    row.append(query_data[gbk_acc]["gbk_sequence_date"])

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

    if args.cazy_class:
        column_names.append("class")
    if args.cazy_family:
        column_names.append("family")
    if args.cazy_subfamily:
        column_names.append("subfamily")
    if args.kingdom:
        column_names.append("kingdom")
    if args.genus:
        column_names.append("genus")
    if args.organism:
        column_names.append("source_organism")
    if args.ec:
        column_names.append("ec_number")
    if args.pdb:
        column_names.append("pdb_accession")
    if args.uniprot:
        column_names.append("uniprot_accession")
        column_names.append("uniprot_name")
    if args.seq_uniprot:
        column_names.append("uniprot_sequence")
        column_names.append("uniprot_sequence_date")
    if args.seq_genbank:
        column_names.append("genbank_sequence")
        column_names.append("genbank_sequence_date")

    return column_names


if __name__ == "__main__":
    main()
