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
import sys

import pandas as pd

from datetime import datetime
from typing import List, Optional

from saintBioutils.utilities.logger import config_logger
from saintBioutils.utilities import file_io

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

    output_file_types = validate_file_types(args)

    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = get_expansion_configuration(args)

    output_dict = {}  # dict to store all output from interrogating the CAZyme db

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

    write_output(query_data, output_file_types, args)


def validate_file_types(args):
    """Validate the output files types chosen by the user
    
    :param args: cmd-line args parser
    
    Return list of output file types
    """
    logger = logging.getLogger(__name__)

    file_types = args.file_types.split(",")

    accepted_file_types = ['csv', 'json']
    
    chosen_types = [ftype.lower() for ftype in file_types if ftype.lower() in accepted_file_types]
    for ftype in file_types:
        if ftype.lower() in accepted_file_types:
            chosen_types.append(ftype.lower())
        else:
            logger.error(f'File type: {ftype} not supported by cazy_webscraper')

    if len(chosen_types) == 0:
        logger.error(
            "No accepted file types provided. Please chose at least one of the following:\n"
            "csv, json\n"
            "Terminating program"
        )
        sys.exit(1)
    
    return chosen_types


def get_query_data(gbk_dict, connection, args):
    """Retrieve additional data as requested per retrieved GenBank accession
    
    :param gbk_dict: dict, {gbk_acc: gbk_id} GenBank accessions matching user criteria
    :param connection: open sqlaclchemy connection for an SQLite db
    :param args: cmd-line args parser
    
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
        # retrieve GEnbank protein sequences from the local CAZyme database
        query_data = get_api_data.get_gbk_seq(gbk_dict, query_data, connection)


if __name__ == "__main__":
    main()
