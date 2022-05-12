#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
#
# Bio.PDB reference:
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
# implemented in Python. Bioinformatics 19: 2308–2310
"""Retrieve the accessions for proteins of interest, and store accessions in the local db"""


import logging

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio import Entrez
from tqdm import tqdm
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger

from cazy_webscraper import closing_message, connect_existing_db
from cazy_webscraper.sql.sql_interface import get_selected_gbks, get_table_dicts
from cazy_webscraper.expand import get_chunks_list
from cazy_webscraper.utilities import parse_configuration
from cazy_webscraper.utilities.parsers.extract_seq_parser import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # parse cmd-line arguments
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        args = build_parser(argv).parse_args()

    if logger is None:
        logger = logging.getLogger(__package__)
        config_logger(args)

    Entrez.email = args.email

    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)

    logger.info(f"Connected to local db: {args.database}")

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "genome_accessions"
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    
    logger.info(f"Using cache dir: {cache_dir}")

    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = parse_configuration.get_expansion_configuration(args)

    gbk_table_dict = get_table_dicts.get_gbk_table_dict(connection)
    # {genbank_accession: 'taxa_id': str, 'gbk_id': int}

    # build dick {gbk_acc: db_id} matching the users specified criteria
    # either via a list in a file or parameters provided via config file and/or command line

    gbk_dict = {}  # {gbk_acc: gbk_id}

    if len(list(gbk_dict.keys())) == 0:
        gbk_dict = get_selected_gbks.get_genbank_accessions(
            class_filters,
            family_filters,
            taxonomy_filter_dict,
            kingdom_filters,
            ec_filters,
            connection,
        )
    
    # retrieve genome accessions from NCBI

    closing_message("extract_sequences", start_time, args)
    

def get_genomic_accessions(protein_ids, args):
    """Retrieve genomic accessions for a set of proteins from NCBI.
    
    :param protein_ids: list, GenBank protein IDs
    :param args: cmd-line args parser
    
    Return dict {protein_id: {'genbank': str, 'refseq': str}}
    """
    logger = logging.getLogger(__name__)

    logger.info("Starting retrieval of genomic accessions")

    failed_batches = []  # store lists that causes issues

    genome_dict = {}  # used for storing retrieved genomic accessions

    # break up long list into smaller batches
    batches = get_chunks_list(protein_ids, args.batch_size)

    for batch in tqdm(batches, desc="Batch quering NCBI"):
        try:
            #
            logger.info()
        except RuntimeError:
            failed_batches.append(batch)

    