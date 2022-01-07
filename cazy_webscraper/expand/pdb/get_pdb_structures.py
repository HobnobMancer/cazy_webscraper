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
"""Retrieve PDB structures from RSCB PDB and write to disk"""


import logging
import os
import sys

import pandas as pd

from datetime import datetime
from typing import List, Optional

import Bio.PDB

from tqdm import tqdm

from cazy_webscraper import cazy_webscraper
from cazy_webscraper.expand import get_chunks_gen
from cazy_webscraper.sql.sql_interface import get_selected_pdbs
from cazy_webscraper.utilities import config_logger, file_io, parse_configuration
from cazy_webscraper.utilities.parsers import pdb_strctre_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)
    # parse cmd-line arguments
    if argv is None:
        parser = pdb_strctre_parser.build_parser()
        args = parser.parse_args()
    else:
        args = pdb_strctre_parser.build_parser(argv).parse_args()

    if logger is None:
        logger = logging.getLogger(__package__)
        config_logger(args)

    # validate PDB file choices
    logger.info("Checking valid file types were provided")
    validate_pdb_file_types(args)

    connection, logger_name, cache_dir = cazy_webscraper.connect_existing_db(args, time_stamp)

    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = parse_configuration.get_expansion_configuration(args)

    pdb_accessions = get_selected_pdbs.get_pdb_accessions(
        class_filters,
        family_filters,
        taxonomy_filter_dict,
        kingdom_filters,
        ec_filters,
        connection,
    )

    if len(pdb_accessions) == 0:
        logger.warning(
            "No PDB accessions matched the criteria provided.\n"
            "Retrieving no protein structure files from PDB"
        )
    else:
        logger.warning(f"Retrieving {len(pdb_accessions)} structure files from PDB")

    # make output and cache dirs
    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "uniprot_data_retrieval"
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)

    download_pdb_structures(pdb_accessions, cache_dir, args)

    cache_path = cache_dir / f"pdb_retrieval_{time_stamp}.txt"
    with open(cache_path, 'a') as fh:
        for acc in pdb_accessions:
            fh.write(f"{acc}\n")

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    if args.verbose:
        logger.info(
            "Finished getting structure files from PDB\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_webscraper.VERSION_INFO}\n"
            f"Citation: {cazy_webscraper.CITATION_INFO}"
        )
    else:
        print(
            "=====================cazy_webscraper=====================\n"
            "Finished getting structure files from PDB\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_webscraper.VERSION_INFO}\n"
            f"Citation: {cazy_webscraper.CITATION_INFO}"
        )


def validate_pdb_file_types(args):
    """Check valid file types for PDB were specified by the user.
    :param args: cmd-line args parser
    Return nothing.
    """
    logger = logging.getLogger(__name__)

    valid_choices = ['mmCif', 'pdb', 'xml', 'mmtf', 'bundle']

    invalid_choices = []

    user_choices = (args.pdb).split(",")

    for choice in user_choices:
        if choice not in valid_choices:
            invalid_choices.append(choice)

    if len(invalid_choices) != 0:
        logger.error(
            f"Invalid file option selected: {invalid_choices}.\n"
            f"The valid choices are: {valid_choices}.\n"
            "Terminating program."
        )
        sys.exit(1)

    return


def download_pdb_structures(pdb_accessions, args):
    """Download protein structure from the RSCB PDB database

    :param pdb_accession: list of PDB accessions
    :param args: cmd-line args parser

    Return nothing.
    """
    pdbl = Bio.PDB.PDBList()

    logger = logging.getLogger(__name__)
    logger.warning("Starting downloading of structure files from PDB")

    if args.outdir is None:
        logger.warning("Downloading to current working directory")
        for accession_list in get_chunks_gen(pdb_accessions, args.batch_limit):
            for file_type in tqdm((args.pdb).split(","), desc="Downloading"):
                pdbl.download_pdb_files(
                    pdb_codes=accession_list,
                    file_format=file_type,
                    overwrite=args.overwrite,
                )

    else:
        logger.warning(f"Downloading structures to {args.outdir}")
        for accession_list in get_chunks_gen(pdb_accessions, args.batch_limit):
            for file_type in tqdm((args.pdb).split(","), desc="Downloading"):
                pdbl.download_pdb_files(
                    pdb_codes=accession_list,
                    file_format=file_type,
                    overwrite=args.overwrite,
                    pdir=args.outdir,
                )

    return


if __name__ == "__main__":
    main()
