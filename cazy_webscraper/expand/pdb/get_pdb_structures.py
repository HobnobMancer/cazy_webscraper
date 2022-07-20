#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
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
import sys

import pandas as pd

from datetime import datetime
from typing import List, Optional

import Bio.PDB
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger

from tqdm import tqdm

from cazy_webscraper import closing_message, connect_existing_db
from cazy_webscraper.expand import get_chunks_gen
from cazy_webscraper.sql.sql_interface.get_data.get_table_dicts import (
    get_gbk_table_dict,
    get_uniprot_table_dict,
)
from cazy_webscraper.sql.sql_interface.get_data.get_selected_pdbs import get_pdb_accessions
from cazy_webscraper.sql.sql_interface.get_data.get_records import (
    get_user_genbank_sequences,
    get_user_uniprot_sequences
)
from cazy_webscraper.utilities import parse_configuration
from cazy_webscraper.utilities.parsers.pdb_strctre_parser import build_parser


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

    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)

    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = parse_configuration.get_expansion_configuration(args)

    gbk_dict = {}  # {gbk_acc: gbk_id}

    gbk_table_dict = get_gbk_table_dict(connection)
    # {genbank_accession: 'taxa_id': str, 'gbk_id': int}

    if args.genbank_accessions is not None:
        logger.warning(
            "Retrieving PDB structures for GenBank accessions "
            f"listed in {args.genbank_accessions}"
        )
        gbk_dict.update(get_user_genbank_sequences(gbk_table_dict, args))

    if args.uniprot_accessions is not None:
        logger.warning(
            "Extracting protein sequences for UniProt accessions "
            f"listed in {args.uniprot_accessions}"
        )
        uniprot_table_dict = get_uniprot_table_dict(connection)
        gbk_dict.update(get_user_uniprot_sequences(gbk_table_dict, uniprot_table_dict, args))

    pdb_accessions = get_pdb_accessions(
        class_filters,
        family_filters,
        taxonomy_filter_dict,
        kingdom_filters,
        ec_filters,
        gbk_table_dict,
        connection,
    )

    if len(pdb_accessions) == 0:
        logger.warning(
            "No PDB accessions matched the criteria provided.\n"
            "Retrieving no protein structure files from PDB"
        )
        sys.exit(1)
    else:
        logger.warning(f"Retrieving {len(pdb_accessions)} structure files from PDB")

    # make output and cache dirs
    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "pdb_retrieval"
        make_output_directory(cache_dir, args.force, args.nodelete_cache)

    download_pdb_structures(pdb_accessions, args)

    cache_path = cache_dir / f"pdb_retrieval_{time_stamp}.txt"
    with open(cache_path, 'a') as fh:
        for acc in pdb_accessions:
            fh.write(f"{acc}\n")

    closing_message("Get PDB structure files", start_time, args)


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
        for accession_list in get_chunks_gen(pdb_accessions, args.batch_size):
            for file_type in tqdm(args.pdb, desc="Downloading"):
                pdbl.download_pdb_files(
                    pdb_codes=accession_list,
                    file_format=file_type,
                    overwrite=args.overwrite,
                )

    else:
        logger.warning(f"Downloading structures to {args.outdir}")
        for accession_list in get_chunks_gen(pdb_accessions, args.batch_size):
            for file_type in tqdm(args.pdb, desc="Downloading"):
                pdbl.download_pdb_files(
                    pdb_codes=accession_list,
                    file_format=file_type,
                    overwrite=args.overwrite,
                    pdir=args.outdir,
                )

    return


if __name__ == "__main__":
    main()
