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

# Bio.PDB reference:
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
# implemented in Python. Bioinformatics 19: 2308–2310
"""Retrieve PDB structures from RSCB PDB and write to disk"""


import logging
import os
import time

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio.PDB import PDBList
from tqdm import tqdm

from scraper import file_io
from scraper.sql.sql_orm import (
    CazyFamily,
    Pdb,
    get_db_session,
)
from scraper.utilities import config_logger, build_pdb_structures_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # parse cmd-line arguments
    if argv is None:
        parser = build_pdb_structures_parser()
        args = parser.parse_args()
    else:
        args = build_pdb_structures_parser(argv).parse_args()

    # build logger
    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)

    # create session to local database
    if os.path.isfile(args.database) is False:
        logger.error(
            "Could not find local CAZy database. Check path is correct. Terminating programme."
        )
    session = get_db_session(args)

    # create output directory
    if args.outdir is None:
        # save structure files to the cwd
        outdir = os.getcwd()
    else:
        outdir = args.outdir
    file_io.make_output_directory(outdir, args.force, args.nodelete)

    # check if any classes or families were specified to retrieve the sequences only for them
    if (args.classes is None) and (args.families is None) and (args.config is None):
        config_dict = None

    else:
        # create dictionary of CAZy classes/families to retrieve sequences for
        file_io_path = file_io.__file__
        excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
            file_io_path,
            args,
        )
        # excluded_classes and cazy_dict are not needed, but required when the func is 
        # called in cazy_webscraper.py

    # retrieve protein structures from PDB

    # retrieve sequences for all CAZymes
    if config_dict is None:
        get_every_cazymes_structures(outdir, session, args)

    # retrieve sequences for specific CAZy classes and/or families
    else:
        get_structures_for_specific_cazymes(outdir, config_dict, session, args)

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    end_time = pd.to_datetime(start_time)
    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    logger.info(
        "Finished dowloading protein structures from PDB."
        "Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}"
    )

    print(
        "=====================cazy_webscraper-expand-pdb_structures=====================\n"
        "Finished dowloading protein structures from PDB."
        "Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}\n"
    )


def get_every_cazymes_structures(outdir, session, args):
    """Get PDB accessions in the db, and coordinate downloading the structure from pdb.

    :param outdir: path to output directory
    :param session: open SQLite db session
    :param args: cmd-line argument parser

    Return nothing.
    """
    # retrieve structures for only primary PDB accessions
    if args.primary is True:
        pdb_query = session.query(Pdb).filter(Pdb.primary == True).all()

    # retrieve structures for all PDB accessions
    else:
        pdb_query = session.query(Pdb).all()

    for query_result in pdb_query:
        pdb_accession = query_result.pdb_accession
        pdb_accession = pdb_accession[:pdb_accession.find("[")]
        download_pdb_structures(pdb_accession, outdir, args)

    return


def get_structures_for_specific_cazymes(outdir, config_dict, session, args):
    """Retrieve primary PDB accessions for CAZymes meeting criteria in the config_dict.

    :param outdir: path to output directory
    :param config_dict: dict, defines CAZy classes and families to retrieve accessions from
    :param session: open SQLite db session
    :param args: cmd-line argument parser

    Return nothing.
    """
    # start with the classes
    if len(config_dict["classes"]) != 0:
        # create a dictionary to convert full class name to abbreviation
        cazy_classes = config_dict["classes"]

        for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):
            # retrieve class name abbreviation
            cazy_class = cazy_class[((cazy_class.find("(")) + 1):((cazy_class.find(")")) - 1)]

            # Retrieve PDB accessions under the current working CAZy class
            if args.primary:
                class_query = session.query(Pdb).\
                    filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
                    filter(Pdb.primary == True).\
                    all()

            else:
                class_query = session.query(Pdb).\
                    filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
                    all()

            for query_result in class_query:
                pdb_accession = query_result.pdb_accession
                pdb_accession = pdb_accession[:pdb_accession.find("[")]
                download_pdb_structures(pdb_accession, outdir, args)

    # retrieve protein structure for specified families
    for key in config_dict:
        if key == "classes":
            continue
        if len(config_dict[key]) is None:
            continue

        for family in tqdm(config_dict[key], desc=f"Parsing families in {key}"):
            if family.find("_") != -1:  # subfamily

                if args.primary:
                    family_query = session.query(Pdb).\
                        filter(CazyFamily.subfamily == family).\
                        filter(Pdb.primary == True).\
                        all()
                else:
                    family_query = session.query(Pdb).\
                        filter(CazyFamily.subfamily == family).\
                        all()

            else:  # family

                if args.primary:
                    family_query = session.query(Pdb).\
                        filter(CazyFamily.subfamily == family).\
                        filter(Pdb.primary == True).\
                        all()
                else:
                    family_query = session.query(Pdb).\
                        filter(CazyFamily.subfamily == family).\
                        all()

            for query_result in family_query:
                pdb_accession = query_result.pdb_accession
                pdb_accession = pdb_accession[:pdb_accession.find("[")]
                download_pdb_structures(pdb_accession, outdir, args)

    return


def download_pdb_structures(pdb_accession, outdir, args):
    """Download protein structure from the RSCB PDB database

    :param pdb_accession: str, accession of record in the PDB database
    :param outdir: path to output directory
    :param args: cmd-line args parser

    Return nothing.
    """
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(f"{pdb_accession}", file_format=args.pdb, pdir=args.outdir)
    time.sleep(2)  # to prevent bombarding the system

    return


if __name__ == "__main__":
    main()
