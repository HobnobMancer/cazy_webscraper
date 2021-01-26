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
"""Retrieve PDB structures from RSCB PDB and write to disk"""


import logging
import time

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio.PDB import PDBList
from tqdm import tqdm

from scraper import file_io
from scraper.sql.sql_orm import (
    Cazyme,
    CazyFamily,
    Pdb,
    get_db_session,
)
from scraper.utilities import build_logger, build_pdb_structures_parser


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
    logger = build_logger("expand.pdb_structures", args)

    # create session to local database
    session = get_db_session(args, logger)

    # check if any classes or families were specified to retrieve the sequences only for them
    if (args.classess is not None) and (args.families is not None):
        # create dictionary of CAZy classes/families to retrieve sequences for
        file_io_path = file_io.__file__
        cazy_dict, std_class_names = file_io.get_cazy_dict_std_names(file_io_path, logger)
        config_dict = file_io.get_cmd_defined_fams_classes(cazy_dict, std_class_names, args, logger)

    else:
        config_dict = None

    if config_dict is None:
        # get sequences for everything
        get_every_cazymes_structures(session, args, logger)

    else:
        # get sequences for only specified classes/families
        if args.primary is True:
            get_specific_proteins_structures_primary_only(config_dict, session, args, logger)
        else:
            get_specific_proteins_structures(config_dict, session, args, logger)

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


def get_every_cazymes_structures(session, args, logger):
    """Get PDB accessions in the db, and coordinate downloading the structure from pdb.

    :param session: open SQLite db session
    :param args: cmd-line argument parser
    :param logger: logger object

    Return nothing.
    """
    # retrieve structures for only primary PDB accessions
    if args.primary is True:
        pdb_query = session.query(Pdb).filter(Pdb.primary == True).all()

    # retrieve structures for all PDB accessions
    else:
        pdb_query = session.query(Pdb).all()
    
    for query_result in tqdm(pdb_query, desc="Downloading structures from PDB"):
        pdb_accession = query_result.pdb_accession
        download_pdb_structures(pdb_accession, args, logger)
    
    return


def get_specific_proteins_structures_primary_only(config_dict, session, args, logger):
    """Retrieve primary PDB accessions for CAZymes meeting criteria in the config_dict.

    :param config_dict: dict, defines CAZy classes and families to retrieve accessions from
    :param session: open SQLite db session
    :param args: cmd-line argument parser
    :param logger: logger object

    Return nothing.
    """
    # start with the classes
    if len(config_dict["classes"]) != 0:
        # create a dictionary to convert full class name to abbreviation
        cazy_classes = config_dict["classes"]
        for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):
            # retrieve class name abbreviation
            cazy_class = cazy_class[cazy_class.find("("):cazy_class.find(")")]

            # retrieve all GenBank accessions catalogued in the CAZy class
            class_query = session.query(Pdb).\
                filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
                filter(Pdb.primary == True).\
                all()

            for query_result in tqdm(class_query, desc=f"Parsing families in {cazy_class}"):
                pdb_accession = query_result.pdb_accession
                download_pdb_structures(pdb_accession, args, logger)

    # retrieve protein structure for specified families
    for key in config_dict:
        if key == "classes":
            continue

        for family in tqdm(config_dict[key], desc=f"Parsing families in {key}"):
            if family.find("_") != -1:  # subfamily
                # Retrieve PDB accessions catalogued under the subfamily
                family_query = (Pdb).\
                    filter(CazyFamily.subfamily == family).\
                    filter(Pdb.primary == True).\
                    all()

            else:  # family
                # Retrieve PDB accessions catalogued under the family
                family_query = (Pdb).\
                    filter(CazyFamily.family == family).\
                    filter(Pdb.primary == True).\
                    all()

            for query_result in tqdm(family_query, desc=f"Parsing PDB accessions in {family}"):
                pdb_accession = query_result.pdb_accession
                download_pdb_structures(pdb_accession, args, logger)

    return


def get_specific_proteins_structures(config_dict, session, args, logger):
    """Retrieve all PDB accessions for CAZymes meeting criteria in the config_dict.

    :param config_dict: dict, defines CAZy classes and families to retrieve accessions from
    :param session: open SQLite db session
    :param args: cmd-line argument parser
    :param logger: logger object

    Return nothing.
    """
    # start with the classes
    if len(config_dict["classes"]) != 0:
        # create a dictionary to convert full class name to abbreviation
        cazy_classes = config_dict["classes"]
        for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):
            # retrieve class name abbreviation
            cazy_class = cazy_class[cazy_class.find("("):cazy_class.find(")")]

            # retrieve all GenBank accessions catalogued in the CAZy class
            class_query = session.query(Pdb).\
                filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
                all()

            for query_result in tqdm(class_query, desc=f"Parsing families in {cazy_class}"):
                pdb_accession = query_result.pdb_accession
                download_pdb_structures(pdb_accession, args, logger)

    # retrieve protein structure for specified families
    for key in config_dict:
        if key == "classes":
            continue

        for family in tqdm(config_dict[key], desc=f"Parsing families in {key}"):
            if family.find("_") != -1:  # subfamily
                # Retrieve PDB accessions catalogued under the subfamily
                family_query = (Pdb).\
                    filter(CazyFamily.subfamily == family).\
                    all()

            else:  # family
                # Retrieve PDB accessions catalogued under the family
                family_query = (Pdb).\
                    filter(CazyFamily.family == family).\
                    all()

            for query_result in tqdm(family_query, desc=f"Parsing PDB accessions in {family}"):
                pdb_accession = query_result.pdb_accession
                download_pdb_structures(pdb_accession, args, logger)

    return


def download_pdb_structures(pdb_accession, args, logger):
    """Download protein structure from the RSCB PDB database

    :param pdb_accession: str, accession of record in the PDB database
    :param args: cmd-line argument parser
    :param logger: logger object

    Return nothing.
    """
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(f"{pdb_accession}", file_format=args.pdb, pdir=args.outdir)

    return
