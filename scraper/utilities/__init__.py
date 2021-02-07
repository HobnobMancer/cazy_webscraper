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
"""Build cmd-line parser and logger."""

import argparse
import logging
import sys

from pathlib import Path
from typing import List, Optional


def build_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="cazy_webscraper.py",
        description="Scrapes the CAZy database",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add optional arguments to parser

    # Add option to specify path to configuration file
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        metavar="config file",
        default=None,
        help="Path to configuration file. Default: None, scrapes entire database",
    )

    # Add option to define complete classes to scrape
    parser.add_argument(
        "--classes",
        type=str,
        default=None,
        help="Classes from which all families are to be scraped. Separate classes by ','"
    )

    # Add option to provide a path to an existing SQL database to add newly scraped data to
    parser.add_argument(
        "-d",
        "--database",
        type=Path,
        metavar="local database path",
        default=None,
        help="path to an existing local CAZy SQL database",
    )

    # Add option to force file over writting
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force file over writting",
    )

    # Add option to specify families to scrape
    parser.add_argument(
        "--families",
        type=str,
        default=None,
        help="Families to scrape. Separate families by commas 'GH1,GH2'"
    )

    # Add option to download FASTA file for protein from GenBank
    parser.add_argument(
        "-g",
        "--genbank",
        type=str,
        metavar="Email address of user",
        default=None,
        help="Enable FASTA files from GenBank, and user email required for Entrez",
    )

    # Add log file name option
    # If not given, no log file will be written out
    parser.add_argument(
        "-l",
        "--log",
        type=Path,
        metavar="log file name",
        default=None,
        help="Defines log file name and/or path",
    )

    # Add option to prevent over writing of existing files
    # and cause addition of files to output directory
    parser.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help="enable/disable deletion of exisiting files",
    )

    # Add option to specify output directory to write output dataframes to
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="output file name",
        default=sys.stdout,
        help="Output filename",
    )

    # Add option to specify ouput directory for writing out fasta files from GenBank to
    parser.add_argument(
        "-genbank_output",
        type=Path,
        metavar="output file name",
        default=sys.stdout,
        help="Output filename",
    )

    # Add option to specift output directory for writing out PDB structure files to
    parser.add_argument(
        "-pdb_output",
        type=Path,
        metavar="output file name",
        default=None,
        help="Output filename",
    )

    # Add option to download FASTA file for protein from GenBank
    parser.add_argument(
        "-p",
        "--pdb",
        choices=[None, "mmCif", "pdb", "xml", "mmtf", "bundle"],
        type=str,
        default=None,
        help="Enable downloading of protein structures in XXXX format from PDB",
    )

    # Add option to enable number of times to retry scraping
    parser.add_argument(
        "-r",
        "--retries",
        type=int,
        default=0,
        help="Number of times to retry scraping a family or class page if error encountered",
    )

    # Add option to enable retrieval of subfamilies
    parser.add_argument(
        "-s",
        "--subfamilies",
        dest="subfamilies",
        action="store_true",
        default=False,
        help="Enable retrieval of subfamilies from CAZy",
    )

    # Add option for more detail (verbose) logging
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Set logger level to 'INFO'",
    )

    if argv is None:
        # parse command-line
        return parser
    else:
        # return namespace
        return parser.parse_args(argv)


def config_logger(args) -> logging.Logger:
    """Configure package wide logger.

    Configure a logger at the package level, from which the module will inherit.
    If CMD-line args are provided, these are used to define output streams, and
    logging level.

    :param args: cmd-line args parser

    Return nothing
    """
    logger = logging.getLogger(__package__)

    # Set format of loglines
    log_formatter = logging.Formatter("[%(levelname)s] [%(name)s]: %(message)s")

    # define logging level
    if args.verbose is True:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    # Setup console handler to log to terminal
    console_log_handler = logging.StreamHandler()
    console_log_handler.setFormatter(log_formatter)
    logger.addHandler(console_log_handler)

    # Setup file handler to log to a file
    if args.log is not None:
        file_log_handler = logging.FileHandler(args.log)
        file_log_handler.setFormatter(log_formatter)
        logger.addHandler(file_log_handler)

    return


def build_genbank_sequences_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for the script 'expand.genbank_sequences.py'."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="genbank_sequences.py",
        description="Populates local CAZy database with protein sequences from GenBank",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add positional/required arguments
    parser.add_argument(
        "database",
        type=Path,
        metavar="local CAZy database",
        help="Path to local CAZy database",
    )

    parser.add_argument(
        "email",
        type=str,
        metavar="user email address",
        help="User email address, requirement of NCBI-Entrez",
    )

    # Add optional arguments to parser

    # Add option to specify path to configuration file
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        metavar="config file",
        default=None,
        help="Path to configuration file. Default: None, scrapes entire database",
    )

    # Add option to define classes to retrieve protein sequences for
    parser.add_argument(
        "--classes",
        type=str,
        default=None,
        help="Classes from which all families are to be scraped. Separate classes by ','"
    )

    # specify the number of accessions posted in single ePost to NCBI
    parser.add_argument(
        "-e",
        "--epost",
        type=int,
        default=150,
        help="Number of accessions posted to NCBI per epost, advice to be max 200. Default=150"
    )

    # Add option to specify families to retrieve protein sequences for
    parser.add_argument(
        "-f",
        "--families",
        type=str,
        default=None,
        help="Families to scrape. Separate families by commas 'GH1,GH2'"
    )

    # Add log file name option
    # If not given, no log file will be written out
    parser.add_argument(
        "-l",
        "--log",
        type=Path,
        metavar="log file name",
        default=None,
        help="Defines log file name and/or path",
    )

    # enable retrieving protein sequences for only primary GenBank accessions
    parser.add_argument(
        "-p",
        "--primary",
        dest="primary",
        action="store_true",
        default=False,
        help="Enable retrieveing protein sequences for only primary GenBank accessions",
    )

    # Add option to update sequences if the retrieved sequence is different
    # If not enabled then sequences will only be retrieved and added for proteins that do not
    # already have a protein sequence
    parser.add_argument(
        "-u",
        "--update",
        dest="update",
        action="store_true",
        default=False,
        help="Enable overwriting sequences in the database if the retrieved sequence is different",
    )

    # Add option to specify output directory to write output dataframes to
    parser.add_argument(
        "-w",
        "--write",
        type=Path,
        metavar="path to FASTA file dire",
        default=None,
        help="Enable writing out protein sequences to a FASTA file",
    )

    # Add option for more detail (verbose) logging
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Set logger level to 'INFO'",
    )

    if argv is None:
        # parse command-line
        return parser
    else:
        # return namespace
        return parser.parse_args(argv)


def build_pdb_structures_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for the script 'expand.genbank_sequences.py'."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="pdb_structures.py",
        description="Download structures from PDB",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add positional/required arguments
    parser.add_argument(
        "database",
        type=Path,
        metavar="local CAZy database",
        help="Path to local CAZy database",
    )

    parser.add_argument(
        "pdb",
        choices=["mmCif", "pdb", "xml", "mmtf", "bundle"],
        type=str,
        help="File format of downloaded structure from PDB",
    )

    # Add optional arguments to parser

    # Add option to specify path to configuration file
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        metavar="config file",
        default=None,
        help="Path to configuration file. Default: None, scrapes entire database",
    )

    # Add option to define classes to retrieve protein sequences for
    parser.add_argument(
        "--classes",
        type=str,
        default=None,
        help="Classes from which all families are to be scraped. Separate classes by ','"
    )

    # Add option to specify families to retrieve protein sequences for
    parser.add_argument(
        "--families",
        type=str,
        default=None,
        help="Families to scrape. Separate families by commas 'GH1,GH2'"
    )

    # enable force writing in an existing directory
    parser.add_argument(
        "-f"
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force file over writting",
    )

    # Add log file name option
    # If not given, no log file will be written out
    parser.add_argument(
        "-l",
        "--log",
        type=Path,
        metavar="log file name",
        default=None,
        help="Defines log file name and/or path",
    )

    # Add option to prevent over writing of existing files
    # and cause addition of files to output directory
    parser.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help="enable/disable deletion of exisiting files",
    )

    # enable specifying an output directory
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        metavar="output directory path",
        help="Path to output directory to which downloaded structures are retrieved",
    )

    # enable retrieving protein sequences for only primary GenBank accessions
    parser.add_argument(
        "-p",
        "--primary",
        dest="primary",
        action="store_true",
        default=False,
        help="Enable retrieveing protein sequences for only primary GenBank accessions",
    )

    # Add option to update sequences if the retrieved sequence is different
    # If not enabled then sequences will only be retrieved and added for proteins that do not
    # already have a protein sequence
    parser.add_argument(
        "-u",
        "--update",
        dest="update",
        action="store_true",
        default=False,
        help="Enable overwriting sequences in the database if the retrieved sequence is different",
    )

    # Add option for more detail (verbose) logging
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Set logger level to 'INFO'",
    )

    if argv is None:
        # parse command-line
        return parser
    else:
        # return namespace
        return parser.parse_args(argv)
