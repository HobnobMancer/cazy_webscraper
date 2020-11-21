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

    # Add option on how to split data
    parser.add_argument(
        "-d",
        "--data_split",
        choices=[None, "class", "family"],
        type=str,
        default=None,
        help="How data is to be split. Not split=all, split by class=class, split by family=family",
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

    # Add option to download FASTA file for protein from GenBank
    parser.add_argument(
        "-g",
        "--genbank",
        type=str,
        help="Email address of user"
        default=False,
        help="Enable downloading of protein sequence in FASTA format from GenBank",
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
        "--genbank-output",
        type=Path,
        metavar="output file name",
        default=sys.stdout,
        help="Output filename",
    )

    # Add option to specift output directory for writing out PDB structure files to
    parser.add_argument(
        "--pdb-output",
        type=Path,
        metavar="output file name",
        default=sys.stdout,
        help="Output filename",
    )

    # Add option to download FASTA file for protein from GenBank
    parser.add_argument(
        "-p",
        "--pdb",
        dest="pdb",
        action="store_true",
        default=False,
        help="Enable downloading of protein structures in XXXX format from PDB",
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


def build_logger(script_name, args) -> logging.Logger:
    """Return a logger for this script.
    Enables logger for script, sets parameters and creates new file to store log.
    :param script_name: str, name of script
    :param args: parser argument
    Return logger object.
    """
    logger = logging.getLogger(script_name)

    # Set format of loglines
    log_formatter = logging.Formatter(
        script_name + ": {} - {}".format("%(asctime)s", "%(message)s")
    )

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

    return logger
