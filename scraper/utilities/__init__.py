#!/usr/bin/env python
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
"""Module for building a logger and argument parser."""

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
        description="Programme to scrape the CAZy database",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add arguments to parser

    # Add option to specify path to configuration file
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        metavar="config file",
        default=None,
        help="Path to configuration file. Default: scrape entire database",
    )

    # Add option specify how data is split
    parser.add_argument(
        "-d",
        "--data_split",
        choices=["class", "family"],
        default=None,
        help=(
            (
                "Define if and how data is split."
                "Default: single dataframe and FASTA file"
            )
        ),
    )

    # Add option to force file over writting
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help=(
            (
                "Force writing in existing output directory."
                "Default: do not write in existing directory"
            )
        ),
    )

    # Add log file name option
    # If not given, no log file will be written out
    parser.add_argument(
        "-l",
        "--log",
        type=Path,
        metavar="log file name",
        default=None,
        help="Defines log file name and/or path. Default: no log file written",
    )

    # Add option to prevent over writing of existing files
    # and cause addition of files to output directory
    parser.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help=(
            (
                "Enable deleting files already present"
                "in the output directory."
                "Default: disabled"
            )
        ),
    )

    # Add option to specify output directory
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="output directory name",
        default=sys.stdout,
        help="Path to output directory",
    )

    # Add option to specify if retrieving subfamily data
    parser.add_argument(
        "-s",
        "--subfamily",
        dest="subfamily",
        action="store_true",
        default=False,
        help=(
            (
                "Enable retrieval of subfamily number."
                "Default: annotate proteins with their"
                "family not subfamily number"
            )
        ),
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

    # Setup console handler to log to terminal
    console_log_handler = logging.StreamHandler()
    if args.verbose is True:
        console_log_handler.setLevel(logging.INFO)
    else:
        console_log_handler.setLevel(logging.WARNING)
    console_log_handler.setFormatter(log_formatter)
    logger.addHandler(console_log_handler)

    # Setup file handler to log to a file
    if args.log is not None:
        file_log_handler = logging.FileHandler(args.log)
        if args.verbose is True:
            file_log_handler.setLevel(logging.INFO)
        else:
            file_log_handler.setLevel(logging.WARNING)
        file_log_handler.setFormatter(log_formatter)
        logger.addHandler(file_log_handler)

    return logger
