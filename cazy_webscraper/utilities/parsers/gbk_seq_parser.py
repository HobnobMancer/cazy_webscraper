#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
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
"""Submodule for building cmd-line parsers"""


import argparse
import sys

from pathlib import Path
from typing import List, Optional


def build_parser(argv: Optional[List] = None):
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

    # Add option for building a BLAST database of retrieved protein sequences
    parser.add_argument(
        "-b",
        "--blastdb",
        type=Path,
        default=None,
        help=(
            "Create BLAST database of retrieved GenBank protein sequences.\n"
            "Give the path to the directory to store the database"
        ),
    )

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

    # Add option to force file over writting
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force file over writting",
    )

    # Add option to specify families to retrieve protein sequences for
    parser.add_argument(
        "--families",
        type=str,
        default=None,
        help="Families to scrape. Separate families by commas 'GH1,GH2'",
    )

    # Add option to enable writing sequences to FASTA file or files, or not at all
    parser.add_argument(
        "--fasta",
        type=str,
        default=None,
        help=(
            "Enable writing out retrieved sequences to FASTA file(s).\n"
            "Writing 'separate' produces a single FASTA file per retrieved protein sequence,\n"
            "else, write the path to the FASTA to add retrieved protein sequences to\n"
            "(this can be a pre-existing or non-existing FASTA file."
        ),
    )

    # Add option to restrict the scrape to specific kingdoms
    parser.add_argument(
        "--kingdoms",
        type=str,
        default=None,
        help=(
            "Kingdoms to scrape. Separate by a single comma.\n"
            "Options= archaea, bacteria, eukaryota, viruses, unclassified (not case sensitive)"
        ),
    )

    # Add option to restrict scrape to specific genera
    parser.add_argument(
        "--genera",
        type=str,
        default=None,
        help="Genera to restrict the scrape to"
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

    # Add option to restrict the scrape to specific species. This will scrape CAZymes from
    # all strains belonging to each listed species
    parser.add_argument(
        "--species",
        type=str,
        default=None,
        help="Species (written as Genus Species) to restrict the scrape to"
    )

    # Add option to restrict scraping to specific strains of organisms
    parser.add_argument(
        "--strains",
        type=str,
        default=None,
        help=(
            "Specific strains of organisms to restrict the scrape to "
            "(written as Genus Species Strain)"
        ),
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