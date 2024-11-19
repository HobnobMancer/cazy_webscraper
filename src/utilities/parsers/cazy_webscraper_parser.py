#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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

    parser.add_argument(
        "email",
        type=str,
        nargs='?',
        help="User email address. Requirement of Entrez, used to get source organsism data. Email is not stored be cazy_webscraper."
    )

    # Add optional arguments to parser

    # Add option to specify path to configuration file
    parser.add_argument(
        "--cache_dir",
        type=Path,
        default=None,
        help="Target path for cache dir to be used instead of default path",
    )

    # Add option to use a pre-downloaded CAZy txt file
    parser.add_argument(
        "--cazy_data",
        type=Path,
        default=None,
        help="Path predownloaded CAZy txt file",
    )

    # Add option to use own CAZy class synoymn dict
    parser.add_argument(
        "--cazy_synonyms",
        type=Path,
        default=None,
        help="Path to JSON file containing CAZy class synoymn names",
    )

    # Add option to define complete classes to scrape
    parser.add_argument(
        "--classes",
        type=str,
        default=None,
        help="Classes from which all families are to be scraped. Separate classes by ','"
    )

    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        metavar="config file",
        default=None,
        help="Path to configuration file. Default: None, scrapes entire database",
    )

    # Add option to display citation
    parser.add_argument(
        "-C",
        "--citation",
        dest="citation",
        action="store_true",
        default=False,
        help="Print cazy_webscraper citation message",
    )

    parser.add_argument(
        "-o",
        "--db_output",
        type=Path,
        default=None,
        help="Target output path to build new SQL database",
    )

    parser.add_argument(
        "-d",
        "--database",
        type=Path,
        default=None,
        help="Path to an existing local CAZy SQL database",
    )

    parser.add_argument(
        "--delete_old_relationships",
        dest="delete_old_relationships",
        action="store_true",
        default=False,
        help=(
            "Delete old GenBank accession - CAZy family relationships (annotations)\n"
            "that are in the local db but are not in CAZy, e.g. when CAZy has moved a\n"
            "protein from one fam to another, delete the old family annotation."
        ),
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
        help="Families to scrape. Separate families by commas 'GH1,GH2' (case sensitive)"
    )

    # Add option to restrict scrape to specific genera
    parser.add_argument(
        "--genera",
        type=str,
        default=None,
        help="Genera to restrict the scrape to"
    )

    parser.add_argument(
        "--kingdoms",
        type=str,
        default=None,
        help="Tax Kingdoms to restrict the scrape to"
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

    parser.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help="When called, content in the existing out dir is NOT deleted",
    )

    parser.add_argument(
        "--ncbi_batch_size",
        type=int,
        default=200,
        help="Number of genbank accessions in each NCBI Taxonomy db batch query"
    )

    # Add option to not delete content in the existing cache dir
    parser.add_argument(
        "--nodelete_cache",
        dest="nodelete_cache",
        action="store_true",
        default=False,
        help="When called, content in the existing cache dir is NOT deleted",
    )

    parser.add_argument(
        "--nodelete_log",
        dest="nodelete_log",
        action="store_true",
        default=False,
        help="When called, content in the existing log dir is NOT deleted",
    )

    # Add option to enable number of times to retry scraping
    parser.add_argument(
        "-r",
        "--retries",
        type=int,
        default=10,
        help="Number of times to retry scraping a family or class page if error encountered",
    )

    parser.add_argument(
        "--skip_ncbi_tax",
        dest="skip_ncbi_tax",
        action="store_true",
        default=False,
        help=(
            "Skip retrieving the latest tax classification from the NCBI Taxonomy db for proteins\n"
            "listed with multiple taxs in CAZy.\n"
            "For these proteins the first taxonomy listed in CAZy is added to the local CAZyme db"
        ),
    )

    # Add option to force file over writting
    parser.add_argument(
        "--sql_echo",
        dest="sql_echo",
        action="store_true",
        default=False,
        help="Set SQLite engine echo to True (SQLite will print its log messages)",
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

    # Add option to define time out limit for trying to connect to CAZy
    parser.add_argument(
        "-t",
        "--timeout",
        type=int,
        default=45,
        help="Connection timeout limit (seconds)"
    )

    parser.add_argument(
        "--validate",
        dest="validate",
        action="store_true",
        default=False,
        help=(
            "Retrieve CAZy fam population sizes from CAZy and use to check\n"
            "the number of family members added to the local database"
        ),
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
    
    # Add option to display version
    parser.add_argument(
        "-V",
        "--version",
        dest="version",
        action="store_true",
        default=False,
        help="Print cazy_webscraper version number",
    )    

    if argv is None:
        # parse command-line
        return parser
    else:
        # return namespace
        return parser.parse_args(argv)
