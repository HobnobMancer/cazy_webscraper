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
        prog="cw_query_database",
        description="Interrogate a local CAZyme database",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "database",
        type=Path,
        help="Path to local CAZyme database"
    )

    parser.add_argument(
        "file_types",
        action="store",
        nargs="+",
        choices=["csv", "json"],
        help="File types to write the query output in [csv,json]"
    )

    parser.add_argument(
        "--cache_dir",
        type=Path,
        default=None,
        help="Target path for cache dir to be used instead of default path",
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
        help="CAZy classes to retrieve UniProt data for. Separate classes by ','"
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

    parser.add_argument(
        "--ec_filter",
        type=str,
        default=None,
        help="Limit retrieval to proteins annotated with the provided EC numbers. Separate EC numbers with single commas"
    )

    # Add option to specify families to scrape
    parser.add_argument(
        "--families",
        type=str,
        default=None,
        help="CAZy families to UniProt data for. Separate families by commas 'GH1,GH2' (case sensitive)"
    )

    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force writing to existing output dir",
    )

    # Add option to restrict scrape to specific genera
    parser.add_argument(
        "--genera",
        type=str,
        default=None,
        help="Genera to UniProt data for"
    )

    parser.add_argument(
        "--include",
        action="store",
        nargs="+",
        choices=["class", "family", "subfamily", "kingdom", "organism", "uniprot_acc", "uniprot_name", "ec", "pdb", "genbank_seq", "uniprot_seq"],
        help="Additional data to include in the output file. Separate with a single space (' ')"
    )

    parser.add_argument(
        "--kingdoms",
        type=str,
        default=None,
        help="Tax Kingdoms to UniProt data for"
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
        help="When called, content in the existing output dir is NOT deleted",
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
        "-o",
        "--output_dir",
        type=Path,
        default=None,
        help="Path to output dir, default: None (writes to cwd)",
    )

    parser.add_argument(
        "--overwrite",
        dest="overwrite",
        action="store_true",
        default=False,
        help="When called, overwrites existing output files",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default=None,
        help="Str to prefix all output files with, default: None",
    )

    # Add option to force file over writting
    parser.add_argument(
        "--sql_echo",
        dest="sql_echo",
        action="store_true",
        default=False,
        help="Set SQLite engine echo to True (SQLite will print its log messages)",
    )

    # Add option to UniProt data for specific species. This will scrape CAZymes from
    # all strains belonging to each listed species
    parser.add_argument(
        "--species",
        type=str,
        default=None,
        help="Species (written as Genus Species) to UniProt data for"
    )

    # Add option to restrict scraping to specific strains of organisms
    parser.add_argument(
        "--strains",
        type=str,
        default=None,
        help=(
            "Specific strains of organisms to UniProt data for "
            "(written as Genus Species Strain)"
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

    if argv is None:
        # parse command-line
        return parser
    else:
        # return namespace
        return parser.parse_args(argv)
