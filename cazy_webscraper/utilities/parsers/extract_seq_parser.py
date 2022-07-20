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
import sys

from pathlib import Path
from typing import List, Optional


class ValidateNames(argparse.Action):
    """Check the user has provided valid database names."""
    def __call__(self, parser, args, values, option_string=None):
        valid_formats = ("genbank", "uniprot")
        invalid = False

        for value in values:
            if value.lower() not in valid_formats:
                invalid = True
                raise ValueError(f'Invalid source "{value.lower()}" provided. Accepted sources: {valid_formats}')

        if invalid:
            sys.exit(1)

        parsed_values = [_.lower() for _ in values]
        setattr(args, self.dest, parsed_values)


def build_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for the script 'expand.genbank_sequences.py'."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="extract_sequences.py",
        description="Extract GenBank and/or UniProt protein seqs and write to FASTA and/or BLAST db",
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
        "source",
        nargs='+',
        action=ValidateNames,
        choices=["genbank", "uniprot"],
        type=str,
        help="Original db source of protein seqs (genbank and/or uniprot)",
    )

    # Add optional arguments to parser
    parser.add_argument(
        "--genbank_accessions",
        type=Path,
        default=None,
        help="Path to text file contining GenBank accessions",
    )
    parser.add_argument(
        "--uniprot_accessions",
        type=Path,
        default=None,
        help="Path to text file contining UniProt accessions",
    )

    # Add option for building a BLAST database of retrieved protein sequences
    parser.add_argument(
        "-b",
        "--blastdb",
        type=Path,
        default=None,
        help=(
            "Create BLAST database of extracted protein sequences.\n"
            "Provide the path to the directory to store the database"
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

    # Add option to define classes to retrieve protein sequences for
    parser.add_argument(
        "--classes",
        type=str,
        default=None,
        help="Classes from which all families are to be scraped. Separate classes by ','"
    )

    parser.add_argument(
        "--ec_filter",
        type=str,
        default=None,
        help="Limit retrieval to proteins annotated with the provided EC numbers. Separate EC numbers with single commas"
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
        "--fasta_dir",
        type=Path,
        default=None,
        help=(
            "Write out each extracted sequence to a separate FASTA file in the provided dir\n"
            "Provide a path to a dir to write out FASTA files"
        ),
    )

    parser.add_argument(
        "--fasta_file",
        type=Path,
        default=None,
        help=(
            "Write out all extracted sequences to a single FASTA file.\n"
            "Provide a path to write out the FASTA to"
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

    parser.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help="Do not delete content in existing output dir",
    )

    parser.add_argument(
        "--nodelete_cache",
        dest="nodelete_cache",
        action="store_true",
        default=False,
        help="Do not delete content in existing cache dir",
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

    parser.add_argument(
        "--sql_echo",
        dest="sql_echo",
        action="store_true",
        default=False,
        help="Set SQLite engine echo to True (SQLite will print its log messages)",
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
