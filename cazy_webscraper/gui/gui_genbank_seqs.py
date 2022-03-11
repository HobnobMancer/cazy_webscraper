#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Create GUI for get_genbank_seqs.py."""

import argparse

from datetime import datetime
from pathlib import Path

from gooey import Gooey, GooeyParser

from cazy_webscraper.expand.genbank import get_genbank_sequences
from cazy_webscraper.gui.assets import build_menus

cw_menu = build_menus(
    'cw_get_genbank_seqs',
    'Retrieve protein sequences from GenBank afor CAZymes in a local CAZyme database. The retrieved protein sequences are stored in the local CAZyme database.'
)

@Gooey(
    program_name="Get GenBank Protein Sequences",
    image_dir="cazy_webscraper/gui/assets/get_genbank_seqs",
    menu=cw_menu,
)
def main():
    # Create parser object
    parser = GooeyParser(
        prog="get_genbank_seqs.py",
        description="Retrieve protein sequences from GenBank and store the sequences in the local CAZyme database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "database",
        metavar="Local CAZyme database",
        widget="FileChooser",
        help="The path to the local CAZyme database to extract protein sequences from",
    )

    parser.add_argument(
        "email",
        metavar="Email address",
        type=str,
        help=(
            "This is a requirement of NCBI Entrez, which is used to get source "
            "organsism data. The email is not stored be cazy_webscraper."
        ),
    )

    # Add optional arguments to parser

    #
    # ACCESSIONS
    #

    accessions_group = parser.add_argument_group(
        "Use accessions", 
        "Provide a list  or lists of protein accessions"
    )

    accessions_group.add_argument(
        "--genbank_accessions",
        metavar="GenBank accessions",
        widget="FileChooser",
        default=None,
        help="Path to text file contining GenBank accessions",
    )

    #
    # CLASS AND FAM FILTERS
    #

    class_group = parser.add_argument_group(
        "Class and family filters", 
        "Define CAZy classes and CAZy families to scrape"
    )

    class_group.add_argument(
        "--classes",
        metavar="CAZy classes",
        type=str,
        default=None,
        help=(
            "Classes from which all families are to be scraped. "
            "Separate classes with a single comma ','"
        ),
    )

    class_group.add_argument(
        "--families",
        metavar="CAZy (sub)families",
        type=str,
        default=None,
        help="Families and subfamilies to scrape. Separate families by commas 'GH1,GH2'. CAZy families are case sensitive"
    )

    #
    # TAX filters
    #

    tax_group = parser.add_argument_group(
        "Taxonomy filters",
        "These are applied after CAZy class and CAZy family filters",
    )

    tax_group.add_argument(
        "--kingdoms",
        metavar="Taxonomy kingdoms",
        type=str,
        default=None,
        help="Tax Kingdoms to restrict the scrape to"
    )

    # Add option to restrict scrape to specific genera
    tax_group.add_argument(
        "--genera",
        metavar="Genera",
        type=str,
        default=None,
        help="Genera to restrict the scrape to. Separate genera with a single comma"
    )

    tax_group.add_argument(
        "--species",
        metavar="Species",
        type=str,
        default=None,
        help="Species (written as Genus Species) to restrict the scrape to. Separate species with a single comma"
    )

    tax_group.add_argument(
        "--strains",
        metavar="Strains",
        type=str,
        default=None,
        help=(
            "Specific strains of organisms to restrict the scrape to "
            "(written as Genus Species Strain). Separate strains with a single comma."
        ),
    )

    #
    # ADDITIONAL FILTERS
    #

    add_filt_group = parser.add_argument_group(
        "Additionally Options", 
        "Further refine the protein sequences retrieved from the local database"
    )

    add_filt_group.add_argument(
        "--ec_filter",
        type=str,
        default=None,
        help="Limit retrieval to proteins annotated with the provided EC numbers. Separate EC numbers with a single comma, e.g. 1.2.3.4,2.3.4.5. The 'EC' prefix is option, and '*' and *-* are accepted for missing digits"
    )

    #
    # CACHE OPTIONS
    #

    cache_group = parser.add_argument_group(
        "Cache Options", 
        "Use cache files and change the cache location"
    )

    cache_group.add_argument(
        "--seq_dict",
        metavar="Cached sequences",
        widget="FileChooser",
        default=None,
        help="Path to a JSON file containing the cached GenBank protein sequences from a previous run. The file is keyed by GenBank accessions and valued by protein sequences",
    )

    # Add option to specify path to configuration file
    cache_group.add_argument(
        "--cache_dir",
        metavar="Cache directory",
        widget="DirChooser",
        default=None,
        help="Target path for cache dir to be used instead of default path",
    )

    cache_group.add_argument(
        "--nodelete_cache",
        metavar="Do no delete existing cache",
        dest="nodelete_cache",
        action="store_true",
        default=False,
        help="When called, content in the existing cache dir is NOT deleted",
    )

    #
    # LOG OPTIONS
    #

    log_group = parser.add_argument_group("Logging Options")

    log_group.add_argument(
        "-l",
        "--log",
        type=Path,
        metavar="Log file name",
        default=None,
        help="Define directory to write out log files",
    )

    log_group.add_argument(
        "-v",
        "--verbose",
        metavar="Verbose logging",
        dest="verbose",
        action="store_true",
        default=False,
        help="Set logger level to 'INFO'",
    )

    #
    # MISC OPTIONS
    #

    misc_group = parser.add_argument_group("Misc Options")

    misc_group.add_argument(
        "--seq_update",
        dest="seq_update",
        metavar="Update local sequences",
        action="store_true",
        default=False,
        help="Enable overwriting sequences in the database if the retrieved sequence is different",
    )

    misc_group.add_argument(
        "--batch_size",
        type=int,
        default=150,
        help="Size of batch queries sent to NCBI Entrez"
    )

    misc_group.add_argument(
        "-r",
        "--retries",
        metavar="Connection retries",
        type=int,
        default=10,
        help="Number of times to retry the connection to NCBI Entrez if the connection fails",
    )

    # Add option to force file over writting
    misc_group.add_argument(
        "--sql_echo",
        metavar="SQL db echo",
        dest="sql_echo",
        action="store_true",
        default=False,
        help="Set SQLite engine echo to True (SQLite will print its log messages)",
    )

    misc_group.add_argument(
        "--cazy_synonyms",
        metavar="CAZy class synonyms",
        widget="FileChooser",
        default=None,
        help=(
            "Path to JSON file containing CAZy class synoymn names "
            "Use your own CAZy class synonyms"
        ),
    )

    misc_group.add_argument(
        "-c",
        "--config",
        metavar="Configuration file",
        widget="FileChooser",
        default=None,
        help="Path to configuration file. Default: None, scrapes entire database",
    )

    gooey_args = parser.parse_args()

    get_genbank_sequences.main(args=gooey_args)


if __name__ == "__main__":
    main()
