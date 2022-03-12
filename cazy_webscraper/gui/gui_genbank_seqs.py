#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
    'Retrieve protein sequences from GenBank for CAZymes in a local CAZyme database. The retrieved protein sequences are stored in the local CAZyme database.'
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
        help="The path to the local CAZyme database",
    )

    parser.add_argument(
        "email",
        metavar="Email address",
        type=str,
        help=(
            "This is a requirement of NCBI Entrez. The email is not stored be cazy_webscraper."
        ),
    )

    # Add optional arguments to parser
    
    #
    # UPDATE
    #

    update_group = parser.add_argument_group(
        "Update Data",
        "Enable updating data in the local CAZyme database if a more recent version of the data is available")
    )
    
    update_group.add_argument(
        "--seq_update",
        dest="seq_update",
        metavar="Update local sequences",
        action="store_true",
        default=False,
        help="Enable overwriting sequences in the database if a most recent version of the sequence is retrieved",
    )
    
    #
    # ACCESSIONS
    #

    accessions_group = parser.add_argument_group(
        "Use Accessions", 
        "Provide a list or lists of protein accessions to explicitly define the proteins to retrieve data for. The provided protein accessions will be used instead of the filters to define proteins of interest"
    )

    accessions_group.add_argument(
        "--genbank_accessions",
        metavar="GenBank accessions",
        widget="FileChooser",
        default=None,
        help="Path to text file contining GenBank accessions of the proteins of interset",
    )

    #
    # CLASS AND FAM FILTERS
    #

    class_group = parser.add_argument_group(
        "Class and Family Filters", 
        "Retrieve data for proteins from specific CAZy classes, families and subfamilies"
    )

    class_group.add_argument(
        "--classes",
        metavar="CAZy classes",
        type=str,
        default=None,
        help=(
            "Classes from which all families will be scraped. "
            "Separate classes with a single comma ',', e.g. 'GH,GT,PL'"
        ),
    )

    class_group.add_argument(
        "--families",
        metavar="CAZy (sub)families",
        type=str,
        default=None,
        help="Families and subfamilies to scrape. Separate families with single commas 'GH1,GH2'. CAZy families are case sensitive"
    )

    class_group.add_argument(
        "-s",
        "--subfamilies",
        metavar="Retrieve subfamilies",
        dest="subfamilies",
        action="store_true",
        default=False,
        help="Enable retrieving of subfamily annotations from CAZy",
    )

    #
    # TAX filters
    #

    tax_group = parser.add_argument_group(
        "Taxonomy Filters",
        "Limit the retrieval of data to proteins derived from specific taxonomies. These are applied after CAZy class and CAZy family filters",
    )

    tax_group.add_argument(
        "--kingdoms",
        metavar="Taxonomy kingdoms",
        type=str,
        default=None,
        help="Retrieve data for proteins derived from organisms from specific kingdoms. Separate kingdoms with a sinlge comma, e.g. 'bacteria,eukaryota'. Excepted kingdoms are archaea, bacteria, eukaryota, viruses, and unclassified. Kingdoms are not case sensitive"
    )

    # Add option to restrict scrape to specific genera
    tax_group.add_argument(
        "--genera",
        metavar="Genera",
        type=str,
        default=None,
        help="Retrieve data for proteins sourced from organisms belong to specific genera. Separate genera with a single comma. A direct string match is used, therefore capitalise the first letter of each genera"
    )

    tax_group.add_argument(
        "--species",
        metavar="Species",
        type=str,
        default=None,
        help="Retrieve data for proteins soruced from specific species (written as Genus Species). Separate species with a single comma. A direct string match is used to select the species of interest"
    )

    tax_group.add_argument(
        "--strains",
        metavar="Strains",
        type=str,
        default=None,
        help=(
            "Retrieve data for proteins from specific strains of organisms "
            "(written as Genus Species Strain). Separate strains with a single comma. A direct string match is used to select organisms of interest"
        ),
    )

    #
    # ADDITIONAL FILTERS
    #

    add_filt_group = parser.add_argument_group(
        "Additionally Options", 
        "Further refine selecting proteins of interest for whom data will be retrieved"
    )

    add_filt_group.add_argument(
        "--ec_filter",
        metavar="EC numbers",
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
        help="Directory to be used for storing the cache instead of the default cache directory",
    )

    cache_group.add_argument(
        "--nodelete_cache",
        metavar="Do not delete existing cache",
        dest="nodelete_cache",
        action="store_true",
        default=False,
        help="When called, content in the existing cache dir is NOT deleted. By default, if the cache directory already exists, cazy_webscraper will delete the contents",
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
        help="Write out the log to a file with the provided name. If not provided, no log file is written",
    )  
    
    log_group.add_argument(
        "--log_dir",
        widget="DirChooser",
        metavar="Log directory",
        default=None,
        help="Define the directory to write out log files",
    )
    
    log_group.add_argument(
        "--nodelete_log",
        dest="nodelete_log",
        action="store_true",
        default=False,
        help="When called, content in the existing log dir is NOT deleted. By default, if the directory already exists, cazy_webscraper will delete the contents",
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

    misc_group = parser.add_argument_group(
        "Operation Options",
        "Additional operations to fine tune how cazy_webscraper operates"
    )

    misc_group.add_argument(
        "--batch_size",
        metavar="Batch size",
        widget="IntegerField",
        default=150,
        help="Size of batch queries sent to NCBI-Entrez"
    )

    misc_group.add_argument(
        "-r",
        "--retries",
        metavar="Connection retries",
        widget="IntegerField",
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
        help="Set SQLite engine echo to True. sqlalchemy (which is used to manage the local database) will print its log messages to the terminal",
    )

    misc_group.add_argument(
        "--cazy_synonyms",
        metavar="CAZy class synonyms",
        widget="FileChooser",
        default=None,
        help="Path to JSON file containing CAZy class synoymn (i.e. list of accepted alternative names for CAZy classes) to be used",
    )

    misc_group.add_argument(
        "-c",
        "--config",
        metavar="Configuration file",
        widget="FileChooser",
        default=None,
        help="Path to configuration file",
    )

    gooey_args = parser.parse_args()

    # compile path for the log file 
    if gooey_args.log is not None and gooey_args.log_dir is not None:
        gooey_args.log = Path(gooey_args.log_dir) / gooey_args.log
 
    get_genbank_sequences.main(args=gooey_args)


if __name__ == "__main__":
    main()
