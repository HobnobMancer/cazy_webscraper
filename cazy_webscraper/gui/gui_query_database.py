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

from cazy_webscraper.api import cw_query_database
from cazy_webscraper.gui.assets import build_menus

cw_menu = build_menus(
    'cw_query_database',
    'Interrogate and retrieve data from a local CAZyme database.'
)

@Gooey(
    program_name="Query The Local Database",
    image_dir="cazy_webscraper/gui/assets/query_database",
    menu=cw_menu,
)
def main():
    # Create parser object
    parser = GooeyParser(
        prog="cw_query_database.py",
        description="Retrieve data from local database matching specific criteria",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "database",
        metavar="Local CAZyme database",
        widget="FileChooser",
        help="The path to the local CAZyme database to extract protein sequences from",
    )

    parser.add_argument(
        "file_types",
        metavar="Output file formats",
        action="store",
        choices=["csv", "json", "csv and json"],
        help="File types to write the query output in [csv,json]"
    )

    #
    # Data retrieval
    #

    retrieval_group = parser.add_argument_group(
        "Output Data", 
        "By default the GenBank accessions of proteins matching the provided criteria are included in the output. Choose additional data can be included in the output files"
    )

    retrieval_group.add_argument(
        "--include_class",
        metavar="CAZy class",
        dest="include_class",
        action="store_true",
        default=False,
        help="Include CAZy class annotations",
    )

    retrieval_group.add_argument(
        "--include_fam",
        metavar="CAZy family",
        dest="include_fam",
        action="store_true",
        default=False,
        help="Include CAZy family annotations",
    )

    retrieval_group.add_argument(
        "--include_subfam",
        metavar="CAZy subfamily",
        dest="include_subfam",
        action="store_true",
        default=False,
        help="Include CAZy subfamily annotations",
    )

    retrieval_group.add_argument(
        "--include_genus",
        metavar="Genus",
        dest="include_genus",
        action="store_true",
        default=False,
        help="Include the genus of the source organism",
    )

    retrieval_group.add_argument(
        "--include_organism",
        metavar="Organism",
        dest="include_organism",
        action="store_true",
        default=False,
        help="Include the scientific name of the source organism",
    )

    retrieval_group.add_argument(
        "--include_gbk_seq",
        metavar="GenBank protein sequences",
        dest="sequence",
        action="store_true",
        default=False,
        help="Including protein sequences retrieved from GenBank",
    )

    retrieval_group.add_argument(
        "--include_uni_acc",
        metavar="UniProt accession",
        dest="include_uni_acc",
        action="store_true",
        default=False,
        help="Include UniProt accessions",
    )

    retrieval_group.add_argument(
        "--include_name",
        metavar="Protein name",
        dest="include_name",
        action="store_true",
        default=False,
        help="Include the protein name retrieved from UniProt",
    )

    retrieval_group.add_argument(
        "--ec",
        dest="ec",
        metavar="EC numbers",
        action="store_true",
        default=False,
        help="Retrieve EC numbers from UniProt",
    )

    retrieval_group.add_argument(
        "--pdb",
        metavar="PDB accessions",
        dest="pdb",
        action="store_true",
        default=False,
        help="Retrieve PDB accessions from UniProt",
    )

    retrieval_group.add_argument(
        "--uniprot_sequence",
        metavar="UniProt protein sequences",
        dest="sequence",
        action="store_true",
        default=False,
        help=(
            "Including protein sequences retrieved from UniProt"
        ),
    )

    #
    # ACCESSIONS
    #

    accessions_group = parser.add_argument_group(
        "Use Accessions", 
        "Provide a list  or lists of protein accessions to retrieve data for. This overrides providing filters to identify proteins matching the specified criteria"
    )

    accessions_group.add_argument(
        "--genbank_accessions",
        metavar="GenBank accessions",
        widget="FileChooser",
        default=None,
        help="Path to text file containing GenBank accessions",
    )

    accessions_group.add_argument(
        "--uniprot_accessions",
        metavar="UniProt accessions",
        widget="FileChooser",
        default=None,
        help="Path to text file containing UniProt accessions",
    )

    #
    # CLASS AND FAM FILTERS
    #

    class_group = parser.add_argument_group(
        "Class and Family Filters", 
        "Define CAZy classes, families and subfamilies to retrieve data from"
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
        "Taxonomy Filters",
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
        "Additional Filters", 
        "Further refine the protein sequences retrieved from the local database"
    )

    add_filt_group.add_argument(
        "--ec_filter",
        metavar="EC number filters",
        type=str,
        default=None,
        help="Retrieve data for proteins annotated with at least one of the provided EC numbers. Separate EC numbers with a single comma, e.g. 1.2.3.4,2.3.4.5. The 'EC' prefix is optional, and '*' and *-* are accepted for missing digits"
    )

    #
    # CACHE OPTIONS
    #

    cache_group = parser.add_argument_group(
        "Cache Options", 
        "Use cache files and change the cache location"
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

    misc_group = parser.add_argument_group(
        "Operation Options",
        "Options to find tune the operation of cazy_webscraper"
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

    # Parse the selected file types into a single args
    if gooey_args.file_types == ["csv and json"]:
        gooey_args.file_types = ["csv", "json"]

    # Parse the data to be included in the output into a single args


    cw_query_database.main(args=gooey_args)


if __name__ == "__main__":
    main()
