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
"""Create GUI for cw_query_database.py."""

import argparse

from datetime import datetime
from pathlib import Path

from gooey import Gooey, GooeyParser

from cazy_webscraper.api import cw_query_database
from cazy_webscraper.gui import build_and_covert_to_paths
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
        help="The path to the local CAZyme database",
    )

    parser.add_argument(
        "file_types",
        metavar="Output file formats",
        action="store",
        choices=["csv", "json", "csv and json"],
        help="File types of output files"
    )

    #
    # Output options
    #

    output_group = parser.add_argument_group(
        "Output options", 
        "By default the the output files are written to the current working directory. Detfine a different output directory and/or add a prefix to all output file names",
    )

    output_group.add_argument(
        "-o",
        "--output_dir",
        widget="DirChooser",
        default=None,
        help="Path to output dir, default: None (writes to cwd)",
    )

    output_group.add_argument(
        "-p",
        "--prefix",
        type=str,
        default=None,
        help="String to prefix all output files with, default: None",
    )

    output_group.add_argument(
        "--overwrite",
        dest="overwrite",
        action="store_true",
        default=False,
        help="Overwrite existing output files",
    )

    output_group.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force writing to existing output directory",
    )

    output_group.add_argument(
        "-n",
        "--nodelete",
        metavar="Do NOT delete content in output directory",
        dest="nodelete",
        action="store_true",
        default=False,
        help="When called, content in the existing out dir is NOT deleted. By default cazy_webscraper deletes content in the existing output dir",
    )

    #
    # Data retrieval
    #

    retrieval_group = parser.add_argument_group(
        "Output Data", 
        "By default the GenBank accessions of proteins matching the provided criteria are included in the output. Choose additional data can be included in the output files"
    )

    retrieval_group.add_argument(
        "--include",
        metavar="CAZy class",
        dest="include",
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
        dest="include_gbk_seq",
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
        "--include_ec",
        dest="include_ec",
        metavar="EC numbers",
        action="store_true",
        default=False,
        help="Retrieve EC numbers from UniProt",
    )

    retrieval_group.add_argument(
        "--include_pdb",
        metavar="PDB accessions",
        dest="include_pdb",
        action="store_true",
        default=False,
        help="Retrieve PDB accessions from UniProt",
    )

    retrieval_group.add_argument(
        "--include_uni_seq",
        metavar="UniProt protein sequences",
        dest="include_uni_seq",
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
        "Provide a list or lists of protein accessions to explicitly define the proteins to retrieve data for. The provided protein accessions will be used instead of the filters to define proteins of interest"
    )

    accessions_group.add_argument(
        "--genbank_accessions",
        metavar="GenBank accessions",
        widget="FileChooser",
        default=None,
        help="Path to text file contining GenBank accessions of the proteins of interset",
    )

    accessions_group.add_argument(
        "--uniprot_accessions",
        metavar="UniProt accessions",
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
        help=("Path to JSON file containing CAZy class synoymn (i.e. list of accepted alternative names for CAZy classes) to be used"),
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

    # Parse the selected file types into a single args
    if gooey_args.file_types == ["csv and json"]:
        gooey_args.file_types = ["csv", "json"]

    # Parse the data to be included in the output into a single args
    included_data = []
    if gooey_args.include:
        included_data.append("class")
    if gooey_args.include_fam:
        included_data.append("family")
    if gooey_args.include_subfam:
        included_data.append("subfamily")  
    if gooey_args.include_genus:
        included_data.append("genus")
    if gooey_args.include_organism:
        included_data.append("organism")
    if gooey_args.include_uni_acc:
        included_data.append("uniprot_acc")
    if gooey_args.include_name:
        included_data.append("uniprot_name")
    if gooey_args.include_ec:
        included_data.append("ec")
    if gooey_args.include_pdb:
        included_data.append("pdb")
    if gooey_args.include_uni_seq:
        included_data.append("uniprot_seq")
    if gooey_args.include_gbk_seq:
        included_data.append("genbank_seq")
    
    gooey_args.include = included_data
    
    gooey_args = build_and_covert_to_paths(gooey_args)

    cw_query_database.main(args=gooey_args)


if __name__ == "__main__":
    main()
