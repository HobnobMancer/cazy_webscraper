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
"""Create GUI for get_pdb_structures.py."""


import argparse
import sys

from pathlib import Path

from gooey import Gooey, GooeyParser

from cazy_webscraper.expand.pdb import get_pdb_structures
from cazy_webscraper.gui import build_and_covert_to_paths
from cazy_webscraper.gui.assets import build_menus


cw_menu = build_menus(
    'cw_get_pdb_structures',
    'Retrieve protein structure files from the RCSB PDB database for PDB accessions in a local CAZyme database, and write the structure files to disk.'
)


@Gooey(
    program_name="Get Protein Structure Files",
    image_dir="cazy_webscraper/gui/assets/get_pdb_structures",
    menu=cw_menu,
)
def main():
    # Create parser object
    parser = GooeyParser(
        prog="get_pdb_structures.py",
        description="Retrieve protein structure files from PDB for CAZymes with PDB accession in a local CAZyme database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(  
        "database",
        metavar="Local CAZyme database",
        widget="FileChooser",
        help="The path to the local CAZyme database",
    )

    structure_group = parser.add_argument_group(
        "Structure files [Required]", 
        "Choose which structure file formats to download from PDB. At least one file type must be selected. More than one file type can be selected and the selected PDB accessions will be downloaded in all selected file formats"
    )

    structure_group.add_argument(
        "--mmcif",
        dest="mmcif",
        metavar="mmCif",
        action="store_true",
        default=False,
        help="Retrieve structure files in mmCif format",
    )
    structure_group.add_argument(
        "--pdb",
        dest="pdb",
        metavar="pdb",
        action="store_true",
        default=False,
        help="Retrieve structure files in pdb format",
    )
    structure_group.add_argument(
        "--xml",
        dest="xml",
        metavar="xml",
        action="store_true",
        default=False,
        help="Retrieve structure files in xml format",
    )
    structure_group.add_argument(
        "--bundle",
        dest="bundle",
        metavar="bundle",
        action="store_true",
        default=False,
        help="Retrieve structure files in bundle format",
    )


    # Add optional arguments to parser

    #
    # OUTPUT OPTIONS
    #

    output_group = parser.add_argument_group(
        "Output Options", 
        "Configure where the output is written"
    )

    output_group.add_argument(
        "-o",
        "--outdir",
        metavar="Output directory",
        widget="DirChooser",
        default=None,
        help=(
            "Directory to write build the new database in. "
            "If not specified the output will be written to the current working directory"
        ),
    )

    output_group.add_argument(
        "--overwrite",
        metavar="Overwrite existing files",
        dest="overwrite",
        action="store_true",
        default=False,
        help=(
            "Overwrite existing structure file with the same PDB accession as "
            "file being downloaded in the output directory. Default: don't overwrite existing file"
        ),
    )

    output_group.add_argument(
        "-f",
        "--force",
        metavar="Force",
        dest="force",
        action="store_true",
        default=False,
        help="Force writting to an existing output directory. If the output directory already exists and force is False, cazy_webscraper will not run",
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
        "--batch_size",
        metavar="Batch size",
        widget="IntegerField",
        default=150,
        help="Size of batch queries sent to RCSB PDB"
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

    # parse the PDB structure file choices
    if gooey_args.mmcif is False and gooey_args.pdb is False and gooey_args.xml is False and gooey_args.bundle is False:
        print("ERROR: At least one structure file format must be selected")
        sys.exit(1)

    gooey_args = build_and_covert_to_paths(gooey_args)

    if gooey_args.output_dir is not None:
        gooey_args.output_dir = Path(gooey_args.output_dir)

    if gooey_args.genbank_accessions is not None:
        gooey_args.genbank_accessions = Path(gooey_args.genbank_accessions)

    if gooey_args.uniprot_accessions is not None:
        gooey_args.uniprot_accessions = Path(gooey_args.uniprot_accessions)

    get_pdb_structures.main(args=gooey_args)


if __name__ == "__main__":
    main()
