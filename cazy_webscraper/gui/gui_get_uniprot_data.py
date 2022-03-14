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
"""Create GUI for get_uniprot_data.py."""

import argparse

from datetime import datetime
from pathlib import Path

from gooey import Gooey, GooeyParser

from cazy_webscraper.expand.uniprot import get_uniprot_data
from cazy_webscraper.gui import build_and_covert_to_paths
from cazy_webscraper.gui.assets import build_menus


cw_menu = build_menus(
    'cw_get_uniprot_data',
    'Retrieve protein data from UniProtKB for CAZymes in a local CAZyme database. The retrieved data is stored in the local CAZyme database.'
)


@Gooey(
    program_name="Get Protein Data From UniProt",
    image_dir="cazy_webscraper/gui/assets/get_uniprot_data",
    menu=cw_menu,
)
def main():
    # Create parser object
    parser = GooeyParser(
        prog="get_uniprot_data.py",
        description="Retrieve protein data from UniProt and store the data in the local CAZyme database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "database",
        metavar="Local CAZyme database",
        widget="FileChooser",
        help="The path to the local CAZyme database",
    )

    # Add optional arguments to parser

    #
    # Data retrieval
    #

    retrieval_group = parser.add_argument_group(
        "Additional Data Options", 
        "Choose additional data to retrieve from UniProt, in addition to the UniProt ID and protein name"
    )

    retrieval_group.add_argument(
        "-e",
        "--ec",
        dest="ec",
        metavar="EC numbers",
        action="store_true",
        default=False,
        help="Retrieve EC numbers from UniProt",
    )

    retrieval_group.add_argument(
        "-p",
        "--pdb",
        metavar="PDB accessions",
        dest="pdb",
        action="store_true",
        default=False,
        help="Retrieve PDB accessions from UniProt",
    )

    retrieval_group.add_argument(
        "-s",
        "--sequence",
        metavar="Protein sequences",
        dest="sequence",
        action="store_true",
        default=False,
        help=(
            "Retrieve protein sequences from UniProt"
        ),
    )

    update_group = parser.add_argument_group(
        "Update Data Options", 
        "Overwrite existing data in the local database if a more recent version of the data is retrieved from UniProt"
    )

    update_group.add_argument(
        "--seq_update",
        dest="seq_update",
        metavar="Update sequence",
        action="store_true",
        default=False,
        help=(
            "When retrieving protein sequences from UniProt overwrite the existing "
            "sequence in the local database if a more recent version of the sequence is retrieved from UniProt"
        ),
    )

    update_group.add_argument(
        "--name_update",
        dest="name_update",
        metavar="Update protein name",
        action="store_true",
        default=False,
        help=(
            "Overwrite the existing protein name in the local database if a different"
            "protein name is retrieve for the same protein from UniProt"
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
            "Classes to retrieve protein data for. "
            "Separate classes with a single comma ',', e.g. 'GH,GT,PL'"
        ),
    )

    class_group.add_argument(
        "--families",
        metavar="CAZy (sub)families",
        type=str,
        default=None,
        help="Families and subfamilies to retrieve protein data for. Separate families with single commas 'GH1,GH3_1'. CAZy families are case sensitive"
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
        "--skip_uniprot_accessions",
        metavar="Use cached UniProt accessions",
        widget="FileChooser",
        default=None,
        help="Path to a JSON file containing the cached UniProt IDs, GenBank accesions and database IDs. (This skips the retrieval of UniProt IDs from UniProt)",
    )

    cache_group.add_argument(
        "--use_uniprot_cache",
        metavar="Use cached UniProt data",
        widget="FileChooser",
        default=None,
        help="Path to a JSON file containing data previously retrieved from UniProt. (This skips retrieving data from UniProt)",
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
    # DELETE DATA
    #

    delete_group = parser.add_argument_group(
        "Delete Old Data",
        "Delete data that is stored in the local database but no longer stored in UniProt. This permantly removes the data from the local database"
    )

    delete_group.add_argument(
        "--delete_old_ec",
        dest="delete_old_ec",
        metavar="Delete old EC annotations",
        action="store_true",
        default=False,
        help=(
            "Delete EC number-GenBank relationships/annotations "
            "that are listed in the local database but are not longer "
            "listed in UniProt for a given protein. This will not delete the "
            "EC number from the local database."
        ),
    )

    delete_group.add_argument(
        "--delete_old_pdbs",
        dest="delete_old_pdbs",
        metavar="Delete old PDBs accessions",
        action="store_true",
        default=False,
        help=(
            "Delete PDB accessions that are stored in the local database "
            "but are no longer stored in UniProt. This permamently removes "
            "the PDB accessions from the local database."
        ),
    )


    #
    # MISC OPTIONS
    #

    misc_group = parser.add_argument_group(
        "Operation Options",
        "Additional operations to fine tune how cazy_webscraper operates"
    )

    misc_group.add_argument(
        "--bioservices_batch_size",
        metavar="Bioservices batch size",
        widget="IntegerField",
        default=150,
        help=(
            "Batch size for queries parsed by bioservices to retrieve "
            "protein data from UniProt"
        ),
    )

    misc_group.add_argument(
        "--uniprot_batch_size",
        metavar="UniProt batch size",
        widget="IntegerField",
        default=150,
        help=(
            "Batch size for queries sent to the UniProt REST API to "
            "retrieve UniProt IDs for a list of GenBank accessions"
        ),
    )

    misc_group.add_argument(
        "-r",
        "--retries",
        metavar="Connection retries",
        widget="IntegerField",
        default=10,
        help="Number of times to retry the connection to UniProt if the connection fails",
    )

    misc_group.add_argument(
        "-t",
        "--timeout",
        metavar="Timeout limit",
        widget="IntegerField",
        default=45,
        help="Connection timeout limit (seconds)",
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

    gooey_args = build_and_covert_to_paths(gooey_args)

    if gooey_args.genbank_accessions is not None:
        gooey_args.genbank_accession = Path(gooey_args.genbank_accessions)

    if gooey_args.skip_uniprot_accessions is not None:
        gooey_args.skip_uniprot_accessions = Path(gooey_args.skip_uniprot_accessions)
    
    if gooey_args.use_uniprot_cache is not None:
        gooey_args.use_uniprot_cache = Path(gooey_args.use_uniprot_cache)

    get_uniprot_data.main(args=gooey_args)


if __name__ == "__main__":
    main()
