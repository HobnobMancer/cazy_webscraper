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

from cazy_webscraper.expand.uniprot import get_uniprot_data
from cazy_webscraper.gui.assets import build_menus

cw_menu = build_menus(
    'cw_get_uniprot_data',
    'Retrieve protein data from UniProtKB afor CAZymes in a local CAZyme database. The retrieved data is stored in the local CAZyme database.'
)

@Gooey(
    program_name="Get Protein Data From UniProt Protein",
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
        help="The path to the local CAZyme database to extract protein sequences from",
    )

    # Add optional arguments to parser

    #
    # Data retrieval
    #

    retrieval_group = parser.add_argument_group(
        "Additional data options", 
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
            "Retrieve protein Aa sequences from UniProt"
        ),
    )

    update_group = parser.add_argument_group(
        "Update data options", 
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
            "sequence in the local database if a newer sequence is retrieved from UniProt"
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
            "protein name is retrieve for the same UniProt accession"
        ),
    )

    #
    # ACCESSIONS
    #

    accessions_group = parser.add_argument_group(
        "Use accessions", 
        "Provide a list  or lists of protein accessions to retrieve data for"
    )

    accessions_group.add_argument(
        "--genbank_accessions",
        metavar="GenBank accessions",
        widget="FileChooser",
        default=None,
        help="Path to text file containing GenBank accessions",
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
        "Options to find tune the operation of cazy_webscraper"
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
        help="Number of times to retry the connection to NCBI Entrez if the connection fails",
    )

    misc_group.add_argument(
        "-t",
        "--timeout",
        metavar="Timeout limit",
        widget="IntegerField",
        default=45,
        help="Connection timeout limit (seconds)",
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

    get_uniprot_data.main(args=gooey_args)


if __name__ == "__main__":
    main()
