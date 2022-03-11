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
"""Create GUI for cazy_scraper.py."""

import argparse

from datetime import datetime
from pathlib import Path

from gooey import Gooey, GooeyParser

from cazy_webscraper.expand.extract_seqs import extract_db_seqs
from cazy_webscraper.gui.assets import build_menus

cw_menu = build_menus(
    'cw_extract_db_seqs',
    'Extract GenBank and/or UniProt protein sequences from the local CAZyme database, and write to a multiple sequence FASTA file, one sequence per FASTA file, and/or a BLAST database.'
)

@Gooey(
    program_name="Extract Db Sequences",
    image_dir="cazy_webscraper/gui/assets/extract_db_seqs",
    menu=cw_menu,
)
def main():
    # Create parser object
    parser = GooeyParser(
        prog="exract_db_seqs.py",
        description="Extract GenBank and/or UniProt protein sequences from the local CAZyme database, and write to a multiple sequence FASTA file, one sequence per FASTA file, and/or a BLAST database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "database",
        metavar="Local CAZyme database",
        widget="FileChooser",
        help="The path to the local CAZyme database to extract protein sequences from",
    )

    parser.add_argument(
        "source",
        metavar="Source of the sequences",
        nargs='+',
        choices=["GenBank", "UniProt", "GenBank and UniProt"],
        help="Original source(s) of the protein sequences",
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

    accessions_group.add_argument(
        "--uniprot_accessions",
        metavar="UniProt accessions",
        widget="FileChooser",
        default=None,
        help="Path to text file contining GenBank accessions",
    )

    #
    # OUTPUT OPTIONS
    #

    output_group = parser.add_argument_group(
        "Output Options", 
        "Configure where and what type of output is written"
    )

    output_group.add_argument(
        "-b",
        "--blastdb",
        metavar="Build BLAST database",
        widget="DirChooser",
        default=None,
        help=(
            "Create BLAST database of extracted protein sequences. "
            "Provide the path to the directory to store the database"
        ),
    )

    output_group.add_argument(
        "--fasta_dir",
        metavar="Individual FASTA file dir",
        widget="DirChooser",
        default=None,
        help="Write out each extracted sequence to a separate FASTA file in the specified dir",
    )

    output_group.add_argument(
        "--single_fasta_file",
        metavar="Multiple sequence FASTA file dir",
        widget="DirChooser",
        default=None,
        help="Write out all extracted sequence to a single, multiple sequence FASTA file in the specified dir",
    )

    output_group.add_argument(
        "--fasta_file",
        metavar="Name of the mutliple sequence FASTA file",
        type=str,
        default=None,
        help="If not provided a default name of the format 'cw_seqs_<date>_<time>.fasta' will be used",
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

    # compile multi-seq FASTA file path
    if gooey_args.single_fasta_file is not None:

        if gooey_args.fasta_file is None:
            time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
            fasta_path = f"cw_seqs_{time_stamp}.fasta"
            fasta_path = Path(gooey_args.single_fasta_file) / fasta_path
        else:
            fasta_path = Path(gooey_args.single_fasta_file) / gooey_args.fasta_file
        
        gooey_args.fasta_file = fasta_path

    # compile source of the protein sequences
    if gooey_args.source == ['GenBank and UniProt']:
        gooey_args.source = ['genbank', 'uniprot']
    elif gooey_args.source == ['GenBank']:
        gooey_args.source = ['genbank']
    else:
        gooey_args.source = ['uniprot']

    extract_db_seqs.main(args=gooey_args)


if __name__ == "__main__":
    main()
