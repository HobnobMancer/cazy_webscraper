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

from cazy_webscraper import cazy_scraper
from cazy_webscraper.gui.assets import build_menus

cw_menu = build_menus(
    'cazy_webscraper',
    'A tool to retrieve data from CAZy and compile a local CAZyme database.'
)

@Gooey(
    program_name="cazy_webscraper",
    image_dir="cazy_webscraper/gui/assets/cazy_webscraper",
    menu=cw_menu,
)
def main():
    # Create parser object
    parser = GooeyParser(
        prog="cazy_scraper.py",
        description="Scrape the CAZy database and compile proteins matching specific criteria into a local CAZyme database",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "email",
        metavar="Email address",
        type=str,
        help=(
            "This is a requirement of NCBI Entrez, which is used to get source organsism data. The email is not stored be cazy_webscraper."
        ),
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
        "--db_output",
        metavar="Database directory",
        widget="DirChooser",
        default=None,
        help="Directory to write build the new database in",
    )

    output_group.add_argument(
        "--new_database_name",
        metavar="New database name",
        type=str,
        default=None,
        help="Name of the new database file. If not given, the default name is used.",
    )

    output_group.add_argument(
        "-d",
        "--database",
        metavar="Existing database",
        widget="FileChooser",
        default=None,
        help="An existing local CAZy database to add data to",
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
            "Classes from which all families are to be scraped.\n"
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

    class_group.add_argument(
        "-s",
        "--subfamilies",
        metavar="Retrieve subfamilies",
        dest="subfamilies",
        action="store_true",
        default=False,
        help="Enable retrieval of subfamilies from CAZy",
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

    # Add option to use a pre-downloaded CAZy txt file
    cache_group.add_argument(
        "--cazy_data",
        metavar="CAZy db dump",
        widget="FileChooser",
        default=None,
        help="Path to predownloaded CAZy txt file. Use data from a previously downloaded CAZy dump txt file",
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
        "--nodelete_log",
        dest="nodelete_log",
        action="store_true",
        default=False,
        help="When called, content in the existing log dir is NOT deleted",
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
        "-r",
        "--retries",
        metavar="Number of retries",
        widget="IntegerField",
        default=10,
        help="Number of times to retry scraping a family or class page if error encountered",
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

    # Add option to define time out limit for trying to connect to CAZy
    misc_group.add_argument(
        "-t",
        "--timeout",
        metavar="Timeout limit",
        widget="IntegerField",
        default=45,
        help="Connection timeout limit (seconds)"
    )

    misc_group.add_argument(
        "--validate",
        metavar="Validate data (primitive)",
        dest="validate",
        action="store_true",
        default=False,
        help=(
            "Retrieve CAZy fam population sizes from CAZy and use to check\n"
            "the number of family members added to the local database"
        ),
    )

    misc_group.add_argument(
        "--cazy_synonyms",
        metavar="CAZy class synonyms",
        widget="FileChooser",
        default=None,
        help=(
            "Path to JSON file containing CAZy class synoymn names\n"
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

    misc_group.add_argument(
        "--delete_old_relationships",
        metavar="Delete old family annotations",
        dest="delete_old_relationships",
        action="store_true",
        default=False,
        help=(
            "Delete old GenBank accession - CAZy family relationships (annotations)\n"
            "that are in the local db but are not in CAZy, e.g. when CAZy has moved a\n"
            "protein from one fam to another, delete the old family annotation."
        ),
    )

    #
    # CITATION AND VERSION
    #

    cit_ver_group = parser.add_argument_group(
        "Citation and version data",
        "Print the citation and/or version. CAZy will not be scraped.",
    )

    cit_ver_group.add_argument(
        "-C",
        "--citation",
        metavar="Citation",
        dest="citation",
        action="store_true",
        default=False,
        help="Print cazy_webscraper citation message",
    )

    cit_ver_group.add_argument(
        "-V",
        "--version",
        metavar="Version",
        dest="version",
        action="store_true",
        default=False,
        help="Print cazy_webscraper version number",
    )    

    gooey_args = parser.parse_args()

    # compile db_output path
    if gooey_args.db_output is not None:

        if gooey_args.new_database_name is None:
            time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
            db_path = f"cazy_webscraper_{time_stamp}.db"
            db_path = Path(gooey_args.db_output) / db_path
        else:
            db_path = Path(gooey_args.db_output) / gooey_args.new_database_name
        
        gooey_args.db_output = db_path

    cazy_scraper.main(args=gooey_args)


if __name__ == "__main__":
    main()
