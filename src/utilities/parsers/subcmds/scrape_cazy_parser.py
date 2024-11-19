#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
# (c) James Hutton Institute 2024
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# ehobbs@ebi.ac.uk
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


from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from src.entry_points import scrape_cazy

def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    parser = subps.add_parser(
        "scrape_cazy",
        description="Scrapes the CAZy database",
        help="Download data from CAZy and build a local SQLite database",
        formatter_class=ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "email",
        type=str,
        help=(
            "User email address.\n"
            "Required NCBI Entrez - used to get source organsism data.\n"
            "The email address is not stored be cazy_webscraper."
        )
    )

    filters_group = parser.add_argument_group("Filtering arguments")
    operational_group = parser.add_argument_group("Operational arguments")
    utilities_group = parser.add_argument_group("Utility arguments")

    filters_group.add_argument(
        "--classes",
        type=str,
        default=None,
        help="Classes from which all families are to be scraped. Separate classes by ','"
    )
    filters_group.add_argument(
        "--families",
        type=str,
        default=None,
        help="Families to scrape. Separate families by commas 'GH1,GH2' (case sensitive)"
    )
    filters_group.add_argument(
        "--kingdoms",
        type=str,
        default=None,
        help="Tax Kingdoms to restrict the scrape to"
    )
    filters_group.add_argument(
        "--genera",
        type=str,
        default=None,
        help="Genera to restrict the scrape to"
    )
    filters_group.add_argument(
        "--species",
        type=str,
        default=None,
        help="Species (written as Genus Species) to restrict the scrape to"
    )
    filters_group.add_argument(
        "--strains",
        type=str,
        default=None,
        help=(
            "Specific strains of organisms to restrict the scrape to "
            "(written as Genus Species Strain)"
        ),
    )

    operational_group.add_argument(
        "-o",
        "--db_output",
        type=Path,
        default=None,
        help="Build a NEW database. Provide the path to build new SQL database",
    )
    operational_group.add_argument(
        "-d",
        "--database",
        type=Path,
        default=None,
        help="Path to an EXISTING local CAZy SQL database. Add data to this database",
    )
    operational_group.add_argument(
        "-s",
        "--subfamilies",
        dest="subfamilies",
        action="store_true",
        default=False,
        help="Enable retrieval of subfamilies from CAZy",
    )
    operational_group.add_argument(
        "-c",
        "--config",
        type=Path,
        metavar="config file",
        default=None,
        help="Path to configuration file. Default: None, scrapes entire database",
    )
    operational_group.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force over writting an EXISTING database",
    )
    operational_group.add_argument(
        "--cazy_data",
        type=Path,
        default=None,
        help="Path predownloaded CAZy txt file",
    )
    operational_group.add_argument(
        "--delete_old_relationships",
        dest="delete_old_relationships",
        action="store_true",
        default=False,
        help=(
            "Delete old GenBank accession - CAZy family relationships (annotations)\n"
            "that are in the local db but are not in CAZy, e.g. when CAZy has moved a\n"
            "protein from one fam to another, delete the old family annotation."
        ),
    )
    operational_group.add_argument(
        "--skip_ncbi_tax",
        dest="skip_ncbi_tax",
        action="store_true",
        default=False,
        help=(
            "Skip retrieving the latest tax classification from the NCBI Taxonomy db for proteins\n"
            "listed with multiple taxs in CAZy.\n"
            "For these proteins the first taxonomy listed in CAZy is added to the local CAZyme db"
        ),
    )

    utilities_group.add_argument(
        "--cache_dir",
        type=Path,
        default=None,
        help="Target path for cache dir to be used instead of default path",
    )
    utilities_group.add_argument(
        "--cazy_synonyms",
        type=Path,
        default=None,
        help="Path to JSON file containing CAZy class synoymn names",
    )
    utilities_group.add_argument(
        "--ncbi_batch_size",
        type=int,
        default=200,
        help="Number of genbank accessions in each NCBI Taxonomy db batch query"
    )
    utilities_group.add_argument(
        "--nodelete_cache",
        dest="nodelete_cache",
        action="store_true",
        default=False,
        help="When called, content in the existing cache dir is NOT deleted",
    )
    utilities_group.add_argument(
        "-r",
        "--retries",
        type=int,
        default=10,
        help="Number of times to retry scraping a family or class page if error encountered",
    )
    utilities_group.add_argument(
        "--sql_verbose",
        dest="sql_echo",
        action="store_true",
        default=False,
        help="SQLite verbose logging",
    )
    utilities_group.add_argument(
        "-t",
        "--timeout",
        type=int,
        default=45,
        help="Connection timeout limit (seconds)"
    )

    parser.set_defaults(func=scrape_cazy.main)
