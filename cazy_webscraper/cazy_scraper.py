#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
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
"""
Web scraper to scrape CAZy website and retrieve all protein data.
"""


import argparse
import logging
import sys

from datetime import datetime
from typing import List, Optional

import pandas as pd

from pathlib import Path

from Bio import Entrez
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger

from cazy_webscraper import (
    closing_message,
    display_citation_info,
    display_version_info,
)
from cazy_webscraper.cazy.download import get_cazy_db_dump
from cazy_webscraper.cazy.filter_data import (
    apply_tax_filters,
    apply_class_and_family_filters,
    drop_subfamilies
)
from cazy_webscraper.cazy.multi_taxa import process_multi_taxa

from cazy_webscraper.database.connect import (
    connect_to_new_db,
    connect_existing_db,
)

from cazy_webscraper.database.scrape_log import add_main_scrape_message
from cazy_webscraper.database.cazy import dump_cazy_txt

from cazy_webscraper.sql.sql_interface.add_data import add_cazy_data
from cazy_webscraper.utilities import (
    parse_configuration,
    sanity_checks
)
from cazy_webscraper.utilities.parsers.cazy_webscraper_parser import build_parser


logger = logging.getLogger(__name__)


def main(argv: Optional[List[str]] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy."""
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    if not argv:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    config_logger(args, logger_name=__name__)

    if args.version:
        display_version_info()
        return

    if args.citation:
        display_citation_info()
        return

    db = sanity_checks.sanity_check_main_input(time_stamp, args)  # path to the local cazyme db

    if args.skip_ncbi_tax:
        logger.warning(
            "skip_ncbi_tax is True\n"
            "The latest taxa from NCBI for proteins with multipe tax in CAZy will not be retrieved\n."
            "The first taxonomy retrieved from CAZy will be used instead.\n"
            "The latest taxonomic data can be retrieved later using any of the three options:\n"
            "(i) cw_get_ncbi_taxs\n"
            "(ii) cw_get_genomics + cw_get_gtdb_taxs\n"
            "(iii) cw_get_uniprot_data with --taxonomy/-t"
        )
    else:
        Entrez.email = args.email

    logger.info("Parsing configuration")
    (
        excluded_classes,
        config_dict,
        cazy_class_synonym_dict,
        class_filters,
        fam_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        taxonomy_filter_set,
    ) = parse_configuration.parse_configuration(args)

    if args.database:
        connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)
    else:
        connection, logger_name, cache_dir = connect_to_new_db(args, time_stamp, start_time)

    if args.cache_dir:  # use user defined cache dir
        cache_dir = args.cache_dir
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    logger.warning("Created cache dir: %s", cache_dir)

    if args.log:  # write additional log files to user specified dir
        logger_name = args.log.name
        if logger_name.endswith(".log"):
            logger_name = logger_name[:-4]
        make_output_directory(args.log, args.force, args.nodelete_log)
    else:
        # write the additional log files to the .cazy_webscraper cache dire
        logger_name = "log"

    add_main_scrape_message(
        kingdom_filters,
        taxonomy_filter_set,
        taxonomy_filter_dict,
        time_stamp,
        config_dict,
        args,
        connection
    )

    logger.info("Starting retrieval of data from CAZy")

    if args.cazy_data:
        logger.warning(
            "Retrieving CAZy data from predownloaded CAZy db dump at:\n%s", args.cazy_data
        )

    get_cazy_data(
        class_filters,
        fam_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        cache_dir,
        logger_name,
        time_stamp,
        args,
        db,
    )

    closing_message("cazy_webscraper", start_time, args)


def get_cazy_data(
    class_filters: set[str],
    fam_filters: set[str],
    kingdom_filters: set[str],
    taxonomy_filter_dict: dict,
    cache_dir: Path,
    logger_name: str,
    time_stamp: str,
    args: argparse.Namespace,
    db: Path,
):
    """Coordinate retrieval of data from the CAZy website.

    This function coordinates the crawling through the CAZy website by calling the appropriate
    functions, and then retrieving the protein data by calling to the appropriate data again.

    :param class_filters: set of CAZy classes to retrieve proteins from
    :param fam_filters: set of CAZy families to retrieve proteins from
    :param kingdom_filters: set of kingdoms to keep
    :param taxonomy_filters: set of genera, species and strains to restrict the scrape to
    :param cache_dir: Path to dir to write out downloaded family txt files
    :param logger_name: str, name used for additional logger files
    :param time_stramp: str, time cazy_webscraper was invoked
    :param args: cmd args parser
    :param db: path to local cazyme db

    Return nothing.
    """
    if args.cazy_data:
        dump_cazy_txt(args.cazy_data, db)
    else:
        cazy_txt_path = get_cazy_db_dump(cache_dir, time_stamp, args)
        dump_cazy_txt(cazy_txt_path, db)

    # filter data in the cazy db dump to only retain records that match the user criteria
    if any([
        kingdom_filters,
        taxonomy_filter_dict['genera'],
        taxonomy_filter_dict['species'],
        taxonomy_filter_dict['strains'],
    ]):
        logger.warning("Applying taxonomic filters")
        apply_tax_filters(
            kingdom_filters,
            taxonomy_filter_dict['genera'],
            taxonomy_filter_dict['species'],
            taxonomy_filter_dict['strains'],
            db
        )

    if any([class_filters, fam_filters]):
        logger.warning("Applying class and family filters")
        apply_class_and_family_filters(class_filters, fam_filters, db)

    if not args.subfamilies:
        logger.warning("Dropping subfamilies")
        drop_subfamilies(db)

    # deal with instances of multiple taxonomies
    process_multi_taxa(
        db,
        cache_dir / f"{logger_name}_{time_stamp}_multi_taxa.log",
        cache_dir / f"{logger_name}_{time_stamp}_replaced_taxa.log",
        args,
    )

    # add data to the database
    add_cazy_data.add_kingdoms(db)
    add_cazy_data.add_source_organisms(db)
    sys.exit(1)

    add_cazy_data.add_cazy_families(cazy_data, connection)

    add_cazy_data.add_genbanks(cazy_data, connection)

    add_cazy_data.add_genbank_fam_relationships(cazy_data, connection, args)

    # delete the TempTable
    ######


if __name__ == "__main__":
    main()
