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

:cmd_args --cache_dir: target path for cache
:cmd_args --cazy_synonms: path to yaml file containing CAZy class name synonms
:cmd_args --classes: specify CAZy classes to scrape
:cmd_args --citation: print citation information
:cmd_args --config: path to configruration file
:cmd_args --database: provide path to a local SQLite database to add additional data to
:cmd_args --db_output: path to write out new SQLite database
:cmd_args --families: specify CAZy families to retrieve CAZymes from
:cmd_args --force: force overwriting existing database
:cmd_args --genera: specify Genera to retrieve CAZymes from
:cmd_args --kingdoms: specify taxonomy Kingdoms to scrape proteins from
:cmd_args --log: path to log file, enables writing out log messages to a log file
:cmd_args --nodelete_cache: do not deleted existing content in cache dir
:cmd_args --nodelete_log: do not deleted existing content in log dir
:cmd_args --output: path to output directory
:cmd_args --retries: specify the number of times to try scraping a page if connection fails
:cmd_args --subfamilies: enable retrieval of subfamilies from CAZy
:cmd_args --species: specify species to retrieve CAZymes from
:cmd_args --strains: specify specific strains of species to retrieve CAZymes from
:cmd_args --timeout: specify the maximum time (in seconds) before determining connection timed out
:cmd_args --validate: retrieve CAZy fam population sizes and check against
:cmd_args --verbose: change logger level from warning to info, verbose logging
:cmd_args --version: print version info
"""


import logging
import json
import os

from datetime import datetime
from typing import List, Optional

import pandas as pd

from Bio import Entrez
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger, build_logger

from cazy_webscraper import (
    CITATION_INFO,
    VERSION_INFO,
    closing_message,
    connect_to_new_db,
    connect_existing_db,
    display_citation_info,
    display_version_info,
)
from cazy_webscraper.cache.cazy import cache_cazy_data
from cazy_webscraper.crawler.get_validation_data import get_validation_data
from cazy_webscraper.cazy import (
    build_taxa_dict,
    get_cazy_txt_file_data,
    parse_all_cazy_data,
    parse_cazy_data_with_filters,
)
from cazy_webscraper.ncbi.taxonomy.multiple_taxa import (
    identify_multiple_taxa,
    replace_multiple_tax,
)
from cazy_webscraper.sql import sql_orm, sql_interface
from cazy_webscraper.sql.sql_interface.add_data import add_cazyme_data
from cazy_webscraper.utilities import (
    parse_configuration,
    termcolour,
)
from cazy_webscraper.utilities.parsers.cazy_webscraper_parser import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy."""
    cazy_home_url = "http://www.cazy.org"

    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # Program preparation
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args, logger_name=__name__)

    # check if printing out version or citation information
    if args.version:
        display_version_info()
        return

    if args.citation:
        display_citation_info()
        return

    if args.email is None:
        logger.error(
            "No email address provided.\n"
            "Email address required by NCBI - which is required to retrieve the latest taxonomic\n"
            "classifications for proteins listed with multiple source organisms in the CAZy database\n"
            "Please provide your email address.\n"
            "Terminating program."
        )
        return

    # check correct output was provided, exit if not operable
    if args.database is not None and args.db_output is not None:
        warning_message = (
            "Target path for a NEW database (--db_output, -d) and\n"
            "a path to an EXISTING database (--database, -D) were provided.\n"
            "Please provide one OR the other.\n"
            "Terminating program."
        )
        logger.warning(termcolour(warning_message, "red"))
        closing_message("cazy_webscraper", start_time, args, early_term=True)
        return

    if args.db_output is not None and args.db_output.exists():
        if args.force:
            logger.warning(
                f"Local db {args.db_output} already exists\n"
                "Force is True\n"
                "Ovewriting existing database."
            )
            os.remove(args.db_output)
        else:
            logger.warning(
                f"Local db {args.db_output} already exists\n"
                "Force is False\n"
                "Not ovewriting existing database\n"
                "Termianting program"
            )
            closing_message("cazy_webscraper", start_time, args, early_term=True)
            return

    if args.skip_ncbi_tax:
        logger.warning(
            "skip_ncbi_tax is True - not retrieving the latest taxa from NCBI for proteins with multipe tax. Will use the first taxa listed in CAZy\n"
            "The latest taxonomic data can be retrieved using any of the three options:\n"
            "(i) cw_get_ncbi_taxs\n"
            "(ii) cw_get_genomics + cw_get_gtdb_taxs\n"
            "(iii) cw_get_uniprot_data with --taxonomy/-t"
        )

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

    scrape_config_message = (
        "Configuration:\n"
        f"Classes to scrape: {config_dict['classes']}\n"
        f"GH fams to scrape: {config_dict['Glycoside Hydrolases (GHs)']}\n"
        f"GT fams to scrape: {config_dict['GlycosylTransferases (GTs)']}\n"
        f"PL fams to scrape: {config_dict['Polysaccharide Lyases (PLs)']}\n"
        f"CE fams to scrape: {config_dict['Carbohydrate Esterases (CEs)']}\n"
        f"AA fams to scrape: {config_dict['Auxiliary Activities (AAs)']}\n"
        f"CBM fams to scrape: {config_dict['Carbohydrate-Binding Modules (CBMs)']}\n"
        f"Scraping subfamilies: {args.subfamilies}"
    )

    if len(taxonomy_filter_set) != 0:
        scrape_config_message += "\nTaxonomy filters applied."

    if len(kingdom_filters) < 5:
        scrape_config_message += f"\nScraping only tax kingdoms: {kingdom_filters}"

    logger.info(termcolour(scrape_config_message, "cyan"))

    if args.database is not None:  # adding data to an EXISTING database
        connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)
    else:  # build a new database
        connection, logger_name, cache_dir = connect_to_new_db(args, time_stamp, start_time)

    logger.info("Adding log of scrape to the local CAZyme database")
    with sql_orm.Session(bind=connection) as session:
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            kingdom_filters,
            taxonomy_filter_dict,
            set(),  # ec_filters not applied when scraping CAZy
            'CAZy',
            'CAZy annotations',
            session,
            args,
        )

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        make_output_directory(cache_dir, args.force, args.nodelete_cache)

    logger.warning(f"Created cache dir: {cache_dir}")

    if args.log is not None:  # write additional log files to user specified dir
        logger_name = args.log.name
        if logger_name.endswith(".log"):
            logger_name = logger_name[:-4]
        make_output_directory(args.log, args.force, args.nodelete_log)
    else:
        # write the additional log files to the .cazy_webscraper cache dire
        logger_name = "log"

    logger.info("Starting retrieval of data from CAZy")

    if args.cazy_data is not None:
        logger.warning(
            f"Retrieving CAZy data from predownloaded CAZy db dump at:\n{args.cazy_data}"
        )

    get_cazy_data(
        cazy_home_url,
        excluded_classes,
        cazy_class_synonym_dict,
        config_dict,
        class_filters,
        fam_filters,
        kingdom_filters,
        taxonomy_filter_set,
        connection,
        cache_dir,
        logger_name,
        time_stamp,
        args,
    )

    closing_message("cazy_webscraper", start_time, args)


def get_cazy_data(
    cazy_home_url,
    excluded_classes,
    cazy_class_synonym_dict,
    config_dict,
    class_filters,
    fam_filters,
    kingdom_filters,
    taxonomy_filters,
    connection,
    cache_dir,
    logger_name,
    time_stamp,
    args,
):
    """Coordinate retrieval of data from the CAZy website.

    This function coordinates the crawling through the CAZy website by calling the appropriate
    functions, and then retrieving the protein data by calling to the appropriate data again.

    :param cazy_home_url: str, url of CAZy home page
    :param excluded_classes: list, list of classes to not scrape from CAZy
    :param cazy_class_synonym_dict: dict of accepted CAZy class name synonyms
    :param config_dict: dict of CAZy families to scrape, or None if args.validate is False
    :param class_filters: set of CAZy classes to retrieve proteins from
    :param fam_filters: set of CAZy families to retrieve proteins from
    :param taxonomy_filters: set of genera, species and strains to restrict the scrape to
    :param connection: sqlalchemy connection obj, connection to SQLite db engine
    :param cache_dir: Path to dir to write out downloaded family txt files
    :param logger_name: str, name used for additional logger files
    :param time_stramp: str, time cazy_webscraper was invoked
    :param args: cmd args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # define paths for additional logs files
    # unless specifed they are added to the logs dir in the cache dir
    connection_failures_logger = build_logger(
        cache_dir, f"{logger_name}_{time_stamp}_connection_failures.log"
    )
    multiple_taxa_logger = build_logger(cache_dir, f"{logger_name}_{time_stamp}_multiple_taxa.log")
    replaced_taxa_logger = build_logger(cache_dir, f"{logger_name}_{time_stamp}_replaced_taxa.log")

    if args.validate:  # retrieve CAZy family population sizes for validating all data was retrieved
        # {fam (str): pop size (int)}
        cazy_fam_populations = get_validation_data(
            cazy_home_url,
            excluded_classes,
            cazy_class_synonym_dict,
            config_dict,
            cache_dir,
            connection_failures_logger,
            time_stamp,
            args,
        )
    else:
        cazy_fam_populations = None

    cazy_txt_lines = get_cazy_txt_file_data(cache_dir, time_stamp, args)

    logger.info(f"Retrieved {len(cazy_txt_lines)} lines from the CAZy db txt file")

    if (len(class_filters) == 0) and \
        (len(fam_filters) == 0) and \
            (len(kingdom_filters) == 0) and \
                (len(taxonomy_filters) == 0):
        cazy_data = parse_all_cazy_data(cazy_txt_lines, cazy_fam_populations, args)

    else:
        cazy_data = parse_cazy_data_with_filters(
            cazy_txt_lines,
            class_filters,
            fam_filters,
            kingdom_filters,
            taxonomy_filters,
            cazy_fam_populations,
            args,
        )

    logger.info(
        f"Retrieved {len((list(cazy_data.keys())))} proteins from the CAZy txt file "
        "matching the scraping criteria"
    )

    # check for GenBank accessions with multiple source organisms in the CAZy data
    multiple_taxa_gbks = identify_multiple_taxa(cazy_data, multiple_taxa_logger)

    if len(multiple_taxa_gbks) != 0:
        # remove the multiple taxa, and retrieve the latest taxa from NCBI
        cazy_data, successful_replacement = replace_multiple_tax(
            cazy_data,
            multiple_taxa_gbks,
            replaced_taxa_logger,
            args,
            invalid_ids=False,
        )

    # add separate kingdom and organism keys to cazy_data for all proteins
    taxa_dict, cazy_data = build_taxa_dict(cazy_data)  # {kingdom: {organisms}}

    # cache cazy_data dict
    cache_cazy_data(cazy_data, cache_dir)

    add_cazyme_data.add_kingdoms(taxa_dict, connection)

    add_cazyme_data.add_source_organisms(taxa_dict, connection)

    add_cazyme_data.add_cazy_families(cazy_data, connection)

    add_cazyme_data.add_genbanks(cazy_data, connection)

    add_cazyme_data.add_genbank_fam_relationships(cazy_data, connection, args)

    return


if __name__ == "__main__":
    main()
