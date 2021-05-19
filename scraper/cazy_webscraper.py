#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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

:cmd_args --config: path to configruration file
:cmd_args --classes: specify CAZy classes to scrape
:cmd_args --database: provide path to a local SQL database to add additional data to
:cmd_args --families: specify CAZy families to retrieve CAZymes from
:cmd_args --force: force overwriting content in exisiting output directory
:cmd_args --genera: specify Genera to retrieve CAZymes from
:cmd_args --kingdoms: specify taxonomy Kingdoms to scrape proteins from
:cmd_args --log: path to log file, enables writing out log messages to a log file
:cmd_args --nodelete: if true does not delete content in pre-existing output directory
:cmd_args --output: path to output directory
:cmd_args --retries: specify the number of times to try scraping a page if connection fails
:cmd_args --subfamilies: enable retrieval of subfamilies from CAZy
:cmd_args --species: specify species to retrieve CAZymes from
:cmd_args --strains: specify specific strains of species to retrieve CAZymes from
:cmd_args --timeout: specify the maximum time (in seconds) before determining connection timed out
:cmd_args --verbose: change logger level from warning to info, verbose logging
"""


import json
import logging
import os
import sys

import pandas as pd

from datetime import datetime
from pathlib import Path
from typing import List, Optional

from tqdm import tqdm

from scraper import crawler
from scraper.crawler.cazy_html_pages import get_cazy_pages, parse_local_pages
from scraper.crawler.parse_cazy_families import scrape_all, scrape_by_kingdom
from scraper.sql import sql_orm, sql_interface
from scraper.utilities import build_logger, config_logger, file_io, parsers, parse_configuration


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy.

    The collected data can be stored as a singel dataframe containing (not split), split into
    separate dataframes by class or by family. Excluded classes are CAZy classes not specified in
    the configuration file and thus, will not be scraped. User_cazy_families is the list of CAZy
    families specified to be scraped in the configration file.
    """
    # Program preparation
    time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    if argv is None:
        parser = parsers.build_parser()
        args = parser.parse_args()
    else:
        parser = parsers.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)

    if args.output is not sys.stdout:
        file_io.make_output_directory(args.output, args.force, args.nodelete)

    cazy_home = "http://www.cazy.org"

    # retrieve configuration data
    parse_configuration_path = parse_configuration.__file__
    (
        excluded_classes,
        config_dict,
        cazy_dict,
        taxonomy_filters_dict,
        kingdoms,
        ec_filters,
    ) = parse_configuration.parse_configuration(parse_configuration_path, args)
    # convert taxonomy_filters to a set for quicker identification of species to scrape
    taxonomy_filters = get_filter_set(taxonomy_filters_dict)

    # Check if retrieving pages from CAZy and writing to disk for scraping later
    if args.get_pages:
        get_cazy_pages.get_cazy_pages(
            args,
            cazy_home,
            time_stamp,
            excluded_classes,
            cazy_dict,
            config_dict,
            kingdoms,
            start_time,
        )

    else:
        # build database and return open database session
        if args.database is not None:  # open session for existing local database
            if args.database == "dict":  # build dictionary of {genbank_accession: CAZy families}
                session = {}

            else:
                if os.path.isfile(args.database) is False:
                    logger.error(
                        "Could not find local CAZy database. Check path is correct.\n"
                        "Had looked for a local SQL database at:"
                        f"{args.database}"
                        "Terminating programme."
                    )
                    sys.exit(1)
                try:
                    session = sql_orm.get_db_session(args)
                except Exception:
                    logger.error("Failed to build SQL database. Terminating program", exc_info=True)
                    sys.exit(1)

        else:  # create a new empty database to populate
            try:
                session = sql_orm.build_db(time_stamp, args)
            except Exception:
                logger.error("Failed to build SQL database. Terminating program", exc_info=True)
                sys.exit(1)

        if args.subfamilies is True:
            logger.warning("Enabled retrieval of subfamily classifications")

        if args.streamline is not None:
            parse_configuration.create_streamline_scraping_warning(args)

        # log scraping of CAZy in local db
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            taxonomy_filters_dict,
            kingdoms,
            ec_filters,
            session,
            args,
        )

        # Check if scraping from local CAZy files
        if args.scrape_files is not None:
            parse_local_pages.parse_local_pages(
                args,
                cazy_home,
                start_time,
                time_stamp,
                session,
                taxonomy_filters,
                ec_filters,
            )

        else:
            get_cazy_data(
                cazy_home,
                excluded_classes,
                config_dict,
                cazy_dict,
                taxonomy_filters,
                kingdoms,
                ec_filters,
                time_stamp,
                session,
                args,
            )

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    logger.info(
        "Finished scraping CAZy. Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}"
    )

    print(
        "=====================cazy_webscraper=====================\n"
        "Finished scraping CAZy\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}\n"
    )


def get_filter_set(taxonomy_filters_dict):
    """Create a set of all taxonomy filters from a dictionary.

    :param taxonomy_filers: dict of taxonomy filters

    Return a set.
    """
    taxonomy_filters = []

    for key in taxonomy_filters_dict:
        try:
            if len(taxonomy_filters_dict[key]) != 0:
                taxonomy_filters += taxonomy_filters_dict[key]
        except TypeError:
            pass

    if len(taxonomy_filters) == 0:
        taxonomy_filters = None

    else:
        taxonomy_filters = set(taxonomy_filters)

    return taxonomy_filters


def get_cazy_data(
    cazy_home,
    excluded_classes,
    config_dict,
    cazy_dict,
    taxonomy_filters,
    kingdoms,
    ec_filters,
    time_stamp,
    session,
    args,
):
    """Coordinate retrieval of data from the CAZy website.

    This function coordinates the crawling through the CAZy website by calling the appropriate
    functions, and then retrieving the protein data by calling to the appropriate data again.

    :param cazy_home: str, url of CAZy home page
    :param excluded_classes: list, list of classes to not scrape from CAZy
    :param config_dict: dict, user defined configuration of the scraper
    :param cazy_dict: dict, dictionary of excepct CAZy synonyms for CAZy classes
    :param taxonomy_filters: set of genera, species and strains to restrict the scrape to
    :param kingdoms: list of taxonomy kingdoms to restrict the scrape to
    :param ec_filters: set of EC numbers to limit the scrape to
    :param time_stamp: str, data and time scrape was initiated
    :param session: session, open database session
    :param args: cmd args parser

    Return nothing.
    """
    if args.output is not sys.stdout:
        out_log_path = args.output
    else:
        out_log_path = None

    connection_failures_logger = build_logger(
        out_log_path, f"CAZy_connection_failures_CW_{time_stamp}.log",
    )
    sql_failures_logger = build_logger(
        out_log_path, f"SQL_errors_CW_{time_stamp}.log",
    )
    format_failures_logger = build_logger(
        out_log_path, f"Format_and_parsing_errors_CW_{time_stamp}.log",
    )

    # retrieve links to CAZy class pages, return list of CazyClass objects
    cazy_classes = crawler.get_cazy_classes(
        cazy_home, excluded_classes, cazy_dict, args,
    )

    # scrape each retrieved class page
    for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):

        # first attempt of scraping, retrieve URLs to CAZy families
        if len(list(cazy_class.failed_families.keys())) == 0:

            # retrieve Family class instances, containing the Family page URL
            class_families, error_message, incorrect_urls = crawler.get_cazy_family_urls(
                cazy_class.url,
                cazy_class.name,
                cazy_home,
                args,
            )

            if incorrect_urls is not None:  # log families for which compiled URL is incorrect
                for url in incorrect_urls:
                    connection_failures_logger.warning(url)

            if class_families is None:  # couldn't retrieve URLs to families for working CAZy class
                cazy_class.tries += 1  # add one to the number of scrapping attempts

                # check if maximum number of attempts to connect have been met
                if cazy_class.tries == (args.retries + 1):
                    connection_failures_logger.warning(
                        f"{cazy_class.url}\t{cazy_class.name}\t"
                        f"No CAZy familes from this class were scraped\t{error_message}"
                    )

                else:
                    cazy_classes.append(cazy_class)  # retry scraping Class page later Classes

                continue

        # Not first try, scrape only the families for which a connections to CAZy previously failed
        else:
            class_families = list(cazy_class.failed_families.keys())

        # Scrape the familes of the current CAZy class, retrieving protein data

        if (config_dict is None) or (config_dict[cazy_class.name] is None):
            # No (sub)families were specified, therefore, scraping all families of the CAZy class

            for family in tqdm(class_families, desc=f"Parsing {cazy_class.name} families"):
                # Scrape the family and add proteins to the local database
                if kingdoms == 'all':
                    (
                        family,
                        retry_scrape,
                        failed_url_connections,
                        family_sql_failures,
                        format_errors,
                        session,
                    ) = scrape_all.parse_family_via_all_pages(
                        family,
                        cazy_home,
                        taxonomy_filters,
                        ec_filters,
                        args,
                        session,
                    )

                else:
                    (
                        family,
                        retry_scrape,
                        failed_url_connections,
                        family_sql_failures,
                        format_errors,
                        session,
                    ) = scrape_by_kingdom.parse_family_by_kingdom(
                        family,
                        cazy_home,
                        taxonomy_filters,
                        kingdoms,
                        ec_filters,
                        args,
                        session,
                    )

                # check if there are pages that were unsuccessfully scraped and have retries left
                if retry_scrape is True:
                    # add one to the number of attempted scrapes for CAZy family
                    try:
                        cazy_class.failed_families[family] += 1
                    except KeyError:
                        cazy_class.failed_families[family] = 1  # first attempt

                    # check if max number of attempts to connect family pages has been met
                    if cazy_class.failed_families[family] == (args.retries + 1):
                        del cazy_class.failed_families[family]  # do not try another scrape
                        continue

                if len(list(cazy_class.failed_families.keys())) != 0:
                    # if there are families with previously failed connection attempts
                    # and remaining tries, retry connection after working through other classes
                    cazy_classes.append(cazy_class)

                # write out errors to their respective log files
                for error in failed_url_connections:
                    connection_failures_logger.warning(error)

                for error in family_sql_failures:
                    sql_failures_logger.warning(error)

                for error in format_errors:
                    format_failures_logger.warning(error)

        else:
            # scrape only (sub)families specified in the config file

            for family in tqdm(class_families, desc=f"Parsing {cazy_class.name} families"):

                # Allow retrieval of subfamilies when only the parent CAZy family was named in the
                # config file, by searching by the family not subfamily in the config file
                if (args.subfamilies is True) and (family.name.find("_") != -1):
                    name_check = family.name[: (family.name.find("_"))]
                else:
                    name_check = family.name

                if name_check not in config_dict[cazy_class.name]:
                    continue

                # Scrape the family and add proteins to the local database
                if kingdoms == 'all':
                    (
                        family,
                        retry_scrape,
                        failed_url_connections,
                        family_sql_failures,
                        format_errors,
                        session,
                    ) = scrape_all.parse_family_via_all_pages(
                        family,
                        cazy_home,
                        taxonomy_filters,
                        ec_filters,
                        args,
                        session,
                    )

                else:
                    (
                        family,
                        retry_scrape,
                        failed_url_connections,
                        family_sql_failures,
                        format_errors,
                        session,
                    ) = scrape_by_kingdom.parse_family_by_kingdom(
                        family,
                        cazy_home,
                        taxonomy_filters,
                        kingdoms,
                        ec_filters,
                        args,
                        session,
                    )

                # check if there are pages that were unsuccessfully scraped and have retries left
                if retry_scrape is True:
                    # add one to the number of attempted scrapes for CAZy family
                    try:
                        cazy_class.failed_families[family] += 1
                    except KeyError:
                        cazy_class.failed_families[family] = 1  # first attempt

                    # check if max number of attempts to connect family pages has been met
                    if cazy_class.failed_families[family] == (args.retries + 1):
                        del cazy_class.failed_families[family]  # do not try another scrape
                        continue

                if len(list(cazy_class.failed_families.keys())) != 0:
                    # if there are families with previously failed connection attempts
                    # and remaining tries, retry connection after working through other classes
                    cazy_classes.append(cazy_class)

                # write out errors to their respective log files
                for error in failed_url_connections:
                    connection_failures_logger.warning(error)

                for error in family_sql_failures:
                    sql_failures_logger.warning(error)

                for error in format_errors:
                    format_failures_logger.warning(error)

    if type(session) is dict:
        if args.output is not sys.stdout:
            output_path = args.output / f"cazy_dict_{time_stamp}.json"
        else:
            output_path = Path(f"{os.getcwd()}/cazy_dict_{time_stamp}.json")
        with open(output_path, 'w') as f:
            json.dump(session, f)

    return


if __name__ == "__main__":
    main()
