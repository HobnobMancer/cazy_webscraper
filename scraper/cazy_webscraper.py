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
from scraper.crawler.parse_cazy_families import scrape_all, scrape_by_kingdom
from scraper.sql import sql_orm, sql_interface
from scraper.utilities import (
    build_logger,
    config_logger,
    file_io,
    parsers,
    parse_configuration,
    termcolour,
)


# Define constants

__version__ = "1.0.0-beta"
VERSION_INFO = [
    termcolour(
        f"cazy_webscraper version: {__version__}",
        "cyan",
    ),
]

CITATION_INFO = [
    termcolour(
        "If you use cazy_webscraper in your work, please cite the following publication:",
        "green",
    ),
    termcolour(
        "\tHobbs, E. E. M., Pritchard, L., Chapman, S., Gloster, T. M.,",
        "yellow",
    ),
    termcolour(
        "\t(2021) cazy_webscraper Microbiology Society Annual Conference 2021 poster. ",
        "yellow",
    ),
    termcolour(
        "\tFigShare. Poster.",
        "yellow",
    ),
    termcolour(
        "\thttps://doi.org/10.6084/m9.figshare.14370860.v7",
        "yellow",
    ),
]


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy.

    The collected data can be stored as a singel dataframe containing (not split), split into
    separate dataframes by class or by family. Excluded classes are CAZy classes not specified in
    the configuration file and thus, will not be scraped. User_cazy_families is the list of CAZy
    families specified to be scraped in the configration file.
    """
    # Program preparation
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
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

    # check if printing out version or citation information
    if args.version:
        sys.stderr.write("\n".join(VERSION_INFO) + "\n")
        return
    
    if args.citation:
        sys.stderr.write("\n".join(CITATION_INFO) + "\n")
        return

    if args.database and args.db_output:
        warning_message = (
            "Path to existing an existing and a target path to a new database were provided.\n"
            "Please provide one OR the other.\n"
            "Terminating program."
        )
        logger.warning(termcolour(warning_message, "red")
        )
    
    if args.no_db and (args.dict_output is None):
        warning_message = (
            "Opted to not write a database or dict/JSON file.\n"
            "Therefore, no method for writing out retrieved data was given.\n"
            "Please provide at least one method for writing out the retrieved data.\n"
            "Terminating program."
        )
        logger.warning(termcolour(warning_message, "red")
        )

    cazy_home_url = "http://www.cazy.org"

    logger.info("Parsing configuration")
    (
        excluded_classes,
        config_dict,
        cazy_class_synonym_dict,
        taxonomy_filter,
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
        f"Scraping subfamilies: {args.subfamilies}\n"
    )

    logger.info(termcolour(scrape_config_message, "cyan"))

    if args.database:
        logger.info("Adding data to an existing local CAZyme database")

        if os.path.isfile(args.database) is False:
            logger.error(
                "Could not find local CAZy database.\n"
                "Check path is correct.\n"
                "Terminating programme."
            )
            sys.exit(1)

        try:
            session = sql_orm.get_db_session(args)
            logger.info("Opened session to local CAZyme database")
        except Exception:
            logger.error("Failed to build SQL database. Terminating program\n", exc_info=True)
            sys.exit(1)
        
        logger.info("Adding log of scrape to the local CAZyme database")
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            taxonomy_filter,
            session,
            args,
        )
    
    else:
        if args.db_output is not None:
            logger.info("Building new local CAZyme database")
            logger.info(f"Output directory: {(args.output).parent}")
            logger.info(f"Force overwriting exiting output file: {args.force}")

            if str((args.db_output).parent) != '.':
                # dirs defined in output put
                file_io.make_db_output_directory(args)
            
            try:
                session = sql_orm.build_db(time_stamp, args)
                logger.info("Built new local CAZyme database")
            except Exception:
                logger.error("Failed to build SQL database. Terminating program", exc_info=True)
                sys.exit(1)
            
            logger.info("Adding log of scrape to the local CAZyme database")
            sql_interface.log_scrape_in_db(
                time_stamp,
                config_dict,
                taxonomy_filter,
                session,
                args,
            )
    
        else:
            if args.no_db:
                logger.warning("Selected to NOT generate a CAZyme database")
                session = None

            else:
                logger.warning(
                    "Default CAZyme database location selected:\n"
                    f"Database output in cwd: cazy_webscraper_{start_time}.db"
                )

                try:
                    session = sql_orm.build_db(f"cazy_webscraper_{start_time}.db", args)
                    logger.info("Built new local CAZyme database")
                except Exception:
                    logger.error("Failed to build SQL database. Terminating program", exc_info=True)
                    sys.exit(1)

                
    if args.dict_output is not None:
        logger.info("Building dict of CAZy family annotations")
        logger.info(f"Output directory: {(args.output).parent}")
        logger.info(f"Force overwriting exiting output file: {args.force}")

        if str((args.dict_output).parent) != '.':
            # dirs defined in output put
            file_io.make_cazy_dict_output_directory(args)
        
        cazy_dict = {}  # {protein_accession: [CAZy fams]}
    
    else:
        cazy_dict = None

    logger.info("Starting retrieval of data from CAZy")
    get_cazy_data(
        cazy_home_url,
        excluded_classes,
        config_dict,
        cazy_class_synonym_dict,
        taxonomy_filter,
        session,
        cazy_dict,
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
        f"Version: {VERSION_INFO}\n"
        f"Citation: {CITATION_INFO}"
    )

    print(
        "=====================cazy_webscraper=====================\n"
        "Finished scraping CAZy\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}"
        f"Version: {VERSION_INFO}\n"
        f"Citation: {CITATION_INFO}"
    )


def get_cazy_data(
    cazy_home_url,
    excluded_classes,
    config_dict,
    cazy_class_synonym_dict,
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

    :param cazy_home_url: str, url of CAZy home page
    :param excluded_classes: list, list of classes to not scrape from CAZy
    :param config_dict: dict, user defined configuration of the scraper
    :param cazy_class_synonym_dict: dict, dictionary of excepct CAZy synonyms for CAZy classes
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
        cazy_home_url, excluded_classes, cazy_class_synonym_dict, args,
    )

    # scrape each retrieved class page
    for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):

        # first attempt of scraping, retrieve URLs to CAZy families
        if len(list(cazy_class.failed_families.keys())) == 0:

            # retrieve Family class instances, containing the Family page URL
            class_families, error_message, incorrect_urls = crawler.get_cazy_family_urls(
                cazy_class.url,
                cazy_class.name,
                cazy_home_url,
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
                        cazy_home_url,
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
                        cazy_home_url,
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
                        cazy_home_url,
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
                        cazy_home_url,
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
            output_path = args.output / f"cazy_class_synonym_dict_{time_stamp}.json"
        else:
            output_path = Path(f"{os.getcwd()}/cazy_class_synonym_dict_{time_stamp}.json")
        with open(output_path, 'w') as f:
            json.dump(session, f)

    return


if __name__ == "__main__":
    main()
