#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
"""
Web scraper to scrape CAZy website and retrieve all protein data.

:cmd_args --config: path to configruration file
:cmd args --data_split: [None, class, family] how data is to be split/separated into dataframes
:cmd_args --force: force overwriting content in exisiting output directory
:cmd_args --log: path to log file, enables writing out log messages to a log file
:cmd_args --nodelete: if true does not delete content in pre-existing output directory
:cmd_args --output: path to output directory
:cmd_args --subfamily: enable retrieval of subfamilies from CAZy
:cmd_args --verbose: change logger level from warning to info, verbose logging

:func main: coordinate scraping of CAZy database
:func get_class_urls: retrieve URLs to CAZy class summary pages from CAZy homepage
:func browser_decorator: decorator for get_page() to coordinate retrying failed connections
:func get_page: connect to webpage and retrieval page as BeautifulSoup4 object

:class Protein: A single protein from CAZy database
:class Family: A single family from CAZy containing proteins
"""

import logging
import re
import sys

from datetime import datetime
from typing import List, Optional

from tqdm import tqdm

from scraper import crawler, file_io, parse, sql, utilities


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy.

    The collected data can be stored as a singel dataframe containing (not split), split into
    separate dataframes by class or by family. Excluded classes are CAZy classes not specified in
    the configuration file and thus, will not be scraped. User_cazy_families is the list of CAZy
    families specified to be scraped in the configration file.
    """
    # Program preparation
    time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")  # used in naming files

    if argv is None:
        parser = utilities.build_parser()
        args = parser.parse_args()
    else:
        args = utilities.build_parser(argv).parse_args()

    if logger is None:
        logger = utilities.build_logger("cazy_webscraper", args)
    logger.info("Run initiated")

    if args.output is not sys.stdout:
        file_io.make_output_directory(args.output, logger, args.force, args.nodelete)

    if args.genbank is not None:
        # create directory to write FASTA files to
        if (args.genbank_output is not sys.stdout) and (args.genbank_output != args.output):
            file_io.make_output_directory(args.genbank_output, logger, args.force, args.nodelete)

    if (args.pdb is not None) and (args.pdb_output != args.output):
        # create directory to write structure files to
        file_io.make_output_directory(args.pdb_output, logger, args.force, args.nodelete)

    if args.subfamilies is True:
        logger.warning("Enabled to retrieve subfamilies")

    max_tries = args.retries + 1  # maximum number of times to try scraping a CAZy

    # build database and return open database session
    try:
        session = sql.build_db(time_stamp, args, logger)
    except Exception:
        logger.error("Failed to build SQL database. Terminating program", exc_info=1)
        sys.exit(1)

    # retrieve configuration data
    file_io_path = file_io.__file__
    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args,
        logger,
    )

    cazy_home = "http://www.cazy.org"

    get_cazy_data(
        cazy_home,
        excluded_classes,
        config_dict,
        cazy_dict,
        max_tries,
        time_stamp,
        session,
        logger,
        args,
    )

    logger.info("Finished program preparation")
    logger.info(
            "Starting retrieval of data from lass_url cwebsite.\n"
            "Thank you for using the cazy_webscraper.py\n"
            "Terminating program"
        )


def get_cazy_data(
    cazy_home,
    excluded_classes,
    config_dict,
    cazy_dict,
    max_tries,
    time_stamp,
    session,
    logger,
    args,
):
    """Coordinate retrieval of data from the CAZy website.

    This function coordinates the crawling through the CAZy website by calling the appropriate
    functions, and then retrieving the protein data by calling to the appropriate data again.

    :param cazy_home: str, url of CAZy home page
    :param excluded_classes: list, list of classes to not scrape from CAZy
    :param config_dict: dict, user defined configuration of the scraper
    :param cazy_dict: dict, dictionary of excepct CAZy synonyms for CAZy classes
    :param max_tries: int, maximum number of times to scrape CAZy if errors are encountered
    :param time_stamp: str, data and time scrape was initiated
    :param session: session, open database session
    :param logger: logger object
    :param args: cmd args parser

    Return nothing.
    """
    # List of urls that were failed to be scraped
    failed_url_scrapes = []

    # List of proteins that were not added to the database and attempting to do so raised an error
    sql_failures = []

    # retrieve links to CAZy class pages, return list of CazyClass objects
    cazy_classes = crawler.get_cazy_class_urls(
        cazy_home,
        excluded_classes,
        max_tries,
        cazy_dict,
        logger,
    )

    logger.info("Starting retrieval of CAZy families")

    # scrape each retrieved class page
    for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):

        # first attempt of scraping, retrieve URLs to CAZy families
        if len(cazy_class.failed_families.keys) == 0:

            # retrieve URLs to families for the current CAZy class, return as Family class objects
            class_families, error_message, incorrect_urls = crawler.get_cazy_family_urls(
                cazy_class.url,
                cazy_class.name,
                cazy_home,
                args,
                logger,
            )

            if incorrect_urls is not None:
                failed_url_scrapes += incorrect_urls

            if class_families is None:  # couldn't retrieve URLs to families for working CAZy class
                # add one to the number of scrapping attempts
                cazy_class.tries += 1

                # check if maximum number of attempts to connect have been met
                if cazy_class.tries == max_tries:
                    failed_url_scrapes += (
                        f"{cazy_class.url}\t"
                        f"{cazy_class.name}\t"
                        "No CAZy familes from this class were scraped\t"
                        f"{error_message}"
                    )
                    continue

                else:
                    cazy_classes += cazy_class  # retry scraping Class page after the other Classes
                    continue

        # Not first try, scrape only the families for which a connections to CAZy previously failed
        else:
            class_families = cazy_class.failed_families.keys()

        # Scrape the familes of the current CAZy class, retrieving protein data

        if (config_dict is None) or (config_dict[cazy_class.name] is None):
            # No (sub)families were specified, therefore, scraping all families of the CAZy class

            for cazy_family in tqdm(class_families, desc=f"Parsing {cazy_class.name} families"):
                # Populate family with Proteins catalogued under the CAZy family
                family, failed_family_page_scrapes, family_sql_failures = crawler.parse_family(
                    cazy_family,
                    cazy_home,
                    logger,
                    session,
                )

                # if failed to scrape some pages for CAZy family
                if failed_family_page_scrapes is not None:
                    # add one to the number of attempted scrapes for CAZy family
                    try:
                        cazy_class.failed_families[cazy_family] += 1
                    except KeyError:
                        cazy_class.failed_families[cazy_family] = 1  # first attempt

                    # check if max number of attempts to connect family pages has been met
                    if cazy_class.failed_families[cazy_family] == max_tries:
                        failed_url_scrapes += failed_family_page_scrapes  # store urls
                        del cazy_class.failed_families[cazy_family]  # do not try another scrape
                        continue

                if len(cazy_class.failed_families.keys) != 0:
                    # if there are families with previously failed connection attempts
                    # and remaining tries, retry connection after working through other classes
                    cazy_classes += cazy_class

                sql_failures += family_sql_failures

        else:
            # scrape only (sub)families specified in the config file

            for cazy_family in tqdm(class_families, desc=f"Parsing {cazy_class.name} families"):

                # Allow retrieval of subfamilies when only the parent CAZy family was named in the
                # config file, by searching by the family not subfamily in the config file
                if (args.subfamilies is True) and (family.name.find("_") != -1):
                    name_check = family.name[: (family.name.find("_"))]
                else:
                    name_check = family.name

                if name_check not in config_dict[cazy_class.name]:
                    continue

                # Populate family with Proteins catalogued under the CAZy family
                family, failed_family_page_scrapes, family_sql_failures = crawler.parse_family(
                    cazy_family,
                    cazy_home,
                    logger,
                    session,
                )

                # if failed to scrape some pages for CAZy family
                if failed_family_page_scrapes is not None:
                    # add one to the number of attempted scrapes for CAZy family
                    try:
                        cazy_class.failed_families[cazy_family] += 1
                    except KeyError:
                        cazy_class.failed_families[cazy_family] = 1  # first attempt

                    # check if max number of attempts to connect family pages has been met
                    if cazy_class.failed_families[cazy_family] == max_tries:
                        failed_url_scrapes += failed_family_page_scrapes  # store urls
                        del cazy_class.failed_families[cazy_family]  # do not try another scrape
                        continue

                if len(cazy_class.failed_families.keys) != 0:
                    # if there are families with previously failed connection attempts
                    # and remaining tries, retry connection after working through other classes
                    cazy_classes += cazy_class

                sql_failures += family_sql_failures

    # write out URLs which failed to be scaped
    if len(failed_url_scrapes) != 0:
        file_io.write_out_failed_scrapes(failed_url_scrapes, time_stamp, args, logger)

    # write out Proteins which failed to be be added to the database
    if len(sql_failures) != 0:
        file_io.write_out_failed_proteins(sql_failures, time_stamp, args, logger)

    return


if __name__ == "__main__":
    main()
