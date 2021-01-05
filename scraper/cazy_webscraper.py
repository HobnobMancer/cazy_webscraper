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

from scraper import crawler, file_io, parse, utilities


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

    if args.pdb is not None:
        # create directory to write structure files to
        if (args.pdb_output is not None) and (args.pdb_output != args.output):
            file_io.make_output_directory(args.pdb_output, logger, args.force, args.nodelete)

    if args.subfamilies is True:
        logger.warning("Enabled to retrieve subfamilies")

    max_tries = args.retries + 1  # maximum number of times to try scraping a CAZy

    # retrieve configuration data
    file_io_path = file_io.__file__
    excluded_classes, config_dict, cazy_dict = file_io.parse_configuration(
        file_io_path,
        args,
        logger,
    )

    logger.info("Finished program preparation")
    logger.info("Starting retrieval of data from CAZy")

    # Crawl through and scrape CAZy website/database
    cazy_home = "http://www.cazy.org"  # the CAZy homepage URL

    # Retrieve data from CAZy database
    get_cazy_data(
        cazy_home,
        excluded_classes,
        config_dict,
        cazy_dict,
        max_tries,
        time_stamp,
        logger,
        args,
    )

    logger.info(
        (
            "Finished scraping the CAZy website.\n"
            "Thank you for using the cazy_webscraper.py\n"
            "Terminating program"
        )
    )


def get_cazy_data(
    cazy_home,
    excluded_classes,
    config_dict,
    cazy_dict,
    max_tries,
    time_stamp,
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
    :param logger: logger object
    :param args: cmd args parser

    Return nothing.
    """
    # List of urls that were failed to be scraped
    failed_url_scrapes = []

    # retrieve links to CAZy class pages
    # nested lists of [class_url, number_of_attempted_class_url_scrapes]
    class_urls = get_class_urls(cazy_home, excluded_classes, max_tries, logger)

    all_data = []  # stores all Family class objects if not splitting the data

    logger.info("Starting retrieval of CAZy families")

    # scrape each retrieved class page
    for class_url in tqdm(class_urls, desc="Parsing CAZy classes"):

        # retrieve class name from url
        class_name = class_url[0][20:-5]
        #  convert retrieved class name to synonym used in configuration file
        for key in cazy_dict:
            if class_name in cazy_dict[key]:
                class_name = key

        # retrieve URLs to families under current working CAZy class
        # nested list of [family_url, number_of_attempts_to_scrape_family_all_page]
        family_urls = crawler.get_cazy_family_urls(
            class_url[0],
            cazy_home,
            class_name,
            args,
            logger,
        )

        if family_urls is None:  # couldn't retrieve URLs to family for working CAZy class
            # add one to the number of scrapping attempts
            class_url[1] += 1

            if class_url[1] == max_tries:  # max number of attempts to scrape has been met
                failed_url_scrapes += (
                    f"{class_url[0]} - no CAZy familes from this class were scraped"
                )
                continue
            else:
                class_urls += class_url  # retry scraping Class page after the other Class URLs
                continue

        families = []  # store Family class objects if splitting data be class

        logger.info("Starting retrieval of protein records of protein records from families")

        # Scrape the familes of the current CAZy class, retrieving protein data

        if (config_dict is None) or (config_dict[class_name] is None):
            # No (sub)families were specified, therefore, scraping all families of the CAZy class

            for family_url in tqdm(family_urls, desc="Parsing CAZy families"):
                # check url format is correct
                try:
                    re.match(
                        r"http://www.cazy.org/(\D{2,3})(\d+|\d+_\d+).html", family_url[0]
                    ).group()
                except AttributeError:
                    logger.warning(
                        f"Formate of URL {family_url[0]} is incorrect.\n"
                        "Will not attempt to scrape this URL."
                    )
                    failed_url_scrapes += f"{family_url[0]} - url format was incorrect"
                    continue

                family = None
                family_name = family_url[0][(len(cazy_home) + 1): -5]
                # build family object, populated by Proteins catalogued under the CAZy family
                family = crawler.parse_family(family_url[0], family_name, cazy_home, logger)

                # [family_object, error]
                if family[1] is not None:  # Scraping family '_all' page was unsuccessful
                    # add one to the number of attempted scrapes the CAZy family's '_all' page
                    family_url[1] += 1

                    if family_url[1] == max_tries:  # max number of scraping attempts reached
                        failed_url_scrapes += (
                            f"{family_url[0]} - the following error was raised {family[2]}"
                        )
                        continue
                    else:
                        family_urls += family_url
                        continue

                # store the family if scraped successfully
                if args.data_split == "family":
                    logger.info(f"Data split by Family. Writing out df for {family_name}")
                    parse.proteins_to_dataframe([family[0]], time_stamp, args, logger)
                else:
                    families.append(family[0])

        else:
            # scrape only (sub)families specified in the config file

            for family_url in tqdm(family_urls, desc="Parsing CAZy families"):
                # check url format is correct
                try:
                    re.match(
                        r"http://www.cazy.org/(\D{2,3})(\d+|\d+_\d+).html", family_url[0]
                    ).group()
                except AttributeError:
                    logger.warning(
                        (
                            f"Formate of URL {family_url[0]} is incorrect.\n"
                            "Will not attempt to scrape this URL."
                        )
                    )
                    continue

                family = None
                family_name = family_url[0][(len(cazy_home) + 1) : -5]

                # Allows retrieval of subfamilies when only the parent CAZy family was named in the
                # config file
                if (args.subfamilies is True) and (family_name.find("_") != -1):
                    name_check = family_name[: (family_name.find("_"))]
                else:
                    name_check = family_name

                if name_check in config_dict[class_name]:
                    # build family object, populated by Proteins catalogued under the CAZy family
                    family = crawler.parse_family(family_url[0], family_name, cazy_home, logger)

                    # [family_object, error]
                    if family[1] is not None:  # Scraping family '_all' page was unsuccessful
                        # add one to the number of attempted scrapes the CAZy family's '_all' page
                        family_url[1] += 1

                        if family_url[1] == max_tries:  # max number of scraping attempts reached
                            failed_url_scrapes += (
                                f"{family_url[0]} - the following error was raised {family[2]}"
                            )
                            continue
                        else:
                            family_urls += family_url
                            continue

                    # store the family if scraped successfully
                    if args.data_split == "family":
                        logger.info(f"Data split by Family. Writing out df for {family_name}")
                        parse.proteins_to_dataframe([family[0]], time_stamp, args, logger)
                    else:
                        families.append(family[0])

        if args.data_split == "class":
            if len(families) != 0:
                logger.info(f"Data split by Class. Writing out df for {class_name}")
                parse.proteins_to_dataframe(families, time_stamp, args, logger)
            else:
                logger.warning(f"Didn't retrieve any families for {class_name}")

        elif args.data_split is None:
            all_data += families

    if args.data_split is None:
        if len(all_data) != 0:
            logger.info("Data was not split. Writing all retrieved data to a single df")
            parse.proteins_to_dataframe(all_data, time_stamp, args, logger)
        else:
            logger.warning("Didn't retrieve any protein data from CAZy")

    # write out URLs which failed to be scaped
    file_io.write_out_failed_scrapes(failed_url_scrapes, time_stamp, args, logger)

    return


def get_class_urls(cazy_home, excluded_classes, max_tries, logger):
    """Retrieve class urls, add storage of number of attempted scrapes of each class URL.

    :param cazy_url: str, URL to the CAZy home page.
    :param excluded_classes: list, list of CAZy classes not to be scraped
    :param max_tries: int, maximum number of times to try scrape if errors are encountered
    :param logger: logger object

    Return list of CAZy class URLs. Each item is a list of [URL, 0]
    - 0 is used to count number of attempted connections.
    """

    class_urls = crawler.get_cazy_class_urls(cazy_home, excluded_classes, max_tries, logger)

    try:
        if len(class_urls) == 0:
            logger.error("Failed to retrieve URLs to CAZy class pages.\nTerminating program")
            sys.exit(1)
    except TypeError:  # rased when class_pages is None
        logger.error("Failed to retrieve URLs to CAZy class pages.\nTerminating program")
        sys.exit(1)

    # add storing the number of times an attempt to scrape the class page as been performed
    index = 0
    for index in range(len(class_urls)):
        # list structure [class_url, number_of_tried_scrapes]
        class_urls[index] = [class_urls[index], 0]

    return class_urls


if __name__ == "__main__":
    main()
