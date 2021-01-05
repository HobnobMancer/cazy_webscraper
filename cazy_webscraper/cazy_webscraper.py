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
import time

import numpy as np

from collections import defaultdict
from typing import List, Optional
from requests.exceptions import ConnectionError, MissingSchema
from urllib3.exceptions import HTTPError, RequestError

import mechanicalsoup

from tqdm import tqdm

from cazy_webscraper import crawler, file_io, parse, utilities


class Protein:
    """A single protein.

    Each protein has a name, source organism (source), and links to external databases. The links to
    external databases are stored in a dictionary, keyed by the external database name ('str') with
    'list' values becuase there may be multiple links per database.

    Multiple 'synonym' GenBank accession numbers maybe listed for a single protein. CAZy only
    hyperlinks the first listed accession number. This accession is the one listed for the protein,
    because is presumed to be the accession used by CAZy in their classification. All other listed
    GenBank accessions are regarded as synonyms, including for example splice variants and identical
    protein sequence submissions.
    """

    def __init__(self, name, family, ec, source, links=None, genbank_synonyms=None):
        self.name = name
        self.family = family
        self.ec = ec
        self.source = source
        if links is None:
            self.links = defaultdict(list)
        else:
            self.links = links
        self.genbank_synonyms = genbank_synonyms

    def __str__(self):
        """Create representative string of class object"""
        return f"{self.name} ({self.family} {self.source}): links to {self.links.keys()}"

    def __repr__(self):
        """Create representative object"""
        return (
            f"<Protein: {id(self)}: {self.name}, {self.family} "
            f"({self.source}), {len(self.links)} to external databases>"
        )

    def get_protein_dict(self):
        """Return a dictionary containing all the data of the protein."""
        protein_dict = {"Protein_name": [self.name], "CAZy_family": [self.family]}

        if len(self.ec) == 0:
            protein_dict["EC#"] = [np.nan]
        elif len(self.ec) == 1:
            protein_dict["EC#"] = self.ec
        else:
            ec_string = "\n".join(self.ec)
            protein_dict["EC#"] = [ec_string]

        protein_dict["Source_organism"] = [self.source]

        if type(self.links) is dict:
            for database in ["GenBank", "UniProt", "PDB/3D"]:
                try:
                    if len(self.links[database]) == 1:
                        protein_dict[database] = self.links[database]
                    else:
                        accession_string = ",\n".join(self.links[database])
                        protein_dict[database] = [accession_string]
                except KeyError:
                    protein_dict[database] = [np.nan]
        else:
            for database in ["GenBank", "UniProt", "PDB/3D"]:
                protein_dict[database] = [np.nan]
        return protein_dict


class Family:
    """A single CAZy family."""

    members = set()  # holds Protein instances

    def __init__(self, name, cazy_class):
        self.name = name
        self.cazy_class = cazy_class

    def __str__(self):
        return f"CAZy family {self.name}: {len(self.members)} protein members"

    def __repr__(self):
        return f"<Family: {id(self)}: {self.name}, {len(self.members)} protein members"

    def get_proteins(self):
        """Return a list of all protein members of the CAZy family."""
        return self.members

    def get_family_name(self):
        """Return family name"""
        return self.name


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy.

    The collected data can be stored as a singel dataframe containing (not split), split into
    separate dataframes by class or by family. Excluded classes are CAZy classes not specified in
    the configuration file and thus, will not be scraped. User_cazy_families is the list of CAZy
    families specified to be scraped in the configration file.
    """
    # Program preparation
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
    get_cazy_data(cazy_home, excluded_classes, config_dict, cazy_dict, max_tries, logger, args)

    logger.info(
        (
            "Finished scraping the CAZy website.\n"
            "Thank you for using the cazy_webscraper.py\n"
            "Terminating program"
        )
    )


def get_cazy_data(cazy_home, excluded_classes, config_dict, cazy_dict, max_tries, logger, args):
    """Coordinate retrieval of data from the CAZy website.

    This function coordinates the crawling through the CAZy website by calling the appropriate
    functions, and then retrieving the protein data by calling to the appropriate data again.

    :param cazy_home: str, url of CAZy home page
    :param excluded_classes: list, list of classes to not scrape from CAZy
    :param config_dict: dict, user defined configuration of the scraper
    :param cazy_dict: dict, dictionary of excepct CAZy synonyms for CAZy classes
    :param max_tries: int, maximum number of times to scrape CAZy if errors are encountered
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
                    parse.proteins_to_dataframe([family[0]], args, logger)
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
                        parse.proteins_to_dataframe([family[0]], args, logger)
                    else:
                        families.append(family[0])

        if args.data_split == "class":
            if len(families) != 0:
                logger.info(f"Data split by Class. Writing out df for {class_name}")
                parse.proteins_to_dataframe(families, args, logger)
            else:
                logger.warning(f"Didn't retrieve any families for {class_name}")

        elif args.data_split is None:
            all_data += families

    if args.data_split is None:
        if len(all_data) != 0:
            logger.info("Data was not split. Writing all retrieved data to a single df")
            parse.proteins_to_dataframe(all_data, args, logger)
        else:
            logger.warning("Didn't retrieve any protein data from CAZy")

    # write out URLs which failed to be scaped
    file_io.write_out_failed_scrapes(failed_url_scrapes, args, logger)

    return


def get_class_urls(cazy_home, excluded_classes, max_tries, logger):
    """Retrieve class urls, add storage of number of attempted scrapes of each class URL.

    :param cazy_url: str, URL to the CAZy home page.
    :param excluded_classes: list, list of CAZy classes not to be scraped
    :param max_tries: int, maximum number of times to try scrape if errors are encountered
    :param logger: logger object

    Return list of CAZy class URLs.
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


def browser_decorator(func):
    """Decorator to retry the wrapped function up to 'retries' times."""

    def wrapper(*args, retries=10, **kwargs):
        tries, success, err = 0, False, None
        while not success and (tries < retries):
            try:
                response = func(*args, **kwargs)
            except (
                ConnectionError,
                HTTPError,
                OSError,
                MissingSchema,
                RequestError,
            ) as err_message:
                success = False
                response = None
                err = err_message
            if response is not None:  # response was successful
                success = True
            # if response from webpage was not successful
            tries += 1
            time.sleep(10)
        if (not success) or (response is None):
            return [None, err]
        else:
            return [response, None]

    return wrapper


@browser_decorator
def get_page(url):
    """Create browser and use browser to retrieve page for given URL.

    :param url: str, url to webpage

    Return browser response object (the page).
    """
    # create browser object
    browser = mechanicalsoup.Browser()
    # create response object
    page = browser.get(url)
    page = page.soup

    return page


if __name__ == "__main__":
    main()
