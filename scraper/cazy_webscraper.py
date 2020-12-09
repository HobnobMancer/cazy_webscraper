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
:func get_cazy_class_urls: retrieve URLs to CAZy class summary pages from CAZy homepage
:func get_cazy_family_urls: retrieve URLs to families on CAZy class summary page
:func get_subfamily_links: retrieve URLs to subfamilies on CAZy class summary page
:func parse_family: build Family class object to represent CAZy family
:func parse_family_pages: retrieve all URLs to pages containing protein records for CAZy family
:func parse_proteins: retrieve protein records from protein table page
:func row_to_protein: parse the protein record to build a Protein class object
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

from scraper import utilities, file_io, parse


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
        family_urls = get_cazy_family_urls(class_url[0], cazy_home, class_name, args, logger)

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
                    re.match(rf"{cazy_home}/(\D{2,3})(\d+|\d+_\d+).html", family_url[0]).group()
                except AttributeError:
                    logger.warning(
                        f"Formate of URL {family_url[0]} is incorrect.\n"
                        "Will not attempt to scrape this URL."
                    )
                    failed_url_scrapes += f"{family_url[0]} - url format was incorrect"
                    continue

                family = None
                family_name = family_url[0][(len(cazy_home) + 1) : -5]
                # build family object, populated by Proteins catalogued under the CAZy family
                family = parse_family(family_url[0], family_name, cazy_home, logger)

                # [family_object, error]
                if len(family[1]) is not None:  # Scraping family '_all' page was unsuccessful
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
                    re.match(rf"{cazy_home}/(\D{2,3})(\d+|\d_\d+).html", family_url).group()
                except AttributeError:
                    logger.warning(
                        (f"Formate of URL {family_url} is incorrect.\n" "Will not attempt to scrape this URL.")
                    )
                    continue

                family = None
                family_name = family_url[(len(cazy_home) + 1) : -5]

                # Allows retrieval of subfamilies when only the parent CAZy family was named in the
                # config file
                if (args.subfamilies is True) and (family_name.find("_") != -1):
                    name_check = family_name[: (family_name.find("_"))]
                else:
                    name_check = family_name

                if name_check in config_dict[class_name]:
                    # build family object, populated by Proteins catalogued under the CAZy family
                    family = parse_family(family_url[0], family_name, cazy_home, logger)

                    # [family_object, error]
                    if len(family[1]) is not None:  # Scraping family '_all' page was unsuccessful
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

    return


def get_class_urls(cazy_home, excluded_classes, max_tries, logger):
    """Retrieve class urls, add storage of number of attempted scrapes of each class URL.

    :param cazy_url: str, URL to the CAZy home page.
    :param excluded_classes: list, list of CAZy classes not to be scraped
    :param max_tries: int, maximum number of times to try scrape if errors are encountered
    :param logger: logger object

    Return list of CAZy class URLs.
    """

    class_urls = get_cazy_class_urls(cazy_home, excluded_classes, max_tries, logger)

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


def get_cazy_class_urls(cazy_home, excluded_classes, max_tries, logger):
    """Returns a list of CAZy class main/home page URLs for each specified class as the CAZy site.

    :param cazy_url: str, URL to the CAZy home page.
    :param excluded_classes: list, list of CAZy classes not to be scraped
    :param max_tries: int, maximum number of times to try scrape if errors are encountered
    :param logger: logger object

    Return list of URLs.
    """
    logger.info("Retrieving URLs to summary CAZy class pages")

    # define items to be excluded from returned class list, ALWAYS exlide links to genomes
    if excluded_classes is not None:
        exclusions = tuple(["<strong>Genomes</strong>"] + excluded_classes)
    else:
        exclusions = "<strong>Genomes</strong>"

    home_page = [None, None]
    tries = 0  # number of the try/attempt

    while (home_page[0] is None) and (tries < max_tries):
        home_page = get_page(cazy_home)

        if (home_page[0] is None) and (tries < max_tries):
            logger.error(
                f"Failed to connect to CAZy homepage after 10 attempts,\n"
                f"for attempt# {(tries+1)}/{max_tries}\n"
                "The following error was raised:\n"
                f"{home_page[1]}\n"
                f"Reattempting for attempt# {(tries+2)} in 10s."
            )
            time.sleep(10)
            tries += 1

    if home_page[0] is None:
        logger.error(
            (
                "Failed to connect to CAZy home-page after multiple attempts.\n"
                "The following error was raised:\n"
                f"{home_page[1]}"
                "Could not retrieve URLs to CAZy classes.\n"
                "Check the network connection.\n"
                "Terminating program."
            )
        )
        sys.exit(1)

    try:
        return [
            f"{cazy_home}/{_['href']}"
            for _ in home_page[0].find_all("a", {"class": "spip_out"})
            if (not _["href"].startswith("http")) and (str(_.contents[0]) not in exclusions)
        ]
    except AttributeError:  # raise if can't find results with find_all("a", {"class": "spip_out"})
        return None


def get_cazy_family_urls(class_url, cazy_home, class_name, args, logger):
    """Retrieve all protein members of each CAZy family within the given CAZy class.

    :param class_url: str, URL to the CAZy class
    :param cazy_home: str, URL to CAZy home page
    :param class_name: str, name of CAZy class
    :param args: args parser object
    :param logger: logger object

    Returns list of URLs to family pages.
    """
    logger.info(f"Retrieving URLs to families under {class_name}")

    # scrape the class page
    class_page = get_page(class_url)

    if class_page[0] is None:
        logger.error(
            (
                f"Failed to connect to {class_url} after 10 attempts.\n"
                "The following error was raised:\n"
                f"{class_page[1]}"
                "Could not retrieve URLs to CAZy famileis for this class.\n"
                "This class will be skipped during the scraping process."
            )
        )
        return None

    # retrieve the <h3> element that titles the div section containing the tables of family links
    family_h3_element = [
        _
        for _ in class_page[0].find_all("h3", {"class": "spip"})
        if str(_.contents[0]).strip() == "Tables for Direct Access"
    ][0]

    # retrieve all tables within the parent div section of the <h3> element
    tables = family_h3_element.parent.find_all("table")

    # tables[0] is the table containing links to CAZy families
    # tables[1] is the table containing the link to unclassified proteins

    family_urls = family_urls = [f"{cazy_home}/{_['href']}" for _ in tables[0].find_all("a")]
    family_urls.append(f"{cazy_home}/{tables[1].a['href']}")
    if args.subfamilies is True:
        subfam_urls = get_subfamily_links(family_h3_element, cazy_home, logger)
        if subfam_urls is not None:
            family_urls += subfam_urls

    if family_urls is not None:
        # add storage of number of attempts to scrape the Family _all page
        index = 0
        for index in range(len(family_urls)):
            family_urls[index] = [family_urls[index], 0]

    return family_urls


def get_subfamily_links(family_h3_element, cazy_home, logger):
    """Retrieve URL links to CAZy subfamilies.

    :param family_h3_element: bs4.element.Tag, h3 element titling the page div
    :param cazy_home: str, URL to CAZy home_page
    :param logger: logger object

    Return list of URLs to subfamilies.
    """
    parent_div = family_h3_element.parent
    all_links = parent_div.find_all("a")

    pattern = re.compile(r"\D+?\d+?_\d+?\.html")

    urls = []  # empty list to store subfamily URLs

    for link in all_links:
        try:
            search_result = re.search(pattern, link["href"])
            urls.append(f"{cazy_home}/{search_result.group()}")
        except (KeyError, AttributeError) as error:
            # KeyError raised if link does not have ['href']
            # AttributeError error raised if search_result is None becuase not subfam link
            pass

    if len(urls) == 0:
        return
    else:
        return urls


def parse_family(family_url, family_name, cazy_home, logger):
    """Returns a Family object with Protein members, scraped from CAZy.

    :param family_url: str, URL to CAZy family summary page
    :param family_name: str, name of CAZy family
    :param cazy_home: str, URL to CAZy home page
    :param logger: logger object

    Return list of Family object encapsulating data from the passed family
    summary page, and an error code if failed to scrape the family 'all' page.
    """
    logger.info(f"Starting retrieval of proteins for {family_name} from {family_url}")

    # retrieve class from family name
    pattern = re.compile(r"\D+")
    search_result = re.match(pattern, family_name)
    cazy_class = search_result.group()

    family = Family(family_name, cazy_class)
    family.members = set()

    # compile URL to first family page of protein records
    first_pagination_url = family_url.replace(".html", "_all.html")

    # check url formating
    try:
        re.match(rf"{cazy_home}/\D{2,3}\d+_?\d+_all.html", first_pagination_url).group()
    except AttributeError:
        logger.warning(
            f"Incorrect formatting of first protein table page URL: {first_pagination_url}\n"
            "Will not try and connect to this URL."
        )
        return family

    first_pagination_page = get_page(first_pagination_url)

    if first_pagination_page[0] is None:
        logger.warning(
            (
                f"Could not connect to {family_url}\n"
                "The following error was raised:\n"
                f"{first_pagination_page[1]}"
                f"No protein records for CAZy family {family_name} will be retried."
            )
        )
        return [family, first_pagination_page[1]]

    protein_page_urls = get_protein_page_urls(first_pagination_url, first_pagination_page[0], cazy_home)

    for protein in tqdm(
        (y for x in (parse_proteins(url, family_name, logger) for url in protein_page_urls) for y in x),
        total=len(protein_page_urls),
        desc=f"Scraping protein pages for {family_name}",
    ):
        family.members.add(protein)

    return [family, None]


def get_protein_page_urls(first_pagination_url, first_pagination_page, cazy_home):
    """Retrieve the URLs to all pages containing proteins for the current working family.

    :param first_pagination_url: str, URL to first page contaiing proteins
    :param first_pagination_page: BS4 object, first page containing proteins
    :param cazy_home: str, URL of CAZy homepage

    Return list of URLs.
    """
    protein_page_urls = [first_pagination_url]

    # retrieve the URL to the final page of protein records in the pagination listing
    try:
        last_pagination_url = first_pagination_page.find_all("a", {"class": "lien_pagination", "rel": "nofollow"})[-1]
    except IndexError:  # there is no pagination; a single-query entry
        last_pagination_url = None

    if last_pagination_url is not None:
        url_prefix = last_pagination_url["href"].split("PRINC=")[0] + "PRINC="
        last_princ_no = int(last_pagination_url["href"].split("PRINC=")[-1].split("#pagination")[0])
        url_suffix = "#pagination" + last_pagination_url["href"].split("#pagination")[-1]

        # Build list of urls to all pages in the pagination listing, increasing the PRINC increment
        protein_page_urls.extend(
            [f"{cazy_home}/{url_prefix}{_}{url_suffix}" for _ in range(1000, last_princ_no + 1000, 1000)]
        )

    return protein_page_urls


def parse_proteins(protein_page_url, family_name, logger):
    """Returns generator of Protein objects for all protein rows on a single CAZy family page.

    :param protein_page_url, str, URL to the CAZy family page containing protein records
    :param family_name: str, name of CAZy family
    :param logger: logger object

    Return generator object.
    """
    logger.info(f"Retrieving proteins from {protein_page_url}")

    protein_page = get_page(protein_page_url)

    if protein_page[0] is None:
        logger.warning(
            (
                f"Could not connect to {protein_page_url}\n"
                "The following error was raised:\n"
                f"{protein_page[1]}\n"
                f"No protein records from this page will be retried."
            )
        )
        return

    # retrieve protein record table
    protein_table = protein_page[0].find_all("table", {"class": "listing"})[0]
    protein_rows = [
        _ for _ in protein_table.descendants if (_.name == "tr") and ("id" not in _.attrs) and ("class" not in _.attrs)
    ]

    for row in protein_rows:
        yield row_to_protein(row, family_name)


def row_to_protein(row, family_name):
    """Returns a Protein object representing a single protein row from a CAZy family protein page.

    Each row, in order, contains the protein name, EC number, source organism, GenBank ID(s),
    UniProt ID(s), and PDB accession(s).

    :param row: tr element from CAZy family protein page
    :param family_name: str, name of CAZy family

    Return Protein instance.
    """
    # retrieve list of cells ('td' elements) in row
    tds = list(row.find_all("td"))

    protein_name = tds[0].contents[0].strip()
    source_organism = tds[2].a.get_text()
    ec_numbers = []
    links = {}

    all_links = None
    try:
        all_links = tds[1].find_all("a")
    except TypeError:  # raised when no EC numbers are listed
        pass

    if all_links is not None:
        for link in all_links:
            ec_numbers.append(link.text)
    else:
        ec_numbers = [np.nan]

    # test for len(tds[x].contents) in case there is no link,
    # the check of .name then ensures link is captured
    if len(tds[3].contents) and tds[3].contents[0].name == "a":
        links["GenBank"] = [f"{_.get_text()} {_['href']}" for _ in tds[3].contents if _.name == "a"]
    if len(tds[4].contents) and tds[4].contents[0].name == "a":
        links["UniProt"] = [f"{_.get_text()} {_['href']}" for _ in tds[4].contents if _.name == "a"]
    if len(tds[5].contents) and tds[5].contents[0].name == "a":
        links["PDB/3D"] = [f"{_.get_text()} {_['href']}" for _ in tds[5].contents if _.name == "a"]

    # Retrieve GenBank accession synonms
    try:
        if len(links["GenBank"]) != 0:
            # retrieve HTML link GenBank accession (CAZy only hyperlinks the first listed accession)
            key = links["GenBank"][0]
            key = key.split(" ")[0]  # remove the URL address

            # retrieve all GenBank accessions listed in td element
            all_genbank_accessions = tds[3].get_text(separator=" ")
            all_genbank_accessions = all_genbank_accessions.split(" ")

            if key in all_genbank_accessions:
                all_genbank_accessions.remove(key)

            # remove duplicates
            all_genbank_accessions = list(dict.fromkeys(all_genbank_accessions))

            if len(all_genbank_accessions) == 0:
                genbank_synonyms = None
            else:
                genbank_synonyms = {key: all_genbank_accessions}
        else:
            genbank_synonyms = None
    except KeyError:
        genbank_synonyms = None

    return Protein(protein_name, family_name, ec_numbers, source_organism, links, genbank_synonyms)


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
