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

:cmd_args:...

:func ...:...

Produces a dataframe containing protein data and write protein
sequences to FASTA files.
"""

import logging
import re
import sys

from typing import List, Optional

import mechanicalsoup

from scraper.parse import parse_cazy_protein_data
from scraper.file_io import make_output_directory
from scraper.utilities import build_parser, build_logger


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate retrieval of links to CAZy website pages,
    required for retrieval of protein data.

    :param argv: (optional) args_parser.Namespace object
    """
    # Program preparation

    # Build args parser
    # Check if Namespace is passed, if not parse cmd-line
    if argv is None:
        # parse cmd-line
        parser = build_parser()
        args = parser.parse_args()
    else:
        args = build_parser(argv).parse_args()

    # Build logger
    # Log file only created if specified at cmd-line
    if logger is None:
        logger = build_logger("cazy_webscraper", args)
    logger.info("Run initated")

    # If specified output directory, create output directory
    if args.outdir is not sys.stdout:
        make_output_directory(args.outdir, logger, args.force, args.nodelete)

    # page to start browser at: the CAZy homepage
    base_url = "http://www.cazy.org"

    # create browser object
    browser = mechanicalsoup.Browser()

    # Retrieve all links from the CAZy homepage
    links_dict = get_all_homepage_links(browser, base_url)

    # tuple of CAZy classes and abbrievations
    cazy_classes = [
        ["Glycoside-Hydrolases", "GH"],
        ["GlycosylTransferases", "GT"],
        ["Polysaccharide-Lyases", "PL"],
        ["Carbohydrate-Esterases", "CE"],
        ["Auxiliary-Activities", "AA"],
        ["Carbohydrate-Binding Modules", "CBM"],
    ]

    # Navigate through each CAZy class main page
    index = 0

    for index in range(len(cazy_classes)):  # for each specificed class
        cazy_class = cazy_classes[index][0]
        class_abbrev = cazy_classes[index][1]

        # compile URL for the class main page
        class_url = base_url + "/" + links_dict[cazy_class]
        # retrieve URLs for each family's main/summary page
        family_links = get_family_links(browser, class_url, class_abbrev, args)

        for family in family_links:
            family_url = base_url + "/" + family
            # get the link to the page for 'all' proteins catalogued in the family
            get_family_table_links(browser, family_url, base_url)

    # site is now populated with class, family and protein_table pages
    # parse protein tables and write to .csv files, with data separation specified by user
    # handled by parse module
    parse_cazy_protein_data(site, cazy_classes, args, logger)


def get_all_homepage_links(browser, base_url):
    """Retrieve all links from the homepage.

    Return dictionary of link.text : link key/value pairs.
    """
    # create response object, <Response [200]> is a successful connection
    home_page = browser.get(base_url)

    # obtain links on homepage
    all_links = home_page.soup.select("a")
    # empty dictionary to store links in
    # as text : url_address key/value pairs
    link_dict = {}

    for link in all_links:
        try:
            link_dict[link.text] = link["href"]
        except KeyError:
            pass

    return link_dict


def get_family_links(browser, class_url, class_abbreviation, args):
    """Navigate CAZy class page, iterating through pages listing CAZymes.

    Return list of links to CAZy family pages for given class.
    """
    class_page = browser.get(class_url)

    # obtain links on class main page
    all_links = class_page.soup.select("a")

    # empty lists to store links to family pages
    family_links = []

    # search pattern to determine if link is for CAZy family or not
    pattern = re.compile(rf"{class_abbreviation}\d+?.*?\.html")

    # retieve all links from CAZy class main page
    for link in all_links:
        # retrieve the link from the bs4 object
        try:
            # check link is for a family/subfamily page
            address = link["href"]
            search_result = re.match(pattern, link.text)
            if search_result is True:  # link is a CAZy family/subfamily page link

                if args.subfamily is False:  # do not retrieve subfamily links
                    # subfamily names contain '_'
                    subfam_index = address.find("_")
                    if subfam_index != -1:  # link is for family page
                        family_links.append(address)

                else:  # retrieve all family and subfamily page links
                    family_links.append(address)
                    # duplicate protein results caused by retrieving family and
                    # subfamily data are removed when parsing the pages
                    # becuase some proteins catalogued under the family with
                    # no subfamily
        except KeyError:
            pass

    return family_links


def get_family_table_links(browser, family_url, base_url):
    """Retrieves the links for all tables containing proteins for the family.

    :param browser: Beautifulsoup browser objecct
    :param family_link: Beautifulsoup link object

    Return a list of urls for all pages containing protein data.
    """
    family_page = browser.get(family_url)

    # obtain all links on family main/summary page
    all_links = family_page.soup.select("a")

    # empty dictionary to store links:
    all_links_dict = {}

    for link in all_links:
        try:
            all_links_dict[link.text] = link["href"]
        except KeyError:
            pass

    # retrieves full url not only the suffix
    family_all_url = all_links_dict["All"]

    # navigate to family 'all' page
    family_all_page = browser.get(family_all_url)

    # obtain all links on family 'all' page
    all_links = family_all_page.soup.select("a")

    all_link_dict = {}
    for link in all_links:
        try:
            all_link_dict[link.text] = link["href"]
        except KeyError:
            pass

    # retrieve all links to the pages containing protein tables/data
    # for the family
    table_pages_urls = []
    table_pages_urls.append(family_all_url)  # contains first table

    # iterate through the page numbers to retrieve the links
    # of all the pages containing protein tables for the family
    # These cannot be retrieved in on go as each page contains
    # a limited number of links to the other table pages.
    page_count = 2
    error_check = 0
    while error_check == 0:
        try:
            new_url = base_url + "/" + all_link_dict[str(page_count)]
            table_pages_urls.append(new_url)

            # open new url to retrieve next page and
            # retrieve all links from next page
            # so can repeat and retrieve link for next page
            new_page = browser.get(new_url)
            all_links = new_page.soup.select("a")

            all_link_dict = {}
            for link in all_links:
                try:
                    all_link_dict[link.text] = link["href"]
                except KeyError:
                    pass

            page_count += 1

        except KeyError as error:
            error_check = -1
    return table_pages_urls


if __name__ == "__main__":
    main()
