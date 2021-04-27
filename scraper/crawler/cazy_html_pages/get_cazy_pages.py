#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
"""Retrieve webpages from CAZy and write to disk."""


import logging
import os
import re
import sys
import time

import pandas as pd

from datetime import datetime
from tqdm import tqdm

from scraper import crawler
from scraper.utilities import build_logger


def get_cazy_pages(
    args, cazy_home, time_stamp, excluded_classes, cazy_dict, config_dict, kingdoms, start_time,
):
    """Coordinate retrieving CAZy HTML webpages from CAZy and write to disk.

    :param args: cmd-line args parser
    :param cazy_home: str, CAZy homepage website
    :param time_stamp: str, time programme was invoked
    :param excluded_classes: list, list of CAZy classes not to scrape
    :param cazy_dict: dict of CAZy classes name synonyms
    :param config_dict: dict of user defined configuration of the scraper
    :param kingdoms: list of Taxonomy Kingdoms to retrieve pages for
    :param start_time: pd object, time programme was invoked

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    if args.output is not sys.stdout:
        out_log_path = args.output
    else:
        out_log_path = None

    connection_failures_logger = build_logger(
        out_log_path, f"CAZy_connection_failures_CW_{time_stamp}.log",
    )
    format_failures_logger = build_logger(
        out_log_path, f"Format_and_parsing_errors_CW_{time_stamp}.log",
    )

    # retrieve links to the CAZy class pages, returned as CazyClass objects

    cazy_classes = crawler.get_cazy_classes(cazy_home, excluded_classes, cazy_dict, args)

    for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):

        # first attempt of scraping, retrieve URLs to CAZy families
        if len(list(cazy_class.failed_families.keys())) == 0:

            # retrieve Family class instances, containing the Family page URL
            class_families, error_message, incorrect_urls = crawler.get_cazy_family_urls(
                cazy_class.url, cazy_class.name, cazy_home, args,
            )

            if incorrect_urls is not None:
                for url in incorrect_urls:
                    connection_failures_logger.warning(url)

            if class_families is None:  # couldn't retrieve URLs to families for working CAZy class
                # add one to the number of scrapping attempts
                cazy_class.tries += 1

                # check if maximum number of attempts to connect have been met
                if cazy_class.tries == (args.retries + 1):
                    connection_failures_logger.warning(
                        f"{cazy_class.url}\t{cazy_class.name}\t"
                        f"No CAZy familes from this class were scraped\t{error_message}"
                    )
                    continue

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
                    family, retry_scrape, failed_url_connections, failed_formats = parse_all_family(
                        family, cazy_home, args,
                    )

                else:
                    (
                        family, retry_scrape, failed_url_connections, failed_formats,
                    ) = parse_family_by_kingdom(
                        family, cazy_home, kingdoms, args,
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

                for error in failed_formats:
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
                    family, retry_scrape, failed_url_connections, failed_formats = parse_all_family(
                        family, cazy_home, args,
                    )

                else:
                    (
                        family, retry_scrape, failed_url_connections, failed_formats,
                    ) = parse_family_by_kingdom(
                        family, cazy_home, kingdoms, args,
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

                for error in failed_formats:
                    format_failures_logger.warning(error)

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
    sys.exit(1)


def parse_all_family(family, cazy_home, args):
    """Retrieve webpages for CAZy family 'all' tabs.

    CAZy families have separate series of HTML tables for each taxonomy Kingdom, and a single series
    of HTML tables containing all proteins from all Kingdoms, called 'all'. These 'all' pages
    containing 1000 proteins per page and thus scraping these pages is significantly faster than
    scraping all the specific Kingdom HTML tables (which hold only 100 proteins each, and thus
    require more calls to CAZy).

    :param family: Family class object, representation of CAZy family
    :param cazy_home: str, URL to CAZy home page
    :param args: cmd-line args parser

    Return:
        Family object,
        Boolean whether to scrape the family again,
        List of URLs which couldn't connect to CAZy (and no have re-try attempts left),
        List of pages whose parsing did match the expected format
    """
    logger = logging.getLogger(__name__)

    # define lists for storing error messages during scraping and parsing proteins
    failed_scrapes = []  # URLs of pages for which maximum number of scrape attempts is MET
    format_failures = []  # pages where no proteins were retrieved

    # check if there were pagination pages to which a connection could not be made previously

    if len(list(family.failed_pages.keys())) != 0:  # rescraping selected pages of the family
        protein_page_urls = list(family.failed_pages.keys())

        if len(protein_page_urls) == 1:
            first_url_template = family.url.replace(".html", "_all.html")
            if protein_page_urls[0] == first_url_template:
                protein_page_urls = get_pagination_pages(
                    protein_page_urls[0], family, cazy_home, args,
                )

    else:  # scraping for the first time
        # compile URL to first family page of protein records
        first_pagination_url = family.url.replace(".html", "_all.html")

        # check url formating of first pagination url
        try:
            re.match(
                r"http://www.cazy.org/\D{2,3}(\d+|\d+_\d+)_all.html", first_pagination_url
            ).group()
        except AttributeError:
            logger.warning(
                f"Incorrect formatting of first protein table page URL: {first_pagination_url}\n"
                "Will not try and connect to this URL."
            )
            format_failures.append(
                    f"{first_pagination_url}\tIncorrect URL format therefore could not retrieve "
                    f"proteins from CAZy family {family.name}\n"
            )
            retry_scrape = False
            return family, retry_scrape, failed_scrapes, format_failures

        # retrieve pages to other pagination pages of protein tables for the family
        protein_page_urls = get_pagination_pages(
            first_pagination_url, family, cazy_home, args,
        )

    # iterate through the pages and retrieve proteins, adding to the local CAZy data

    for protein_table_page in tqdm(
        (y for x in (
            get_html_page(url, family.name, args) for url in protein_page_urls
        ) for y in x),
        total=len(protein_page_urls),
        desc=f"Retrieving Family pages for {family.name}",
    ):
        if protein_table_page["url"] is not None:
            # Protein not retrieved because couldn't connect to CAZy
            try:
                family.failed_pages[protein_table_page["url"]] += 1
            except KeyError:
                family.failed_pages[protein_table_page["url"]] = 1  # First failed attempt

            if family.failed_pages[protein_table_page["url"]] == args.retries:
                # maximum attempts to connect have been reached no more attempts made, write to file
                failed_scrapes.append(
                    f"{protein_table_page['url']}\t{family.cazy_class}\t"
                    f"Failed to connect to this page of proteins for {family.name}, "
                    f"and raised the following error message:\n{protein_table_page['error']}"
                )
                # ... and do no attempt to scrape again
                del family.failed_paged[protein_table_page["url"]]

    # check if any pages for the family still have attempts left for trying for a successful scrape
    if len(list(family.failed_pages.keys())) == 0:
        retry_scrape = False
    else:
        retry_scrape = True

    return family, retry_scrape, failed_scrapes, format_failures


def get_pagination_pages(first_pagination_url, family, cazy_home, args):
    """Parse the first pagination page and retrieve URLs to all pagination page for the Family.

    :param first_pagination_url: str, URL to the fist page of the family
    :param family: Family class instance, represents a unique CAZy family
    :param cazy_home: str, URL to CAZy home page
    :param args: cmd-line args parser

    Return dict of error messages if errors arise OR list of URLs to pagination pages, and the
    total number of proteins in the family.
    """
    logger = logging.getLogger(__name__)

    # retrieve a list of all page urls of protein tables for the CAZy family
    first_pagination_page, error_message = crawler.get_page(
        first_pagination_url, args, max_tries=(args.retries + 1),
    )

    if first_pagination_page is None:
        logger.warning(
                f"Could not connect to {first_pagination_url} after {args.retries} attempts\n"
                f"The following error was raised:\n{error_message}\nTherefore, could not "
                "retrieve all pagination pages URLs, therefore, cannot scrape proteins from "
                f"{family.name}"
        )

        return(
            {
                "url": (
                    f"{first_pagination_url}\t{family.cazy_class}\t"
                    f"Failed to connect to first pagination page for {family.name}, therefore "
                    f"could not retrieve URLs to all pagination pages\t{error_message}"
                ),
                "format": None,
            },
            None,
        )

    # Get the URLS to all pages of proteins and the total number of proteins in the (sub)family
    protein_page_urls = get_pagination_page_urls(
        first_pagination_url, first_pagination_page, cazy_home, family.name,
    )

    return protein_page_urls


def get_pagination_page_urls(first_pagination_url, first_pagination_page, cazy_home, family_name):
    """Retrieve the URLs to all pages containing proteins for the current working family.

    Also retrieve the total number of proteins catagloued under the family.

    :param first_pagination_url: str, URL to first page contaiing proteins
    :param first_pagination_page: BS4 object, first page containing proteins
    :param cazy_home: str, URL of CAZy homepage
    :param family_name: str, name of the current working CAZy family

    Return list of URLs
    """
    protein_page_urls = [first_pagination_url]

    # retrieve the URL to the final page of protein records in the pagination listing
    try:
        last_pagination_url = first_pagination_page.find_all(
            "a", {"class": "lien_pagination", "rel": "nofollow"}
        )[-1]
    except IndexError:  # there is no pagination; a single-query entry
        last_pagination_url = None

    if last_pagination_url is not None:
        url_prefix = last_pagination_url["href"].split("PRINC=")[0] + "PRINC="
        last_princ_no = int(last_pagination_url["href"].split("PRINC=")[-1].split("#pagination")[0])
        url_suffix = "#pagination" + last_pagination_url["href"].split("#pagination")[-1]
        # Build list of urls to all pages in the pagination listing, increasing the PRINC increment
        protein_page_urls.extend(
            [f"{cazy_home}/{url_prefix}{_}{url_suffix}" for _ in range(
                1000, last_princ_no + 1000, 1000
            )]
        )

    return protein_page_urls


def get_html_page(url, family_name, args):
    """Get HTML page from CAZy and write to file on disk.

    :param url: str, URL to CAZy page
    :param family_name: str, name of the curernt working CAZy family
    :param args: cmd-line args parser

    Return dict of errors."""
    logger = logging.getLogger(__name__)
    # create file name from the URL
    filename = url
    for c in ['/', '.', '?', '#', ':']:
        filename = filename.replace(c, "_")
    if args.output is sys.stdout:
        outdir = os.getcwd()
    else:
        outdir = args.output
    filename = outdir / filename

    # connect to the page
    page, error = crawler.get_page(url, args, max_tries=args.retries)

    if page is None:
        logger.warning(
            f"Could not connection to {url} after {args.retries} attemtps\n"
            f"The following error was raised:\n{error}\nNo page written out to file"
        )
        return {"url": url, "error": error}

    # double check if page was fully returned, sometimes the protein table page is not populated
    cazyme_table = page.select("table")[1]

    if len(cazyme_table.select("tr")) == 1:
        # check if deleted family
        activities = page.select("table")[0].select("tr")[0].select("td")[0].contents

        if type(activities) is list:
            if activities[0] == 'Deleted family!':
                logger.warning(f"{family_name} is a deleted Family in CAZy. Downloaded page anyway")

            else:  # check if got any proteins
                data = page.find_all("div", {"class": "pos_choix"})
                try:
                    protein_total = int(re.findall(
                        r"all \(\d+\)", data[0].text, flags=re.IGNORECASE,
                    )[0].split("(")[1][:-1])
                except IndexError as e:
                    logger.warning(
                        f"Could not retrieve the number of proteins for {url}.\n"
                        "Therefore, could not verify page was properly returned. Saved file anyway."
                    )
                    protein_total = 1000

                if protein_total == 0:
                    logger.warning(
                        f"Family {url} is not listed as a deleted Family in CAZy but is an empty "
                        "family,\npage was written to file."
                    )

                else:
                    logger.warning(
                        f"A connection to {url} was made but the protein table was not populated "
                        "when it should be.\n Retrying connection to the page until fully retrieved"
                    )
                    report_dict = retry_failed_connections(url, family_name, args, filename)
                    return report_dict

        else:
            if activities == 'Deleted family!':
                logger.warning(f"{family_name} is a deleted Family in CAZy. Downloaded page anyway")

            else:
                data = page.find_all("div", {"class": "pos_choix"})
                try:
                    protein_total = int(re.findall(
                        r"all \(\d+\)", data[0].text, flags=re.IGNORECASE,
                    )[0].split("(")[1][:-1])
                except IndexError as e:
                    logger.warning(
                        f"Could not retrieve the number of proteins for {url}.\n"
                        "Therefore, could not verify page was properly returned. Saved file anyway."
                    )

                if protein_total == 0:
                    logger.warning(
                        f"Family {url} is not listed as a deleted Family in CAZy but is an empty "
                        "family,\npage was written to file."
                    )

                else:
                    logger.warning(
                        f"A connection to {url} was made but the protein table was not populated "
                        "when it should be.\n Retrying connection to the page until fully retrieved"
                    )
                    report_dict = retry_failed_connections(url, family_name, args, filename)
                    return report_dict

    with open(filename, "w") as fh:
        fh.write(str(page))
        time.sleep(5)  # adds delay to minimise burden on the CAZy server

    return {"url": None, "error": None}


def retry_failed_connections(url, family_name, args, filename):
    """Retry connection untill full CAZy page is returned.

    Sometimes HTML code is returned by the protein table is not populated when it should be. Keep
    retrying the connection until the full page is returned.

    :param url: str, URL to CAZy page
    :param family_name: str, name of the curernt working CAZy family
    :param args: cmd-line args parser
    :param filename: path to write out HTML file

    Return dict of errors.
    """
    page = None
    while page is None:
        page, error = crawler.get_page(url, args, max_tries=(args.retries + 1))

    # check if protein table page is present
    cazyme_table = page.select("table")[1]

    if len(cazyme_table.select("tr")) == 1:  # len of 1 is when the table is empty
        retry_failed_connections(url, family_name, args, filename)  # retry to scrape

    with open(filename, "w") as fh:
        fh.write(str(page))
        time.sleep(5)  # adds delay to minimise burden on the CAZy server

    return {"url": None, "error": None}


def parse_family_by_kingdom(family, cazy_home, args, kingdoms):
    """Retrieve webpages for specified Kingdoms.

    CAZy families have separate series of HTML tables for each taxonomy Kingdom, and a single series
    of HTML tables containing all proteins from all Kingdoms, called 'all'.

    :param family: Family class object, representation of CAZy family
    :param cazy_home: str, URL to CAZy home page
    :param args: cmd-line args parser
    :param kingdoms: set of taxonomy Kingdoms to retrieve webpages for

    Return:
        Family object,
        Boolean whether to scrape the family again,
        List of URLs which couldn't connect to CAZy (and no have re-try attempts left),
        List of pages whose parsing did match the expected format
    """
    logger = logging.getLogger(__name__)

    # define lists for storing error messages during scraping and parsing proteins
    failed_scrapes = []  # URLs of pages for which maximum number of scrape attempts is MET
    format_failures = []  # pages where no proteins were retrieved

    # check if there were pagination pages to which a connection could not be made previously

    if len(list(family.failed_pages.keys())) != 0:  # rescraping selected pages of the family
        protein_page_urls = list(family.failed_pages.keys())

        if len(protein_page_urls) == 1:
            first_url_template = family.url.replace(".html", "_all.html")
            if protein_page_urls[0] == first_url_template:
                protein_page_urls = get_pagination_pages_kingdom(
                    protein_page_urls[0], family, cazy_home, args, kingdoms,
                )

    else:  # scraping for the first time
        for kingdom in kingdoms:
            # compile URL to first family page of protein records
            first_pagination_url = family.url.replace(".html", f"_{kingdom}.html")

            # check url formating of first pagination url
            try:
                re.match(
                    rf"http://www.cazy.org/(\D\D|\D\D\D)(\d+_{kingdom}|\d+_\d+_{kingdom}).html",
                    first_pagination_url,
                ).group()
            except AttributeError:
                logger.warning(
                    "Incorrect formatting of first protein table page URL: "
                    f"{first_pagination_url}\nWill not try and connect to this URL."
                )
                format_failures.append(
                    f"{first_pagination_url}\tIncorrect URL format for the first protein table "
                    f"page for {family.name} in kingdom {kingdom}, could not scrape any CAZymes "
                    f"for this Kingdom for family {family.name}"
                )
                continue  # could not scrape any CAZymes for the Kingdom

        # retrieve pages to other pagination pages of protein tables for the family
        protein_page_urls = get_pagination_pages_kingdom(
            first_pagination_url, family, cazy_home, args, kingdom,
        )

        # iterate through the pages and retrieve proteins, adding to the local CAZy data

        for protein_table_page in tqdm(
            (y for x in (
                get_html_page(url, family.name, args) for url in protein_page_urls
            ) for y in x),
            total=len(protein_page_urls),
            desc=f"Retrieving Family pages for {family.name} {kingdom}",
        ):
            if protein_table_page["url"] is not None:  # Protein not retrieved, couldn't connect to CAZy
                # protein["url"][0] = URL
                # protein["url"][0] = Error message
                try:
                    family.failed_pages[protein_table_page["url"]] += 1
                except KeyError:
                    family.failed_pages[protein_table_page["url"]] = 1  # First failed attempt

                if family.failed_pages[protein_table_page["url"]] == args.retries:
                    # maximum attempts to connect have been reached no more attempts made, write to file
                    failed_scrapes.append(
                        f"{protein_table_page['url']}\t{family.cazy_class}\t"
                        f"Failed to connect to this page of proteins for {family.name}, "
                        f"and raised the following error message:\n{protein_table_page['error']}"
                    )
                    # ... and do no attempt to scrape again
                    del family.failed_paged[protein_table_page["url"]]

    # check if any pages for the family still have attempts left for trying for a successful scrape
    if len(list(family.failed_pages.keys())) == 0:
        retry_scrape = False
    else:
        retry_scrape = True

    return family, retry_scrape, failed_scrapes, format_failures


def get_pagination_pages_kingdom(first_pagination_url, family, cazy_home, args, kingdom):
    """Parse the first pagination page and retrieve URLs to all pagination page for the Family.

    :param first_pagination_url: str, URL to the fist page of the family
    :param family: Family class instance, represents a unique CAZy family
    :param cazy_home: str, URL to CAZy home page
    :param args: cmd-line args parser

    Return dict of error messages if errors arise OR list of URLs to pagination pages, and the
    total number of proteins in the family.
    """
    logger = logging.getLogger(__name__)

    # retrieve a list of all page urls of protein tables for the CAZy family
    first_pagination_page, error_message = crawler.get_page(
        first_pagination_url, args, max_tries=(args.retries + 1),
    )

    if first_pagination_page is None:
        logger.warning(
                f"Could not connect to {first_pagination_url} after {args.retries} attempts\n"
                f"The following error was raised:\n{error_message}\nTherefore, could not "
                "retrieve all pagination pages URLs, therefore, cannot scrape proteins from "
                f"{family.name} {kingdom}"
        )
        return(
            {
                "url": (
                    f"{first_pagination_url}\t{family.cazy_class} {kingdom} {family.name}\t"
                    f"Failed to connect to first pagination page for {family.name}, therefore "
                    f"could not retrieve URLs to all pagination pages\t{error_message}"
                ),
                "format": None,
            },
        )

    # Get the URLS to all pages of proteins and the total number of proteins in the (sub)family
    protein_page_urls = get_tax_page_urls(
        first_pagination_url, first_pagination_page, cazy_home, family.name,
    )

    return protein_page_urls


def get_tax_page_urls(first_pagination_url, first_pagination_page, cazy_home, family_name):
    """Retrieve the URLs to all pages containing proteins for the current working family.

    Also retrieve the total number of proteins catagloued under the family.

    :param first_pagination_url: str, URL to first page contaiing proteins
    :param first_pagination_page: BS4 object, first page containing proteins
    :param cazy_home: str, URL of CAZy homepage
    :param family_name: str, name of the current working CAZy family

    Return list of URLs
    """
    protein_page_urls = [first_pagination_url]

    try:  # retrieve the URL to the final page of protein records in the pagination listing
        last_pagination_url = first_pagination_page.find_all(
            "a", {"class": "lien_pagination", "rel": "nofollow"}
        )[-1]
    except IndexError:  # there is no pagination; there's only one page of proteins
        last_pagination_url = None

    if last_pagination_url is not None:
        url_prefix = last_pagination_url["href"].split("TAXO=")[0] + "TAXO="
        last_url_num = int(last_pagination_url["href"].split("TAXO=")[-1].split("#pagination")[0])
        url_suffix = "#pagination" + last_pagination_url["href"].split("#pagination")[-1]

        # Build list of urls to all pages in the pagination listing, increasing the PRINC increment
        protein_page_urls.extend(
            [f"{cazy_home}/{url_prefix}{_}{url_suffix}" for _ in range(
                100, last_url_num + 100, 100
            )]
        )

    return protein_page_urls
