#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Module for crawling through the CAZy website and parsing HTML from the CAZy website."""


import logging
import re
import sys
import time

from tqdm import tqdm
from requests.exceptions import ConnectionError, MissingSchema
from urllib3.exceptions import HTTPError, RequestError

import mechanicalsoup

from scraper.sql import sql_interface


class CazyClass:
    """A single CAZy class.

    Used to keep track of if specific families need to be scraped again.
    """

    def __init__(self, name, url, tries, failed_families=None):
        self.name = name
        self.url = url
        self.tries = tries
        if failed_families is None:
            self.failed_families = {}  # keyed by Family instance, valued by attempted scrapes (int)
        else:
            self.failed_families = failed_families

    def __str__(self):
        return f"<CAZy class: {self.name} id={id(self)}>"

    def __repr__(self):
        return(
            f"<CAZy class: {self.name} id={id(self)} url={self.url} "
            f"attempted connections={self.tries}>"
        )


class Family:
    """A single CAZy family.

    Used to keep track if family needs to be scraped again.
    """

    def __init__(self, name, cazy_class, url, failed_pages=None):
        self.name = name
        self.cazy_class = cazy_class
        self.url = url
        if failed_pages is None:
            # {kingdom: {paginiation_page_url: number of tries}}
            self.failed_pages = {}
        else:
            self.failed_pages = failed_pages

    def __str__(self):
        return f"CAZy family {self.name}"

    def __repr__(self):
        return f"<Family: {id(self)}: {self.name}>"


def get_cazy_classes(cazy_home, excluded_classes, cazy_dict, args):
    """Returns a list of CAZy class main/home page URLs for each specified class as the CAZy site.

    :param cazy_url: str, URL to the CAZy home page.
    :param excluded_classes: list, list of CAZy classes not to be scraped
    :param cazy_dict: dictionary of offical CAZy class names
    :param args: cmd line args parser

    Return list of CazyClass instances, or None and an error message.
    """
    logger = logging.getLogger(__name__)
    logger.info("Retrieving URLs to summary CAZy class pages")

    # define items to be excluded from returned class list, ALWAYS exlide links to genomes
    if excluded_classes is not None:
        exclusions = tuple(["<strong>Genomes</strong>"] + excluded_classes)
    else:
        exclusions = "<strong>Genomes</strong>"

    home_page = None
    tries = 0  # number of the attempted connections to CAZy

    while (home_page is None) and (tries < args.retries + 1):
        home_page, error = get_page(cazy_home, args, max_tries=(args.retries + 1))

        if (home_page is None) and (tries < (args.retries + 1)):
            logger.error(
                f"Failed to connect to CAZy homepage after 10 attempts,\n"
                f"On attempt# {(tries+1)}/{args.retries + 1} the following error was raised:\n"
                f"{error}\n"
                f"Reattempting for attempt# {(tries+2)} in 10s."
            )
            time.sleep(10)
            tries += 1

    if home_page is None:
        logger.error(
            (
                "Failed to connect to CAZy home-page after multiple attempts.\n"
                "The following error was raised:\n"
                f"{error}"
                "Could not retrieve URLs to CAZy classes.\n"
                "Check the network connection.\n"
                "Terminating program."
            )
        )
        sys.exit(1)

    try:
        class_urls = [
            f"{cazy_home}/{_['href']}"
            for _ in home_page.find_all("a", {"class": "spip_out"})
            if (not _["href"].startswith("http")) and (str(_.contents[0]) not in exclusions)
        ]
    except AttributeError:  # raise if can't find results with find_all("a", {"class": "spip_out"})
        logger.error(
            (
                "Failed retrieve URLs to CAZy classes from the CAZy homepage.\n"
                "Therefore, cannot scrape CAZy classes, or families\n"
                "Terminating program."
            ),
            exc_info=1,
        )
        sys.exit(1)

    if len(class_urls) == 0:
        logger.error(
            (
                "Failed retrieve URLs to CAZy classes from the CAZy homepage.\n"
                "Therefore, cannot scrape CAZy classes, or families\n"
                "Terminating program."
            ),
            exc_info=1,
        )
        sys.exit(1)

    # create CAZyClass objects
    cazy_classes = []

    for url in class_urls:
        # retrieve class name and standardise it
        class_name = url[20:-5]
        for key in cazy_dict:
            if class_name in cazy_dict[key]:
                class_name = key

        cazy_class = CazyClass(class_name, url, 0)
        cazy_classes.append(cazy_class)

    return cazy_classes


def get_cazy_family_urls(class_url, class_name, cazy_home, args):
    """Retrieve all protein members of each CAZy family within the given CAZy class.

    :param class_url: str, URL to the CAZy class
    :param class_name: str, name of CAZy class
    :param cazy_home: str, URL to CAZy home page
    :param args: args parser object

    Returns list Family class objects, error message when connecting to CAZy and list of incorrectly
    formated URLs
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Retrieving URLs to families under {class_name}")

    # scrape the class page
    class_page, error = get_page(class_url, args, max_tries=(args.retries + 1))

    if class_page is None:
        logger.error(
                f"Couldn't connect to {class_url} after {(args.retries + 1)} attempts. "
                f"The following error was raised:\n{error}"
        )
        return None, error, None

    # Retrieve URLs to the CAZy family pages

    # retrieve the <h3> element that titles the div section containing the tables of family links
    family_h3_element = [
        _
        for _ in class_page.find_all("h3", {"class": "spip"})
        if str(_.contents[0]).strip() == "Tables for Direct Access"
    ][0]

    # retrieve all tables within the parent div section of the <h3> element
    tables = family_h3_element.parent.find_all("table")

    # tables[0] is the table containing links to CAZy families
    # tables[1] is the table containing the link to unclassified proteins

    family_urls = [f"{cazy_home}/{_['href']}" for _ in tables[0].find_all("a")]
    try:
        family_urls.append(f"{cazy_home}/{tables[1].a['href']}")
    except TypeError:
        family_urls = None

    if (args.subfamilies is False) and (family_urls is None):
        logger.warning(f"Failed to retrieve URLs to CAZy families for {class_name}\n")
        return None, f"Failed to retrieve URLs to CAZy families for {class_name}\n", None

    # retrieve URLs to subfamilies
    if args.subfamilies is True:
        subfam_urls = get_subfamily_links(family_h3_element, cazy_home)

        if (family_urls is None) and (subfam_urls is None):
            logger.warning(f"Failed to retrieve URLs to CAZy subfamilies for {class_name}")
            return(
                None,
                f"Could not retrieve family and subfamily URLs from class page of {class_name}",
                None
            )

        elif family_urls is None:
            family_urls = subfam_urls
            logger.warning(
                f"Failed to retrieve URLs to CAZy families for {class_name}\n"
                f"But successfully retrieved the URLs to the CAZy subfamilies for {class_name}"
            )

        else:
            family_urls += subfam_urls

    # create Family class objects
    cazy_families = []
    incorrect_urls = []

    for url in family_urls:
        # check URL format
        try:
            re.match(
                r"http://www.cazy.org/(\D{2,3})(\d+|\d+_\d+).html", url
            ).group()
        except AttributeError as error:
            logger.warning(
                f"Format of URL {url} is incorrect from {class_name}.\n"
                "Will not attempt to scrape this URL."
            )
            incorrect_urls.append(
                f"{url}\t"
                f"{class_name}\t"
                "Format of the URL is incorrect\t"
                f"{error}"
            )
            continue

        family_name = url[(len(cazy_home) + 1): -5]

        family = Family(family_name, class_name, url)
        family.members = set()  # later used to store Protein members
        cazy_families.append(family)

    if len(incorrect_urls) == 0:
        incorrect_urls = None

    return cazy_families, None, incorrect_urls


def get_subfamily_links(family_h3_element, cazy_home):
    """Retrieve URL links to CAZy subfamilies.

    :param family_h3_element: bs4.element.Tag, h3 element titling the page div
    :param cazy_home: str, URL to CAZy home_page

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


def row_to_protein(row, family_name, taxonomy_filters, kingdom, ec_filters, session, args):
    """Parses row and collates data to be added to records in the local CAZy database.

    Each row, in order, contains the protein name, EC number, source organism, GenBank ID(s),
    UniProt ID(s), and PDB accession(s).

    :param row: tr element from CAZy family protein page
    :param family_name: str, name of CAZy family
    :param taxonomy_filters: set of genera, species and strains to restrict the search to
    :param kingdom: str, Taxonomy Kingdom of CAZymes currently being scraped
    :param ec_filters: set of EC filters to limit scrape to
    :param session: open sqlalchemy database session

    Returns a dictionary to store and reflect any errors that may have occurred during the parsing
    of the proteins.

    Return dictionary, and session
    """
    logger = logging.getLogger(__name__)

    # retrieve list of cells ('td' elements) in row
    tds = list(row.find_all("td"))

    protein_name = tds[0].contents[0].strip()

    source_organism = tds[2].a.get_text()

    if taxonomy_filters is not None:  # apply taxonomy filter
        if any(filter in source_organism for filter in taxonomy_filters) is False:
            # CAZyme does not match any taxonomy filters
            return {"url": None, "error": None, "sql": None, "format": None}

    ec_numbers = []
    all_links = None
    try:
        all_links = tds[1].find_all("a")
    except TypeError:  # raised when no EC numbers are listed
        pass

    if all_links is not None:
        for link in all_links:
            ec_numbers.append(link.text)
    else:
        ec_numbers = None

    if len(ec_filters) != 0:  # user had defined EC numbers to limit the scrape to, apply EC# filter
        if ec_numbers is None:
            return {"url": None, "error": None, "sql": None, "format": None}
        scrape = False
        for ec in ec_numbers:
            if ec in ec_filters:
                scrape = True
                break
        if scrape is False:
            return {"url": None, "error": None, "sql": None, "format": None}

    # retrieve the BeautifulSoup elements of the cell containing accessions for the respective db
    gbk_bs_elements = [_ for _ in row.select("td")[3].contents if getattr(_, "name", None) != "br"]
    uni_bs_elements = [_ for _ in row.select("td")[4].contents if getattr(_, "name", None) != "br"]
    pdb_bs_elements = [_ for _ in row.select("td")[5].contents if getattr(_, "name", None) != "br"]

    # Retrieve primary GenBank and UniProt accessions (identified by being written in bold)
    # CAZy defines the 'best' GenBank and UniProt models by writting them in bold
    gbk_primary = [_.text for _ in row.select("td")[3].find_all('b')]
    uni_primary = [_.text for _ in row.select("td")[4].find_all('b')]

    # Retrieve all accessions listed for each database
    # At this stage the non-primary accession lists containg ALL accessions in the cell
    gbk_nonprimary = get_all_accessions(gbk_bs_elements)
    uni_nonprimary = get_all_accessions(uni_bs_elements)
    pdb_accessions = get_all_accessions(pdb_bs_elements)

    # create dict for storing error messages for writing to the failed_to_scrape output file
    report_dict = {"url": None, "error": None, "sql": None}

    # Ensure a single primary GenBank accession is retrieved
    if len(gbk_primary) == 0:

        if len(gbk_nonprimary) == 0:
            warning = (
                f"NO GenBank accessions retrieved for {protein_name} in {family_name}.\n"
                "Adding protein with the GenBank accession: 'NA' as a new protein to the db."
            )
            logger.warning(warning)
            gbk_primary = ["NA"]
            report_dict["error"] = warning
            report_dict["sql"] = protein_name

        else:
            warning = (
                f"GenBank accessions retrieved for {protein_name} in {family_name} but none were "
                f"written as primary in CAZy.\nThe first accession {gbk_nonprimary[0]} written "
                "as primary in the db "
            )
            gbk_primary.append(gbk_nonprimary[0])
            gbk_nonprimary.remove(gbk_nonprimary[0])

    elif len(gbk_primary) == 1:
        # Remove the primary accession from the non-primary accession list
        if gbk_primary[0] in gbk_nonprimary:
            gbk_nonprimary.remove(gbk_primary[0])

    else:
        warning = (
            f"Multiple primary GenBank acccessions retrieved for {protein_name} in "
            f"{family_name}.\nOnly the first listed accession will be written as primary, the "
            "others are written as non-primary. These are:"
        )
        # remove all but first listed primary accession
        for accession in gbk_primary[1:]:
            warning += f"\nGenBank accession: {accession}"
            gbk_primary.remove(accession)
        # remove primary accession from the non-primary accession list
        if gbk_primary[0] in gbk_nonprimary:
            gbk_nonprimary.remove(gbk_primary[0])

        logger.warning(warning)
        report_dict["error"] = warning
        report_dict["sql"] = protein_name

    if (len(uni_primary) == 0) and (len(uni_nonprimary) != 0):
        # move the first listed UniProt accession to the primary list
        uni_primary.append(uni_nonprimary[0])
        uni_nonprimary.remove(uni_primary[0])

    elif len(uni_primary) > 1:
        warning = (
            f"Multiple UniProt primary accessions retrieved for {protein_name} in "
            f"{family_name}.\nAll listed as primary UniProt accessions in the local db, including:"
        )
        for accession in uni_primary:
            warning += f"\nUniProt accession: {accession}"
        logger.warning(warning)
        report_dict["error"] = warning
        report_dict["sql"] = protein_name

    # Remove primary UniProt accessions from the non-primary accessions list
    for accession in uni_primary:
        if accession in uni_nonprimary:
            uni_nonprimary.remove(accession)

    # add protein to database
    try:
        sql_interface.add_protein_to_db(
            protein_name,
            family_name,
            source_organism,
            kingdom,
            gbk_primary[0],
            session,
            args,
            ec_numbers,
            gbk_nonprimary,
            uni_primary,
            uni_nonprimary,
            pdb_accessions,
        )

    except Exception as error_message:
        warning = (
            f"Failed to add {protein_name} to SQL database, "
            f"the following error was raised:\n{error_message}"
        )
        logger.warning(warning, exc_info=1)
        report_dict["error"] = warning
        report_dict["sql"] = protein_name

    return report_dict


def row_to_protein_in_dict(
    row,
    family_name,
    taxonomy_filters,
    ec_filters,
    session,
):
    """Parses row and collates data to be added to records to dictionary of CAZymes.

    Each row, in order, contains the protein name, EC number, source organism, GenBank ID(s),
    UniProt ID(s), and PDB accession(s).

    :param row: tr element from CAZy family protein page
    :param family_name: str, name of CAZy family
    :param taxonomy_filters: set of genera, species and strains to restrict the search to
    :param kingdom: str, Taxonomy Kingdom of CAZymes currently being scraped
    :param ec_filters: set of EC filters to limit scrape to
    :param session: open sqlalchemy database session

    Returns a dictionary to store and reflect any errors that may have occurred during the parsing
    of the proteins.

    Return dictionary.
    """
    logger = logging.getLogger(__name__)

    # retrieve list of cells ('td' elements) in row
    tds = list(row.find_all("td"))

    protein_name = tds[0].contents[0].strip()

    source_organism = tds[2].a.get_text()

    if taxonomy_filters is not None:  # apply taxonomy filter
        if any(filter in source_organism for filter in taxonomy_filters) is False:
            # CAZyme does not match any taxonomy filters
            return {"url": None, "error": None, "sql": None}

    ec_numbers = []
    all_links = None
    try:
        all_links = tds[1].find_all("a")
    except TypeError:  # raised when no EC numbers are listed
        pass

    if all_links is not None:
        for link in all_links:
            ec_numbers.append(link.text)
    else:
        ec_numbers = None

    if len(ec_filters) != 0:  # user had defined EC numbers to limit the scrape to, apply EC# filter
        if ec_numbers is None:
            return {"url": None, "error": None, "sql": None}
        scrape = False
        for ec in ec_numbers:
            if ec in ec_filters:
                scrape = True
                break
        if scrape is False:
            return {"url": None, "error": None, "sql": None}

    # retrieve the BeautifulSoup elements of the cell containing accessions for the respective db
    gbk_bs_elements = [_ for _ in row.select("td")[3].contents if getattr(_, "name", None) != "br"]

    # Retrieve primary GenBank and UniProt accessions (identified by being written in bold)
    # CAZy defines the 'best' GenBank and UniProt models by writting them in bold
    gbk_primary = [_.text for _ in row.select("td")[3].find_all('b')]

    # Retrieve all accessions listed for each database
    # At this stage the non-primary accession lists containg ALL accessions in the cell
    gbk_nonprimary = get_all_accessions(gbk_bs_elements)

    # create dict for storing error messages for writing to the failed_to_scrape output file
    report_dict = {"url": None, "error": None, "sql": None}

    # Ensure a single primary GenBank accession is retrieved
    if len(gbk_primary) == 0:

        if len(gbk_nonprimary) == 0:
            warning = (
                f"NO GenBank accessions retrieved for {protein_name} in {family_name}.\n"
                "Adding protein with the GenBank accession: 'NA' as a new protein to the db."
            )
            logger.warning(warning)
            gbk_primary = ["NA"]
            report_dict["error"] = warning
            report_dict["protein"] = protein_name

        else:
            warning = (
                f"GenBank accessions retrieved for {protein_name} in {family_name} but none were "
                f"written as primary in CAZy.\nThe first accession {gbk_nonprimary[0]} written "
                "as primary in the db "
            )
            gbk_primary.append(gbk_nonprimary[0])
            gbk_nonprimary.remove(gbk_nonprimary[0])

    elif len(gbk_primary) == 1:
        # Remove the primary accession from the non-primary accession list
        if gbk_primary[0] in gbk_nonprimary:
            gbk_nonprimary.remove(gbk_primary[0])

    else:
        warning = (
            f"Multiple primary GenBank acccessions retrieved for {protein_name} in "
            f"{family_name}.\nOnly the first listed accession will be written as primary, the "
            "others are written as non-primary. These are:"
        )
        # remove all but first listed primary accession
        for accession in gbk_primary[1:]:
            warning += f"\nGenBank accession: {accession}"
            gbk_primary.remove(accession)
        # remove primary accession from the non-primary accession list
        if gbk_primary[0] in gbk_nonprimary:
            gbk_nonprimary.remove(gbk_primary[0])

        logger.warning(warning)
        report_dict["error"] = warning
        report_dict["protein"] = protein_name

    try:
        session[gbk_primary[0]]
        session[gbk_primary[0]].append(family_name)
    except KeyError:
        session[gbk_primary[0]] = [family_name]

    for acc in gbk_nonprimary:
        try:
            session[acc]
            session[acc].append(family_name)
        except KeyError:
            session[acc] = [family_name]

    return report_dict, session


def get_all_accessions(bs_element_lst):
    """Retrieve all accessions listed in a cell from a CAZyme table.

    :param bs_element_list: list of BeautifulSoup element from cell in HTML table

    Return list of accessions."""
    accessions = []

    for bs_element in bs_element_lst:
        try:
            if bs_element.name == "a":  # Hyperlinked, extract accession and add to primary
                accessions.append(bs_element.text)
            elif bs_element.strip() != "":  # There is text in element
                accessions.append(bs_element)
        except TypeError:
            pass

    return accessions


def browser_decorator(func):
    """Decorator to re-invoke the wrapped function up to 'args.retries' times."""

    def wrapper(*args, **kwargs):
        logger = logging.getLogger(__name__)
        tries, success, err = 0, False, None

        while not success and (tries < kwargs['max_tries']):
            try:
                response = func(*args, **kwargs)
            except (
                ConnectionError,
                HTTPError,
                OSError,
                MissingSchema,
                RequestError,
            ) as err_message:
                if (tries < kwargs['max_tries']):
                    logger.warning(
                        f"Failed to connect to CAZy on try {tries}/{kwargs['max_tries']}.\n"
                        f"Error: {err_message}"
                        "Retrying connection to CAZy in 10s"
                    )
                success = False
                response = None
                err = err_message
            if response is not None:  # response was successful
                success = True
            # if response from webpage was not successful
            tries += 1
            time.sleep(10)
        if (not success) or (response is None):
            logger.warning(f"Failed to connect to CAZy.\nError: {err}")
            return None, err
        else:
            return response, None

    return wrapper


@browser_decorator
def get_page(url, args, **kwargs):
    """Create browser and use browser to retrieve page for given URL.

    :param url: str, url to webpage
    :param args: cmd-line args parser
    :kwargs max_tries: max number of times connection to CAZy can be attempted

    Return browser response object (the page).
    """
    # create browser object
    browser = mechanicalsoup.Browser()
    # create response object
    page = browser.get(url, timeout=args.timeout)
    page = page.soup

    return page
