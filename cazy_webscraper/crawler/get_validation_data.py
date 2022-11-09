#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Retrieve CAZy family populations from CAZy website to check the expected number of gamily members
were added to the local CAZyme database while scraping CAZy."""


import logging
import re
import time

from urllib.error import HTTPError

from tqdm import tqdm
from requests.exceptions import ConnectionError, MissingSchema
from saintBioutils.utilities import file_io
from saintBioutils.utilities.file_io import make_output_directory
from urllib3.exceptions import HTTPError, RequestError

import mechanicalsoup


class CazyClass:
    """A single CAZy class.

    Used to keep track of specific families that need to be scraped again.
    """

    def __init__(self, name, url, tries, failed_families=None):
        self.name = name
        self.url = url
        self.tries = tries  # number of attempts to be scraped
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


def get_validation_data(
    cazy_home_url,
    excluded_classes,
    cazy_synonym_dict,
    config_dict,
    cache_dir,
    connection_failures_logger,
    time_stamp,
    args,
):
    """Coordinate retrieving the population sizes of CAZy familes from the CAZy website.

    :param cazy_home_url: str, URL to CAZy home page
    :param excluded_classes: list of CAZy classes NOT to scrape
    :param cazy_synonym_dict: dict of accepted CAZy class name synonyms
    :param config_dict: dict keyed by CAZy classes, values by set of CAZy families to scrape
    :param cache_dir: path to cache dir
    :param connection_failures_logger: logger, logg incorrect URLs and URLs to which a connection 
            could not be made
    :param time_stamp: str, time cazy_webscraper was invoked
    :param args: cmd-line args parser

    Return dict, keyed by CAZy family (str) and valued by population size (int)
    """
    # make dir fo caching HTML files
    cache_dir = cache_dir / "html"
    make_output_directory(cache_dir, args.force, args.nodelete_cache)

    cazy_fam_populations = {}  # {fam(str): population(int)}

    # retrieve list of CAZy class instances, one instance per class to be scrapped
    cazy_classes = get_cazy_classes(
        cazy_home_url,
        excluded_classes,
        cazy_synonym_dict,
        cache_dir,
        time_stamp,
        args,
    )
    if cazy_classes is None:
        return

    for cazy_class in tqdm(cazy_classes, desc="Retrieving CAZy family population sizes"):

        # first attempt of scraping, retrieve URLs to CAZy families
        if len(list(cazy_class.failed_families.keys())) == 0:
            fam_pops_to_retrieve = config_dict[cazy_class.name]  # retrieve user specified fams
        else: 
            fam_pops_to_retrieve = list(cazy_class.failed_families.keys())  # retry failed connections

        family_populations, err_message, incorrect_urls, failed_families = get_cazy_family_pops(
            cazy_class.name,
            cazy_class.url,
            cazy_home_url,
            fam_pops_to_retrieve,
            cache_dir,
            time_stamp,
            args,
        )

        if incorrect_urls is not None:  # log families for which compiled URL is incorrect
            [connection_failures_logger.warning(url_message) for url_message in incorrect_urls]

        if family_populations is None:  # couldn't retrieve family populations
            cazy_class.tries += 1

            # check if maximum number of attempts to connect have been met
            if cazy_class.tries == (args.retries + 1):  # Maximum number of tries met
                connection_failures_logger.warning(
                    f"{cazy_class.url}\t"
                    f"{cazy_class.name}\t"
                    f"CAZy family populations not retrieved from {cazy_class.name}\t"
                    f"{err_message}"
                )

            else:
                for fam in failed_families:
                    try:
                        cazy_class.failed_families[fam] += 1
                        if cazy_class.failed_families[fam] == (args.retries + 1):
                            # max number of attemptes made, do not retry connection
                            del cazy_class.failed_families[fam]
                    except KeyError:
                        cazy_class.failed_families[fam] = 1

                cazy_classes.append(cazy_class)  # retry retriving family populations later

            continue  # go onto next CAZy class

        else:  # retrieved CAZy family populations
            cazy_fam_populations.update(family_populations)

    # log any errors that meant no family population could be retrieved
    for cazy_class in cazy_classes:
        for fam in list(cazy_class.failed_families.keys()):
            connection_failures_logger.warning(
                f"{fam}\t"
                "Retrieved no family population for data retrieval validation\n"
                f"Failed to conencted to CAZy after {(args.retries + 1)*(args.retries +1)} attempts"
            )

    return cazy_fam_populations


def get_cazy_classes(
    cazy_home_url,
    excluded_classes,
    cazy_synonym_dict,
    cache_dir,
    time_stamp,
    args,
    unit_test=False,
):
    """Returns a list of CAZy class instances.

    :param cazy_url: str, URL to the CAZy home page.
    :param excluded_classes: list, list of CAZy classes not to be scraped
    :param cazy_synonym_dict: dictionary of offical CAZy class names
    :param cache_dir: path to cache dir
    :param time_stamp: str, time cazy_webscraper was invoked
    :param args: cmd line args parser

    Return list of CazyClass instances, or None and an error message.
    """
    logger = logging.getLogger(__name__)
    logger.info("Retrieving URLs to summary CAZy class pages")

    # define items to be excluded from returned class list, ALWAYS exlide links to genomes
    if excluded_classes is not None:
        exclusions = tuple(excluded_classes)
    else:
        exclusions = tuple()

    homepage, error = get_page(cazy_home_url, args, max_tries=(args.retries + 1))

    if homepage is None:
        logger.error(
            (
                f"Failed to connect to CAZy home page after {args.retries} attempts.\n"
                "The following error was raised:\n"
                f"{error}"
                "Could not retrieve URLs to CAZy classes.\n"
                "Check the network connection.\n"
                "Terminating program."
            )
        )
        return

    cache_name = cazy_home_url.replace('.', '_')
    cache_path = cache_dir / f"{cache_name}_{time_stamp}.html"
    if unit_test is False:
        with open(cache_path, "w") as cache:
            cache.write(homepage)

    # retrieve the h3 elements with class spip
    h3_spip_elements = homepage.find_all("h3", {"class": "spip"})

    # retrieve the div section containing the h3 element for Enzyme classes catalgoued by CAZy
    try:
        enzyme_classes_div = [
            _ for _ in h3_spip_elements if (
                str(_.contents[0].strip()).replace(u'\xa0', ' ')
                ) == 'Enzyme Classes currently covered'][0].parent

        # Retreive the enzyme class page URLs suffixs
        enzyme_class_urls = [
            f"{cazy_home_url}/{_['href']}" for _ in enzyme_classes_div.find_all("a") 
            if (not _["href"].startswith("http"))
            and (str(_.contents[0]) not in exclusions)
        ]

        # retrieve the div section containing the h3 element for Associated Module catalgoued by CAZy
        associated_module_div = [
            _ for _ in h3_spip_elements if (
                str(_.contents[0].strip()).replace(u'\xa0', ' ')
                ) == 'Associated Modules currently covered'][0].parent

        # Retreive the enzyme class page URLs suffixs
        associated_module_urls = [
            f"{cazy_home_url}/{_['href']}" for _ in associated_module_div.find_all("a") 
            if (not _["href"].startswith("http"))
            and (str(_.contents[0]) not in exclusions)
        ]

    except (AttributeError, IndexError) as err:
        logger.error(
            (
                "Error raised during retrieving of CAZy class URLs.\n"
                "Therefore, cannot validate data retrieval. \n"
                "Will proceed with scraping CAZy. Error message:\n"
            ),
            exc_info=1,
        )
        return

    # compile the full CAZy class URLs from the homepage url and class suffixes

    if len(enzyme_class_urls) == 0 and len(associated_module_urls) == 0:
        logger.error(
            (
                "Failed retrieve URLs to CAZy classes from the CAZy homepage.\n"
                "Therefore, cannot validate data retrieval. \n"
                "Will proceed with scraping CAZy"
            ),
            exc_info=1,
        )
        return

    # create CAZyClass objects
    cazy_class_urls = enzyme_class_urls + associated_module_urls
    cazy_classes = []

    for url in cazy_class_urls:
        # retrieve class name and standardise it
        class_name = url[20:-5]
        for key in cazy_synonym_dict:
            if class_name in cazy_synonym_dict[key]:
                class_name = key

        cazy_class = CazyClass(class_name, url, 0)
        cazy_classes.append(cazy_class)

    logger.info(
        "Retrieved URLs for:"
        f"{len(enzyme_class_urls)} Enzyme Classes and\n"
        f"{len(associated_module_urls)} Associated Modules classes"
    )

    return cazy_classes


def get_cazy_family_pops(
    class_name,
    class_url,
    cazy_home_url,
    fam_pops_to_retrieve,
    cache_dir,
    time_stamp,
    args,
    unit_test=False,
):
    """Retrieve all protein members of each CAZy family within the given CAZy class.

    :param class_name: str, name of CAZy class
    :param class_url: str, URL to CAZy class webpage
    :param cazy_home_url: str, URL to CAZy home page
    :param fam_pops_to_retrieve: list of CAZy families to scrape
    :param cache_dir: str representing Path to dir to write out downloaded family file to
    :param time_stamp: str, date and time cazy_webscraper was invoked
    :param args: args parser object

    Returns:
    A dict of CAZy family populations (fam: pop)
    An error message from when retrieving CAZy family URLs
    A list of incorrectly formated URLs
    A list of URLs of families from which a connection could not be made
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Retrieving URLs to families under {class_name}")

    failed_connections = []
    incorrect_urls = []
    family_populations = {}

    # get the html code of the class page
    class_page, error = get_page(class_url, args, max_tries=(args.retries + 1))

    if class_page is None:
        logger.error(
                f"Couldn't connect to {class_url} after {(args.retries + 1)} attempts.\n"
                f"The following error was raised:\n{error}"
        )
        return None, error, incorrect_urls, failed_connections

    cache_name = class_url.replace('.', '_')
    cache_path = cache_dir / f"{cache_name}_{time_stamp}.html"

    if unit_test is False:
        with open(cache_path, "w") as cache:
            cache.write(class_page)

    family_urls, url_err_message, incorrect_urls = get_families_urls(cazy_home_url, class_page, args)

    if family_urls is None:
        return None, url_err_message, incorrect_urls, failed_connections

    for fam_url in tqdm(family_urls, desc=f"Retrieing fam populations for {class_name}"):
        fam_name = fam_url.replace(cazy_home_url, "").split(".")[0]

        if (fam_pops_to_retrieve is not None) and (fam_name not in fam_pops_to_retrieve):
            continue

        family_page, err = get_page(fam_url, args, max_tries=(args.retries + 1))
        if err is not None:
            logger.warning(
                f"Failed to connect to {fam_name} webpage at {fam_url}\n"
                f"to retrieve the population size, after trying {args.retries + 1} times\n"
            )
            failed_connections.append(fam_url)
            continue

        cache_name = fam_url.replace('.', '_')
        cache_path = cache_dir / f"{cache_name}_{time_stamp}.html"
        if unit_test is False:
            with open(cache_path, "w") as cache:
                cache.write(family_page)

        # retrieve the table containing the Family data
        family_data = family_page.find_all("div", {"class": "pos_choix"})

        try:
            fam_pop = int(re.findall(
                r"Download \D{2,3}\d+? \(\d+?\)",
                family_data[0].text,
                flags=re.IGNORECASE,
            )[0].split("(")[1].replace(")", ""))

        except (IndexError, AttributeError) as err:
            fam_pop = 0

        if fam_pop == 0:  # check if an empty or deleted fam
            try:
                family_activities_cell = family_page.select("table")[
                    0].select("tr")[0].select("td")[0].contents[0].strip()

                if family_activities_cell == 'Deleted family!':
                    logger.warning(f"{fam_name} is a deleted family in CAZy")
                else:
                    logger.warning(f"{fam_name} is an empty CAZy family")
            except Exception:
                logger.warning(f"Could not retrieve family population for {fam_name}")
                fam_pop = 'Failed Retrieval'

            logger.warning(
                f"{fam_name}\t"
                f"{fam_url}\t"
                f"Failed to retrieve population for {fam_name}\t"
                f"{err}"
            )

        family_populations[fam_name] = fam_pop        # handle errors

    logger.info(f"Retrieved URLs for {len(family_urls)} from {class_name} class page")

    return family_urls, url_err_message, incorrect_urls, failed_connections


def get_families_urls(cazy_home_url, class_name, class_page, args):
    """Retrieve the URLs to CAZy family pages.

    :param cazy_home_url: str, CAZ home page URL
    :param class_name: str, name of CAZy class
    :param class_page: bs4 soup object, CAZy class summary page
    :param args: cmd-line args parser

    Return:
    List of CAZy family URLs
    Str, message if any errors arose
    List of inccorectly formated CAZy family URLs
    """
    logger = logging.getLogger(__name__)
    incorrect_urls = []
    err_message = None

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

    family_urls = [f"{cazy_home_url}/{_['href']}" for _ in tables[0].find_all("a")]
    try:
        family_urls.append(f"{cazy_home_url}/{tables[1].a['href']}")
    except TypeError:
        family_urls = None

    if (args.subfamilies is False) and (family_urls is None):
        err_message = f"Failed to retrieve URLs to CAZy families for {class_name}"
        logger.warning(err_message)
        return None, err_message, incorrect_urls

    # retrieve URLs to subfamilies
    if args.subfamilies is True:
        subfam_urls = get_subfamily_links(family_h3_element, cazy_home_url)

        if (family_urls is None) and (subfam_urls is None):
            err_message = f"Failed to retrieve URLs to CAZy subfamilies for {class_name}"
            logger.warning(err_message)
            return None, err_message, incorrect_urls

        elif family_urls is None:
            family_urls = subfam_urls
            err_message = (
                f"Failed to retrieve URLs to CAZy families for {class_name}\n"
                f"But successfully retrieved the URLs to the CAZy subfamilies for {class_name}"
            )
            logger.warning(err_message)

        else:
            family_urls += subfam_urls

    return family_urls, err_message, incorrect_urls


def get_subfamily_links(family_h3_element, cazy_home_url):
    """Retrieve URL links to CAZy subfamilies.

    :param family_h3_element: bs4.element.Tag, h3 element titling the page div
    :param cazy_home_url: str, URL to CAZy homepage

    Return list of URLs to subfamilies.
    """
    parent_div = family_h3_element.parent
    all_links = parent_div.find_all("a")

    pattern = re.compile(r"\D+?\d+?_\d+?\.html")

    urls = []  # empty list to store subfamily URLs

    for link in all_links:
        try:
            search_result = re.search(pattern, link["href"])
            urls.append(f"{cazy_home_url}/{search_result.group()}")
        except (KeyError, AttributeError) as error:
            # KeyError raised if link does not have ['href']
            # AttributeError error raised if search_result is None becuase not subfam link
            pass

    if len(urls) == 0:
        return
    else:
        return urls


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
