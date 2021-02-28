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
"""Module for crawling through the CAZy website and parsing HTML from the CAZy website.

:func get_cazy_family_urls: retrieve URLs to families on CAZy class summary page
:func get_subfamily_links: retrieve URLs to subfamilies on CAZy class summary page
:func parse_family: build Family class object to represent CAZy family
:func parse_family_pages: retrieve all URLs to pages containing protein records for CAZy family
:func parse_proteins: retrieve protein records from protein table page
:func row_to_protein: parse the protein record to build a Protein class object
"""


import logging
import re
import sys
import time

from tqdm import tqdm
from requests.exceptions import ConnectionError, MissingSchema
from urllib3.exceptions import HTTPError, RequestError

import bs4
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
            self.failed_pages = {}  # keyed by URL, valued by number of attempted scrapes
        else:
            self.failed_pages = failed_pages

    def __str__(self):
        return f"CAZy family {self.name}: {len(self.members)} protein members"

    def __repr__(self):
        return f"<Family: {id(self)}: {self.name}, {len(self.members)} protein members"


def get_cazy_classes(cazy_home, excluded_classes, max_tries, cazy_dict):
    """Returns a list of CAZy class main/home page URLs for each specified class as the CAZy site.

    :param cazy_url: str, URL to the CAZy home page.
    :param excluded_classes: list, list of CAZy classes not to be scraped
    :param max_tries: int, maximum number of times to try scrape if errors are encountered
    :param cazy_dict: dictionary of offical CAZy class names

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

    while (home_page is None) and (tries < max_tries):
        home_page, error = get_page(cazy_home)

        if (home_page is None) and (tries < max_tries):
            logger.error(
                f"Failed to connect to CAZy homepage after 10 attempts,\n"
                f"On attempt# {(tries+1)}/{max_tries} the following error was raised:\n"
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
    class_page, error = get_page(class_url)

    if class_page is None:
        logger.error(
            (
                f"Failed to connect to {class_url} after 10 attempts.\n"
                "The following error was raised:\n"
                f"{error}"
                "Could not retrieve URLs to CAZy famileis for this class.\n"
                "This class will be skipped during the scraping process."
            )
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


def parse_family(family, cazy_home, taxonomy_filters, max_tries, session):
    """Returns a Family object with Protein members, scraped from CAZy.

    Returns a Family object populated with Proteins and URLs of paginiation pages for which the
    scrape failed, including the number of times a connection to CAZy has been attempted. Also
    returns a list of strings containing URLs and associated error message of URLs for which a
    a connection to CAZy could not be made.

    :param family: Family class object, representation of CAZy family
    :param cazy_home: str, URL to CAZy home page
    :param taxonomy_filters: set of genera, species and strains to restrict the scrape to
    :param max_tries: int, number of times to reattempt scrape if an error is raised
    :param session: open SQL database session

    Return Family object, Boolean whether to scrape the family again, a list of URLs which couldn't
    connect to CAZy (and no have re-try attempts left) and list of proteins that could not be added
    to the SQL database.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Starting retrieval of proteins for {family.name} from {family.url}")

    # compile URL to first family page of protein records
    first_pagination_url = family.url.replace(".html", "_all.html")

    # check url formating of first paginiation url
    try:
        re.match(
            r"http://www.cazy.org/\D{2,3}(\d+|\d+_\d+)_all.html", first_pagination_url
        ).group()
    except AttributeError:
        logger.warning(
            f"Incorrect formatting of first protein table page URL: {first_pagination_url}\n"
            "Will not try and connect to this URL."
        )
        return(
            family,
            f"Incorrect URL format for the first protein table page for {family.name}",
            None,
        )

    first_pagination_page, error_message = get_page(first_pagination_url)

    if first_pagination_page is None:
        logger.warning(
            (
                f"Could not connect to {first_pagination_url} after 10 attempts\n"
                "The following error was raised:\n"
                f"{error_message}"
            )
        )
        return(
            family,
            [
                f"{first_pagination_url}\t"
                f"{family.cazy_class}\t"
                f"Failed to connect to first page of proteins for {family.name}\t"
                f"{error_message}"
            ],
            [],
        )

    # Get the URLS to all pages of proteins and the total number of proteins in the (sub)family
    protein_page_urls, total_proteins = get_protein_page_urls(
        first_pagination_url,
        first_pagination_page,
        cazy_home,
    )

    if len(protein_page_urls) == 0:
        return(
            family,
            [
                f"{first_pagination_url}\t"
                f"{family.cazy_class}\t"
                f"Failed to retrieve URLs to protein table pages for {family.name}\t"
                f"No specific error message availble. "
                "Failed check on line 411 of scraper.crawler.__init__.py"
            ],
            [],
        )

    failed_scrapes = []  # URLs of pages for which maximum number of scrape attempts is MET
    sql_failures = []

    for protein in tqdm(
        (y for x in (
            parse_proteins(
                url,
                family.name,
                taxonomy_filters,
                session,
            ) for url in protein_page_urls
        ) for y in x),
        total=total_proteins,
        desc=f"Parsing protein pages for {family.name}",
    ):
        if protein["url"] is not None:
            # Protein was note retrieved because could not connect to CAZy
            try:
                family.failed_pages[protein["url"]] += 1
            except KeyError:
                family.failed_pages[protein["url"]] = 1  # First failed attempt to connect to page

            if family.failed_pages[protein["url"]] == max_tries:
                # maximum attempts to connect have been reached no more attempts made
                # save the URL for writing out to a file...
                failed_scrapes.append(
                    f"{protein['url']}\t"
                    f"{family.cazy_class}\t"
                    f"Failed to connect to this page of proteins for {family.name}\t"
                    f"{protein['error']}"
                )
                # ... and do no attempt to scrape again
                del family.failed_paged[protein["url"]]

        elif protein["sql"] is not None:
            # Error occured when adding Protein to SQL database
            sql_failures.append(
                f"{protein['sql']} was not added to the database\n"
                "and raised the following error when atempting to do so:\n"
                f"{protein['error']}"
            )

    # check if any pages for the family still have attempts left for trying for a successful scrape
    if len(list(family.failed_pages.keys())) == 0:
        retry_scrape = False
    else:
        retry_scrape = True

    return family, retry_scrape, failed_scrapes, sql_failures


def get_protein_page_urls(first_pagination_url, first_pagination_page, cazy_home):
    """Retrieve the URLs to all pages containing proteins for the current working family.

    Also retrieve the total number of proteins catagloued under the family.

    :param first_pagination_url: str, URL to first page contaiing proteins
    :param first_pagination_page: BS4 object, first page containing proteins
    :param cazy_home: str, URL of CAZy homepage

    Return list of URLs, and number of proteins in the family.
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

    # Retrieve the number of proteins in the family
    data = first_pagination_page.find_all(
        "div", {"class": "pos_choix"}
    )
    search_results = re.findall(r"\(\d+\)", data[0].text)
    for index in range(len(search_results)):
        search_results[index] = search_results[index][1:-1]
        search_results[index] = int(search_results[index])
    search_results.sort()
    protein_total = search_results[-1]

    return protein_page_urls, protein_total


def parse_proteins(protein_page_url, family_name, taxonomy_filters, session):
    """Returns generator of Protein objects for all protein rows on a single CAZy family page.

    Returns a dictionary containing any errors that arose. If it an attempt to connect to CAZy
    failed the URL is stored under "url". The "sql" key is used to store the name of the protein
    which raised an SQL error further down the pipeline.

    :param protein_page_url, str, URL to the CAZy family page containing protein records
    :param family_name: str, name of CAZy family
    :param taxonomy_filters: set of genera, species and strains to restrict the scrape to
    :param session: open SQL database session

    Return generator object.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Retrieving proteins from {protein_page_url}")

    protein_page, error = get_page(protein_page_url)

    if protein_page is None:
        logger.warning(
            (
                f"Could not connect to {protein_page_url} after 10 attempts\n"
                "The following error was raised:\n"
                f"{error}\n"
                f"No protein records from this page will be retried."
            )
        )
        return {"url": protein_page_url, "error": error, "sql": None}

    # retrieve protein record table
    protein_table = protein_page.find_all("table", {"class": "listing"})[0]
    protein_rows = [
        _ for _ in protein_table.descendants if (
            (_.name == "tr") and ("id" not in _.attrs) and ("class" not in _.attrs)
        )
    ]

    for row in protein_rows:
        yield row_to_protein(row, family_name, taxonomy_filters, session)


def row_to_protein(row, family_name, taxonomy_filters, session):
    """Returns a Protein object representing a single protein row from a CAZy family protein page.

    Each row, in order, contains the protein name, EC number, source organism, GenBank ID(s),
    UniProt ID(s), and PDB accession(s).

    :param row: tr element from CAZy family protein page
    :param family_name: str, name of CAZy family
    :param taxonomy_filters: set of genera, species and strains to restrict the search to
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

    links = {"GenBank": []}
    # test for len(tds[x].contents) in case there is no link,
    # the check of .name then ensures link is captured
    if len(tds[3].contents) and tds[3].contents[0].name == "a":
        links["GenBank"] = [f"{_.get_text()}" for _ in tds[3].contents if _.name == "a"]
    if len(tds[4].contents) and tds[4].contents[0].name == "a":
        links["UniProt"] = [f"{_.get_text()}" for _ in tds[4].contents if _.name == "a"]
    if len(tds[5].contents) and tds[5].contents[0].name == "a":
        links["PDB/3D"] = [f"{_.get_text()}" for _ in tds[5].contents if _.name == "a"]

    # Retrieve non-primary GenBank accession
    # these are the accessions that are not hyerlinked in CAZy
    try:
        genbank_synonyms = tds[3].find('br').next_siblings
        for i in genbank_synonyms:
            if type(i) is bs4.element.NavigableString:
                links["GenBank"].append(i)
    except AttributeError:
        pass

    # check if UniProt or PDB accessions were retrieved. If not store as empty lists
    # this avoids KeyErros when invoking add_protein_to_db
    try:
        uniprot_accessions = links["UniProt"]
    except KeyError:
        uniprot_accessions = []
    try:
        pdb_accessions = links["PDB/3D"]
    except KeyError:
        pdb_accessions = []

    if len(links["GenBank"]) == 0:
        logger.warning(
            f"Did not retrieve any GenBank accessions for {protein_name} in {family_name}.\n"
            "The primary GenBank accession determines what unique protein the current working "
            "protein is.\n"
            "Adding protein with the GenBank accession: 'NA' as a new protein to the db."
        )
        links["GenBank"] = ["NA"]

        try:
            sql_interface.add_protein_to_db(
                protein_name,
                family_name,
                source_organism,
                'NA',
                session,
                ec_numbers=ec_numbers,
                uniprot_accessions=uniprot_accessions,
                pdb_accessions=pdb_accessions,
            )

        except Exception as error_message:
            logger.warning(f"Failed to add {protein_name} to SQL database", exc_info=1)
            return {
                "url": None,
                "error": (
                    f"Failed to add to SQL database, this error was raised: {error_message},\n"
                    f"and no GenBank listed in CAZy for this protein {protein_name}"
                ),
                "sql": protein_name,
            }

        return {
            "url": None,
            "error": (
                f"No GenBank accession listed for protein {protein_name} in CAZy,\n"
                "Protein added with the GenBank accesion 'NA'"
            ),
            "sql": f"{protein_name}",
        }

    # add protein to database
    try:
        sql_interface.add_protein_to_db(
            protein_name,
            family_name,
            source_organism,
            links["GenBank"][0],
            session,
            ec_numbers,
            links["GenBank"][1:],
            uniprot_accessions,
            pdb_accessions,
        )

    except Exception as error_message:
        logger.warning(f"Failed to add {protein_name} to SQL database", exc_info=1)
        return {
            "url": None,
            "error": f"Failed to add to SQL database. {error_message}",
            "sql": protein_name,
        }

    return {"url": None, "error": None, "sql": None}


def browser_decorator(func):
    """Decorator to retry the wrapped function up to 'retries' times."""

    def wrapper(*args, retries=10, **kwargs):
        logger = logging.getLogger(__name__)
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
                if (tries < retries):
                    logger.warning(
                        f"Failed to connect to CAZy on try {tries}/{retries}.\n"
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
