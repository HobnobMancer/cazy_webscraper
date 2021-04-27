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
"""Module for crawling through the CAZy website and parsing HTMLs for specific Taxonomy Kingdoms."""

import sys

import logging
import re

from tqdm import tqdm

from scraper.crawler import get_page, row_to_protein, row_to_protein_in_dict
from scraper.sql import sql_interface


def parse_family_by_kingdom(
    family,
    cazy_home,
    taxonomy_filters,
    kingdoms,
    ec_filters,
    args,
    session,
):
    """Parse family and add protein members to the local CAZy database, for specific tax kingdoms.

    :param family: Family class object, representation of CAZy family
    :param cazy_home: str, URL to CAZy home page
    :param taxonomy_filters: set of genera, species and strains to restrict the scrape to
    :param kingdoms: list of taxonomy Kingdoms to restrict the scrape to
        OR a str if the user did not define any Kingdoms to scrape  (enables faster
        scraping via the 'all' pages)
    :param ec_filters: set of EC filters to limit scrape to
    :param args: cmd-line args parser
    :param session: open SQL database session

    Return:
        Family object,
        Boolean whether to scrape the family again,
        List of URLs which couldn't connect to CAZy (and no have re-try attempts left),
        List of proteins that could not be added to the SQL database.
        List of pages whose parsing did match the expected format
    """
    # define lists for storing error messages during scraping and parsing proteins
    failed_scrapes = []  # URLs of pages for which maximum number of scrape attempts is MET
    sql_failures = []  # Errors and protein names for proteins that raised SQL errors
    format_failures = []  # pages where no proteins were retrieved, e.g. URL formated incorrectly

    # scraping for the first time
    if len(list(family.failed_pages.keys())) == 0:
        for kingdom in kingdoms:
            # compile first url to the first pagination page
            first_pagination_url = family.url.replace(".html", f"_{kingdom}.html")

            # retrieve URLs to all pagination pages and retrieve CAZymes from the pages
            (
                family, failed_scrapes, sql_failures, format_failures,
            ) = parse_kingdom_first_pagination_page(
                first_pagination_url,
                family,
                kingdom,
                failed_scrapes,
                sql_failures,
                format_failures,
                cazy_home,
                taxonomy_filters,
                ec_filters,
                args,
                session
            )

    else:
        # Rescrape pages that failed previously
        # family.failed_pages = {kingdon: {url: #_of_tries}}
        for kingdom in family.failed_pages:
            protein_page_urls = family.failed_pages[kingdom]

            page_urls = []  # store pages to retrieve CAZymes from

            for url in protein_page_urls:
                # check if need to retry parsing the first pagination page
                try:  # {kingdon: {url: ["first_pagination", #_of_tries]}}
                    if protein_page_urls[url][0] == "first_pagination":
                        # retry parsing the first pagination page to get URLs to all pag pages
                        first_pagination_url = url
                        (
                            family, failed_scrapes, sql_failures, format_failures,
                        ) = parse_kingdom_first_pagination_page(
                            first_pagination_url,
                            family,
                            kingdom,
                            failed_scrapes,
                            sql_failures,
                            format_failures,
                            cazy_home,
                            taxonomy_filters,
                            ec_filters,
                            args,
                            session
                        )

                except (TypeError, IndexError):  # do not need to reparse the first pagination page
                    page_urls.append(url)

            if len(page_urls) != 0:
                # retrieve CAZymes from pages that could not be scrapped previously
                total_proteins = len(page_urls * 100)  # provides approximation of protein total

                if type(session) is dict:
                    family, failed_scrapes = parse_protein_pages_dict(
                        family,
                        kingdom,
                        page_urls,
                        taxonomy_filters,
                        ec_filters,
                        failed_scrapes,
                        session,
                        args,
                    )
                else:
                    family, failed_scrapes, sql_failures = parse_protein_pages(
                        family,
                        kingdom,
                        protein_page_urls,
                        total_proteins,
                        taxonomy_filters,
                        ec_filters,
                        failed_scrapes,
                        sql_failures,
                        session,
                        args,
                    )


    # check if any pages have attempts left for retrying for a successful scrape
    if len(list(family.failed_pages.keys())) == 0:
        retry_scrape = False
    else:
        retry_scrape = True

    return family, retry_scrape, failed_scrapes, sql_failures, format_failures, session


def parse_kingdom_first_pagination_page(
    first_pagination_url,
    family,
    kingdom,
    failed_scrapes,
    sql_failures,
    format_failures,
    cazy_home,
    taxonomy_filters,
    ec_filters,
    args,
    session,
):
    """Parse the first pagination page for the current working Kingomd for the CAZy family.

    Coordinate getting all pagination page URLs, and then parse all pagination pages.

    :param first_pagination_url: URL to the first pagination page
    :param family: Family class instance
    :param kingdom: str, taxonomic kingdom currently being scraped
    :param failed_scrapes: list of URLs that could not be scraped
    :param sql_failures: list of proteins that raised errors when being added to the db
    :param format_failures: list of URLs with incorrect formating for scraping
    :param cazy_home: str, CAZy home page URL
    :param taxonomy_filters: taxonomic filters to restrict the retrieval of CAZymes to
    :param ec_filters: list of EC numbers to restrict the scrape to
    :param args: cmd-line args parser
    :param session: open SQL db session

    Return Family class instance, Bool to rescrape or not, list of failed scrapes, list of SQL
    errors raised during adding CAZymes to the db, and URLs with incorrect formats.
    """
    logger = logging.getLogger(__name__)
    # check url formating of first pagination url
    try:
        re.match(
            rf"http://www.cazy.org/(\D\D|\D\D\D)(\d+_{kingdom}|\d+_\d+_{kingdom}).html",
            first_pagination_url,
        ).group()
    except AttributeError:
        logger.warning(
            f"Incorrect formatting of first protein table page URL: {first_pagination_url}\n"
            "Will not try and connect to this URL."
        )
        format_failures.append(
            f"{first_pagination_url}\tIncorrect URL format for the first protein table page "
            f"for {family.name} in kingdom {kingdom}, could not scrape any CAZymes for this "
            f"Kingdom for family {family.name}"
        )
        return family, failed_scrapes, sql_failures, format_failures

    # Retrieve the URLs to all pagination pages for the CAZy family
    protein_page_urls, total_proteins = get_kingdom_pagination_data(
        first_pagination_url, family, kingdom, cazy_home, args, session,
    )

    if type(total_proteins) is not int:  # errors were raised and need parsing
        # protein_page_urls is a dict containing reasons for failed retreival of total proteins
        (
            family, failed_scrapes, sql_failures, format_failures,
        ) = parse_kingdom_total_proteins_error(
            first_pagination_url,
            family,
            kingdom,
            failed_scrapes,
            sql_failures,
            format_failures,
            total_proteins,
            args,
        )
        return family, failed_scrapes, sql_failures, format_failures

    # iterate through the pages and retrieve proteins, adding to the local CAZy data
    if type(session) is dict:
        family, failed_scrapes = parse_protein_pages_dict(
            family,
            kingdom,
            page_urls,
            taxonomy_filters,
            ec_filters,
            failed_scrapes,
            session,
            args,
        )
    else:
        family, failed_scrapes, sql_failures = parse_protein_pages(
            family,
            kingdom,
            protein_page_urls,
            total_proteins,
            taxonomy_filters,
            ec_filters,
            failed_scrapes,
            sql_failures,
            session,
            args,
        )
    return family, failed_scrapes, sql_failures, format_failures


def get_kingdom_pagination_data(first_pagination_url, family, kingdom, cazy_home, args, session):
    """Parse the first pagination page and retrieve URLs to all pagination page for the Family.

    :param first_pagination_url: str, URL to the fist page of the family
    :param family: Family class instance, represents a unique CAZy family
    :param kingdom: str, taxonomic kingdom currently being scraped
    :param cazy_home: str, URL to CAZy home page
    :param args: cmd-line args parser
    :param session: open SQL database session

    Return dict of error messages if errors arise OR list of URLs to pagination pages, and the
    total number of proteins in the family.
    """
    logger = logging.getLogger(__name__)

    # retrieve a list of all page urls of protein tables for the CAZy family
    first_pagination_page, error_message = get_page(
        first_pagination_url, args, max_tries=args.retries,
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
                    f"{first_pagination_url}\t{family.cazy_class}\t{kingdom}\t"
                    f"Failed to connect to first pagination page for {family.name}, therefore "
                    f"could not retrieve URLs to all pagination pages\t{error_message}"
                ),
                "format": None,
            },
            None,
        )

    # Get the URLS to all pages of proteins and the total number of proteins in the (sub)family
    protein_page_urls, total_proteins = get_kingdom_page_urls(
        first_pagination_url, first_pagination_page, kingdom, cazy_home, family.name,
    )

    if total_proteins == 'Deleted family!':
        # add family to the database
        logger.warning(
            f'{family.name} listed as "Deleted family" in CAZy.\n'
            'Adding family name to the database'
        )
        sql_interface.add_deleted_cazy_family(family.name, session)
        return(
            {
                "url": None,
                "format": (
                    f"{first_pagination_url}\t{family.cazy_class}\t"
                    f"{family.name} listed as 'Deleted family' in CAZy.\t"
                    "Added family name to the database"
                ),
            },
            None
        )

    elif total_proteins == 'Empty family':
        # add family to the database
        logger.warning(
            f'{family.name} is an empty family, but not listed "Deleted family" in CAZy.\n'
            'Adding family name to the database'
        )
        sql_interface.add_deleted_cazy_family(family.name, session)
        return(
            {
                "url": None,
                "format": (
                    f"{first_pagination_url}\t{family.cazy_class}\t"
                    f"{family.name} is an empty family in CAZy, but not listed as Deleted fam.\t"
                    "Added family name to the database"
                ),
            },
            None
        )

    elif total_proteins == 'Failed Retrieval':
        # add family to the database
        logger.warning(
            f"Could not retrieve total protiens count for {family.name}, appears to be empty\n"
            "but not listed 'Deleted family' in CAZy.\nAdding family name to the database"
        )
        sql_interface.add_deleted_cazy_family(family.name, session)
        return(
            {
                "url": None,
                "format": (
                    f"{first_pagination_url}\t{family.cazy_class}\t{family.name} Could not "
                    "retrieve total protiens count\tAdded family name to the database"
                ),
            },
            None
        )

    elif len(protein_page_urls) == 0:
        logger.warning(f"No protein page URLs found for {family.name}: {first_pagination_url}")
        return(
            {
                "url": None,
                "format": (
                    f"{first_pagination_url}\t{family.cazy_class}\t"
                    f"Failed to retrieve URLs to protein table pages for {family.name}\t"
                    f"No specific error message availble."
                ),
            },
            None
        )

    else:
        return protein_page_urls, total_proteins


def get_kingdom_page_urls(
    first_pagination_url,
    first_pagination_page,
    kingdom,
    cazy_home,
    family_name,
):
    """Retrieve the URLs to all pages containing proteins for the current working family.

    Also retrieve the total number of proteins catagloued under the family.

    :param first_pagination_url: str, URL to first page contaiing proteins
    :param first_pagination_page: BS4 object, first page containing proteins
    :param kingdom: str, taxonomy Kingdom
    :param cazy_home: str, URL of CAZy homepage
    :param family_name: str, name of the CAZy family

    Return list of URLs, and number of proteins in the family.
    """
    logger = logging.getLogger(__name__)

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

    # Retrieve the number of proteins in the family: kingdom from the hyperlinks to the start of
    # each Kingdom's pagination pages
    data = first_pagination_page.find_all("div", {"class": "pos_choix"})
    # retrieve the text for the Kingdom's hyperlink, and retrieve the number within the brackets
    try:
        protein_total = int(re.findall(
            rf"{kingdom} \(\d+\)", data[0].text, flags=re.IGNORECASE,
        )[0].split("(")[1][:-1])
    except IndexError:
        try:
            family_activities_cell = first_pagination_page.select("table")[
                0].select("tr")[0].select('td')[0].contents[0].strip()

            if family_activities_cell == 'Deleted family!':
                protein_total = 'Deleted family!'
            else:
                logger.warning(f"No proteins found for 'all' in {family_name}")
                protein_total = 'Empty family'
        except Exception:
            logger.warning(f"No proteins found for 'all' in {family_name}")
            protein_total = 'Failed Retrieval'

    return protein_page_urls, protein_total


def parse_kingdom_total_proteins_error(
    first_pagination_url,
    family,
    kingdom,
    failed_scrapes,
    sql_failures,
    format_failures,
    errors,
    args,
):
    """Parse when total proteins for CAZy family was not retrieved.

    :param first_pagination_url: str, URL to first pagination page of CAZy family
    :param family: CAZy Family class instance
    :param kingdom: str, current working taxonomic kingdom
    :param failed_scrapes: list of URLs that could not be scraped
    :param sql_failures: list of proteins that raised errors when being added to the db
    :param format_failures: list of URLs with incorrect formating for scraping
    :param total_proteins_errors: dict of errors raised when retrieving total protein count.
    :param args: cmd-line args parser

    Return Family class instance, list of failed scrapes, list of SQL
    errors raised during adding CAZymes to the db, and URLs with incorrect formats.
    """
    if errors["url"] is not None:  # Could not connect to CAZy, and should try connection later
        try:
            family.failed_pages[kingdom]

            try:
                family.failed_pages[kingdom][first_pagination_url][1] += 1

            except KeyError:
                family.failed_pages[kingdom] = {first_pagination_url: ["first_pagination", 1]}

        except KeyError:
            family.failed_pages = {kingdom: {first_pagination_url: ["first_pagination", 1]}}

        if family.failed_pages[kingdom][first_pagination_url][1] >= args.retries:
            failed_scrapes.append(
                f"{first_pagination_url}\t{family.cazy_class}\t{kingdom}\t{family.name}\t"
                f"Failed to connect to this page of proteins for {family.name}, and "
                f"raised the following error message:\n{errors['url']}"
            )
            # do not attempt to scrape again becuase ran out of attempts
            del family.failed_pages[kingdom][first_pagination_url]

    if errors["format"] is not None:  # URL incorrectly formatted, so cannot scrape
        format_failures.append(f"{first_pagination_url}\t{errors['format']}")

    return family, failed_scrapes, sql_failures, format_failures


def parse_protein_pages(
    family,
    kingdom,
    protein_page_urls,
    total_proteins,
    taxonomy_filters,
    ec_filters,
    failed_scrapes,
    sql_failures,
    session,
    args,
):
    """Parse HTML tables of proteins and add to the local CAZyme database.

    :param family: Family class instance
    :param kingdom: str, current working taxonomic kingdom
    :param protein_page_urls: list of URLs of pages to parse
    :param total_proteins: int, number of proteins to parse
    :param taxonomy_filters: set, genera, species and strains to limit the scrape to
    :param ec_filters: list of EC numbers to limit the retrieval of CAZymes to
    :param failed_scrapes: list of scrapes that failed
    :param sql_failures: list of proteins that raised an error when being added to the sQL database.
    :param session: open SQL database session
    :param args: cmdline args parser

    Return family, failed_scrapes, sql_failures"""
    # Retrieve the CAZymes (proteins) and write to the local database
    for protein in tqdm(
        (y for x in (
            parse_kingdom_protein_tables(
                url,
                family.name,
                kingdom,
                taxonomy_filters,
                ec_filters,
                args,
                session,
            ) for url in protein_page_urls
        ) for y in x),
        total=total_proteins,
        desc=f"Parsing protein pages for {family.name}: {kingdom}",
    ):
        if protein["url"] is not None:
            # Protein was not retrieved because could not connect to CAZy
            try:
                family.failed_pages[kingdom]

                try:
                    family.failed_pages[kingdom][protein["url"]] += 1
                except KeyError:  # first failed scrape for the specific pagination page
                    family.failed_pages[kingdom][protein["url"]] = 1

            except KeyError:  # first failed attempt for the family:kingdom
                family.failed_pages[kingdom] = {protein["url"]: 1}

            if family.failed_pages[kingdom][protein["url"]] >= (args.retries + 1):
                # Reached maximum attempts number of attempted connections ...
                failed_scrapes.append(
                    f"{protein['url']}\t{family.cazy_class}\t"
                    f"Failed to connect to this page of proteins for {family.name}\t"
                    f"{protein['error']}"
                )
                # ... and do no attempt to scrape again
                del family.failed_pages[kingdom][protein["url"]]

        if protein["sql"] is not None:  # Error occured when adding Protein to SQL database
            sql_failures.append(
                f"{protein['sql']} was not added to the database, and raised the following "
                f"error when atempting to do so:\n{protein['error']}"
            )
    return family, failed_scrapes, sql_failures


def parse_kingdom_protein_tables(
    protein_page_url,
    family_name,
    kingdom,
    taxonomy_filters,
    ec_filters,
    args,
    session,
):
    """Returns generator of Protein objects for all protein rows on a single CAZy family page.

    Returns a dictionary containing any errors that arose. If it an attempt to connect to CAZy
    failed the URL is stored under "url". The "sql" key is used to store the name of the protein
    which raised an SQL error further down the pipeline.

    :param protein_page_url, str, URL to the CAZy family page containing protein records
    :param family_name: str, name of CAZy family
    :param kingdom: str, Taxonomy Kingdom of CAZymes currently being scraped
    :param taxonomy_filters: set of genera, species and strains to restrict the scrape to
    :param ec_filters: list of EC numbers to restrict the scraping of CAZymes to
    :param args: cmd-line args parser
    :param session: open SQL database session

    Return generator object.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Retrieving proteins from {protein_page_url}")

    protein_page, error = get_page(protein_page_url, args, max_tries=args.retries)

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

    # Get the table on the page corresponding to CAZymes
    try:
        cazyme_table = protein_page.select("table")[1]
    except AttributeError:
        logger.warning("NO PROTEIN TABLE!")
        # raised if there is not table of proteins
        return {"url": 'No protein table', "error": None, "sql": None}

    # check page was returned correctly. Sometimes CAZy returns the page but the table of CAZymes
    # is not populated with proteins. This often occurs after multiple timedout conenctions to CAZy
    if len(cazyme_table.select("tr")) == 1:
        # double check it is not a deleted family
        fam_activities = protein_page.select("table")[0].select("tr")[0].select("td")[0].contents
        if type(fam_activities) is list:
            if fam_activities[0] == "Deleted family!":
                return {"url": None, "sql": None, "format": "Deleted family in CAZy"}
        else:
            if fam_activities == "Deleted family!":
                return {"url": None, "sql": None, "format": "Deleted family in CAZy"}

        # check protein total to be sure there should be proteins
        data = protein_page.find_all("div", {"class": "pos_choix"})
        # retrieve the text for the Kingdom's hyperlink, and retrieve the number within the brackets
        try:
            protein_total = int(re.findall(
                rf"{kingdom} \(\d+\)", data[0].text, flags=re.IGNORECASE,
            )[0].split("(")[1][:-1])
        except IndexError:
            logger.warning(
                f"No proteins in {kingdom} in {family_name}"
            )
            protein_total = 0

        if protein_total == 0:
            logger.warning(f"{protein_page_url} is an empty in CAZy")
            return {"url": None, "sql": None, "format": "Empty family in CAZy, family added to db"}

        logger.warning(
            f"Previously connected to {protein_page_url} but complete page was not returned.\n"
            "Retrying connection to CAZy to retrieve the complete page."
        )

        success = False

        while success is False:
            protein_page, error = get_page(protein_page_url, args, max_tries=args.retries)

            if protein_page is not None:
                cazyme_table = protein_page.select("table")[1]
                if len(cazyme_table.select("tr")) != 1:
                    success = True
                else:
                    logger.warning(
                        f"Previously connected to {protein_page_url} but complete page was not"
                        "returned.\nRetrying connection to CAZy to retrieve the complete page."
                    )

            else:
                logger.warning(
                    (
                        f"Could not connect to {protein_page_url} after {args.retries} attempts\n"
                        "The following error was raised:\n"
                        f"{error}\n"
                        f"Adding page to {family_name} to be rescraped if try attempts are left."
                    )
                )
                {"url": [protein_page_url, error], "sql": None, "format": None}

    # Get all rows in the table, exluding those that are header/navigation roes
    # Each row has an .attrs attribute, and this holds the selector options, such as
    # "class" and "id"; CAZyme rows don't have these attributes
    cazyme_rows = [
        _ for _ in cazyme_table.select("tr") if "class" not in _.attrs and "id" not in _.attrs
    ]

    # Loop overal all rows and add data to the local database
    for row in cazyme_rows:
        yield row_to_protein(row, family_name, taxonomy_filters, kingdom, ec_filters, session, args)


def parse_protein_pages_dict(
    family,
    kingdom,
    protein_page_urls,
    taxonomy_filters,
    ec_filters,
    failed_scrapes,
    session,
    args,
):
    """Parse HTML pages of proteins and add to the local CAZyme database.

    :param family: Family class instance
    :param kingdom: str, current working taxonomic kingdom
    :param protein_page_urls: list of URLs of pages to parse
    :param taxonomy_filters: set, genera, species and strains to limit the scrape to
    :param ec_filters: list of EC numbers to limit the retrieval of CAZymes to
    :param failed_scrapes: list of scrapes that failed
    :param session: open SQL database session
    :param args: cmdline args parser

    Return family, failed_scrapes, sql_failures"""
    logger = logging.getLogger(__name__)

    for url in tqdm(protein_page_urls, desc=f"Parsing {family.name} {kingdom} pages"):
        protein_page, error = get_page(url, args, max_tries=args.retries)

        if protein_page is None:
            warning = (
                f"Could not connect to {url} after 10 attempts\n"
                "The following error was raised:\n"
                f"{error}\n"
                f"No protein records from this page will be retried."
            )
            logger.warning(warning)
            family, failed_scrapes = parse_url_error_kngdm(
                url, family, kingdom, failed_scrapes, warning, args,
            )
            continue

        # Get the table on the page corresponding to CAZymes
        try:
            cazyme_table = protein_page.select("table")[1]
        except (AttributeError, IndexError):
            warning = f"No protein table found on page {url}"
            logger.warning(warning)
            family, failed_scrapes = parse_url_error_kngdm(
                url, family, kingdom, failed_scrapes, warning, args,
            )
            continue

        # check page was returned correctly. Sometimes CAZy returns the page but the table of 
        # CAZymes is not populated. This often occurs after multiple timedout conenctions to CAZy
        if len(cazyme_table.select("tr")) == 1:
            # double check it is not a deleted family
            fam_activities = protein_page.select("table")[0].select("tr")[0].select("td")[0].contents
            if type(fam_activities) is list:
                if fam_activities[0] == "Deleted family!":
                    failed_scrapes.append(f"{family.name} is a Deleted family in CAZy")
                    continue
            else:
                if fam_activities == "Deleted family!":
                    failed_scrapes.append(f"{family.name} is a Deleted family in CAZy")
                    continue

            # check protein total to be sure there should be proteins
            data = protein_page.find_all("div", {"class": "pos_choix"})
            # retrieve the text for the Kingdom's hyperlink, and retrieve the number within the brackets
            try:
                protein_total = int(re.findall(
                    rf"{kingdom} \(\d+\)", data[0].text, flags=re.IGNORECASE,
                )[0].split("(")[1][:-1])
            except IndexError:
                logger.warning(
                    f"No proteins in {kingdom} in {family.name}"
                )
                protein_total = 0

            if protein_total == 0:
                logger.warning(f"{url} is an empty in CAZy")
                failed_scrapes.append(f"{family.name} is an empty family in CAZy")
                continue

            else:
                logger.warning(
                    f"Previously connected to {url} but complete page was not returned.\n"
                    "Retrying connection to CAZy to retrieve the complete page."
                )

                protein_page, error = get_page(url, args, max_tries=args.retries)

                if protein_page is None:
                    warning = (
                        f"Could not connect to {url} after 10 attempts\n"
                        "The following error was raised:\n"
                        f"{error}\n"
                        f"No protein records from this page will be retried."
                    )
                    logger.warning(warning)
                    family, failed_scrapes = parse_url_error_kngdm(
                        url, family, kingdom, failed_scrapes, warning, args,
                    )
                    continue

                cazyme_table = protein_page.select("table")[1]
                if len(cazyme_table.select("tr")) == 1:
                    warning = (
                        f"Previously connected to {url} but complete page was not returned.\nAdded "
                        f"{family.name} to be rescraped later"
                    )
                    logger.warning(warning)
                    family, failed_scrapes = parse_url_error_kngdm(
                        url, family, kingdom, failed_scrapes, warning, args,
                    )
                    continue

        # Get all rows in the table, exluding those that are header/navigation roes
        # Each row has an .attrs attribute, and this holds the selector options, such as
        # "class" and "id"; CAZyme rows don't have these attributes
        cazyme_rows = [
            _ for _ in cazyme_table.select("tr") if "class" not in _.attrs and "id" not in _.attrs
        ]

        # Loop overal all rows and add data to the local database
        for row in cazyme_rows:
            new_report_dict, session = row_to_protein_in_dict(
                row, family.name, taxonomy_filters, ec_filters, session
            )

            if new_report_dict["error"] is not None:
                failed_scrapes.append(
                    f"Protein {new_report_dict['protein']} raised the following error when being "
                    f"scrapped:\n{new_report_dict['error']}"
                )

    return family, failed_scrapes, session


def parse_url_error_kngdm(url, family, kingdom, failed_scrapes, error_message, args):
    """Parse connection errors when trying to parse protein table pages.

    :param url: str, url of HTML webpage containing the current working protein table page
    :param family: Family class instance, represents a CAZy family
    :param kingdom: str, taxonomic kingdom of proteins being scraped
    :param failed_scrapes: list of erros raised when trying to scrape CAZy
    :param error_message: str, error raised when trying to scrape CAZy

    Return failed_scrapes and Family class instance"""
    try:
        family.failed_pages[kingdom]

        try:
            family.failed_pages[kingdom][url] += 1
        except KeyError:  # first failed scrape for the specific pagination page
            family.failed_pages[kingdom][url] = 1

    except KeyError:  # first failed attempt for the family:kingdom
        family.failed_pages[kingdom] = {url: 1}

    if family.failed_pages[kingdom][url] >= (args.retries + 1):
        # Reached maximum attempts number of attempted connections ...
        failed_scrapes.append(
            f"{url}\t{family.cazy_class}\t"
            f"Failed to connect to this page of proteins for {family.name}\t{error_message}"
        )
        # ... and do no attempt to scrape again
        del family.failed_pages[kingdom][url]

    return family, failed_scrapes
