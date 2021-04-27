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
"""Module for crawling through the CAZy website and parsing HTMLs under the 'ALL' tab."""


import logging
import re

from tqdm import tqdm

from scraper.crawler import get_page, row_to_protein, row_to_protein_in_dict
from scraper.sql import sql_interface


def parse_family_via_all_pages(family, cazy_home, taxonomy_filters, ec_filters, args, session):
    """Parse the protein tables from the 'all' pages of CAZy family.

    CAZy families have separate series of HTML tables for each taxonomy Kingdom, and a single series
    of HTML tables containing all proteins from all Kingdoms, called 'all'. These 'all' pages
    containing 1000 proteins per page and thus scraping these pages is significantly faster than
    scraping all the specific Kingdom HTML tables (which hold only 100 proteins each, and thus
    require more calls to CAZy).

    :param family: Family class object, representation of CAZy family
    :param cazy_home: str, URL to CAZy home page
    :param taxonomy_filters: set of genera, species and strains to restrict the scrape to
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

    # check if scraping family for the first time, and there aren't pages to try scraping
    if len(list(family.failed_pages.keys())) == 0:
        first_pagination_url = family.url.replace(".html", "_all.html")
        (
            family, failed_scrapes, sql_failures, format_failures, session,
        ) = parse_first_paginiation_page(
            first_pagination_url,
            family,
            failed_scrapes,
            sql_failures,
            format_failures,
            cazy_home,
            args,
            taxonomy_filters,
            ec_filters,
            session,
        )

    # check if need to retry parsing the first paginiation page to get the URLs the pagination pages
    if len(list(family.failed_pages.keys())) == 1:
        protein_page_urls = list(family.failed_pages.keys())

        try:
            first_url = protein_page_urls[0]
            family.failed_pages[first_url][0]
            if family.failed_pages[first_url][0] == "first_pagination":
                # retry connecting and parsing the first paginiation page, jumps to line 125
                first_pagination_url = protein_page_urls[0]
                (
                    family, failed_scrapes, sql_failures, format_failures, session
                ) = parse_first_paginiation_page(
                    first_pagination_url,
                    family,
                    failed_scrapes,
                    sql_failures,
                    format_failures,
                    cazy_home,
                    args,
                    taxonomy_filters,
                    ec_filters,
                    session,
                )

        except (IndexError, TypeError):
            # else not the first paginiation page, so retry getting CAZyme data
            total_proteins = 1000
            if type(session) is dict:
                family, failed_scrapes, session = parse_protein_table_dict(
                    family,
                    kingdom,
                    protein_page_urls,
                    taxonomy_filters,
                    ec_filters,
                    failed_scrapes,
                    session,
                    args,
                )

            else:
                family, failed_scrapes, sql_failures, format_failures = parse_proteins(
                    protein_page_urls,
                    total_proteins,
                    family,
                    taxonomy_filters,
                    ec_filters,
                    session,
                    args,
                    failed_scrapes,
                    sql_failures,
                    format_failures,
                )

    else:  # rescraping selected pages of the family
        protein_page_urls = list(family.failed_pages.keys())
        total_proteins = len(protein_page_urls * 1000)  # provides approximation of protein total
        if type(session) is dict:
            family, failed_scrapes, session = parse_protein_table_dict(
                family,
                kingdom,
                protein_page_urls,
                taxonomy_filters,
                ec_filters,
                failed_scrapes,
                session,
                args,
            )

        else:
            family, failed_scrapes, sql_failures, format_failures = parse_proteins(
                protein_page_urls,
                total_proteins,
                family,
                taxonomy_filters,
                ec_filters,
                session,
                args,
                failed_scrapes,
                sql_failures,
                format_failures,
            )

    # check if any pages have attempts left for retrying for a successful scrape
    if len(list(family.failed_pages.keys())) == 0:
        retry_scrape = False
    else:
        retry_scrape = True

    return family, retry_scrape, failed_scrapes, sql_failures, format_failures, session


def parse_first_paginiation_page(
    first_pagination_url,
    family,
    failed_scrapes,
    sql_failures,
    format_failures,
    cazy_home,
    args,
    taxonomy_filters,
    ec_filters,
    session,
):
    """Parse the first paginiation page.

    Coordinate getting all paginiation page URLs, and then parse all paginiation pages.

    :param first_pagination_url: URL to the first paginiation page
    :param family: Family class instance
    :param failed_scrapes: list of URLs that could not be scraped
    :param sql_failures: list of proteins that raised errors when being added to the db
    :param format_failures: list of URLs with incorrect formating for scraping
    :param cazy_home: str, CAZy home page URL
    :param args: cmd-line args parser
    :param taxonomy_filters: taxonomic filters to restrict the retrieval of CAZymes to
    :param ec_filters: list of EC numbers to restrict the scrape to
    :param session: open SQL db session

    Return Family class instance, list of failed scrapes, list of SQL
    errors raised during adding CAZymes to the db, and URLs with incorrect formats.
    """
    logger = logging.getLogger(__name__)
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
        format_failures.append(
                f"{first_pagination_url}\tIncorrect URL format therefore could not retrieve "
                f"proteins from CAZy family {family.name}\n"
        )
        return family, failed_scrapes, sql_failures, format_failures, session

    # retrieve pages to other pagination pages of protein tables for the family
    protein_page_urls, total_proteins = get_paginiation_data(
        first_pagination_url, family, cazy_home, args, session,
    )

    if type(total_proteins) is not int:  # errors were raised and need parsing
        # protein_page_urls is a dict containing reasons for failed retreival of total proteins
        (
            family, failed_scrapes, sql_failures, format_failures
        ) = parse_total_proteins_error(
            first_pagination_url,
            family,
            failed_scrapes,
            sql_failures,
            format_failures,
            protein_page_urls,
            args,
        )
        return family, failed_scrapes, sql_failures, format_failures, session

    # iterate through the pages and retrieve proteins, adding to the local CAZy data
    if type(session) is dict:
        family, failed_scrapes, session = parse_protein_table_dict(
            family,
            kingdom,
            protein_page_urls,
            taxonomy_filters,
            ec_filters,
            failed_scrapes,
            session,
            args,
        )

    else:
        family, failed_scrapes, sql_failures, format_failures = parse_proteins(
            protein_page_urls,
            total_proteins,
            family,
            taxonomy_filters,
            ec_filters,
            session,
            args,
            failed_scrapes,
            sql_failures,
            format_failures,
        )
    return family, failed_scrapes, sql_failures, format_failures, session


def get_paginiation_data(first_pagination_url, family, cazy_home, args, session):
    """Parse the first paginiation page and retrieve URLs to all pagination page for the Family.

    :param first_pagination_url: str, URL to the fist page of the family
    :param family: Family class instance, represents a unique CAZy family
    :param cazy_home: str, URL to CAZy home page
    :param args: cmd-line args parser
    :param session: open SQL database session

    Return dict of error messages if errors arise OR list of URLs to paginiation pages, and the
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
                    f"{first_pagination_url}\t{family.cazy_class}\t"
                    f"Failed to connect to first pagination page for {family.name}, therefore "
                    f"could not retrieve URLs to all paginiation pages\t{error_message}"
                ),
                "format": None,
            },
            None,
        )

    # Get the URLS to all pages of proteins and the total number of proteins in the (sub)family
    protein_page_urls, total_proteins = get_paginiation_page_urls(
        first_pagination_url, first_pagination_page, cazy_home, family.name,
    )

    if total_proteins == 'Deleted family!':
        # add family to the database
        logger.warning(
            f'{family.name} listed as "Deleted family" in CAZy.\nAdding family name to the database'
        )
        sql_interface.add_deleted_cazy_family(family.name, session)
        return(
            {
                "url": None,
                "format": (
                    f"{first_pagination_url}\t{family.cazy_class}\t{family.name} listed as "
                    "'Deleted family' in CAZy.\tAdded family name to the database"
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
                    f"{first_pagination_url}\t{family.cazy_class}\t{family.name} is empty in CAZy, "
                    "but not listed as Deleted fam.\tAdded family name to the database"
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
                    f"{first_pagination_url}\t{family.cazy_class}\tFailed to retrieve URLs to "
                    f"protein table pages for {family.name}\tNo specific error message availble."
                ),
            },
            None
        )

    else:
        return protein_page_urls, total_proteins


def get_paginiation_page_urls(first_pagination_url, first_pagination_page, cazy_home, family_name):
    """Retrieve the URLs to all pages containing proteins for the current working family.

    Also retrieve the total number of proteins catagloued under the family.

    :param first_pagination_url: str, URL to first page contaiing proteins
    :param first_pagination_page: BS4 object, first page containing proteins
    :param cazy_home: str, URL of CAZy homepage
    :param family_name: str, name of the current working CAZy family

    Return list of URLs, and the number of proteins in the family.
    """
    logger = logging.getLogger(__name__)

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

    # retrieve the data element that contains the links the sets of HTML tables
    data = first_pagination_page.find_all("div", {"class": "pos_choix"})

    # retrieve the number of proteins listed after 'all' in the data
    try:
        protein_total = int(re.findall(
            r"all \(\d+\)", data[0].text, flags=re.IGNORECASE,
        )[0].split("(")[1][:-1])
    except IndexError:
        # check if the family has been deleted
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


def parse_total_proteins_error(
    first_pagination_url,
    family,
    failed_scrapes,
    sql_failures,
    format_failures,
    errors,
    args,
):
    """Parse when total proteins for CAZy family was not retrieved.

    :param first_pagination_url: str, URL to first pagination page of CAZy family
    :param family: CAZy Family class instance
    :param failed_scrapes: list of URLs that could not be scraped
    :param sql_failures: list of proteins that raised errors when being added to the db
    :param format_failures: list of URLs with incorrect formating for scraping
    :param total_proteins_errors: dict of errors raised when retrieving total protein count.
    :param args: cmd-line args parser

    Return Family class instance, list of failed scrapes, list of SQL
    errors raised during adding CAZymes to the db, and URLs with incorrect formats.
    """
    if errors["url"] is not None:  # means could not connect to CAZy
        try:
            family.failed_pages[first_pagination_url][1] += 1
        except KeyError:
            family.failed_pages[first_pagination_url] = ["first_pagination", 1]  # first attempt

        if family.failed_pages[first_pagination_url][1] >= args.retries:
            failed_scrapes.append(
                f"{first_pagination_url}\t{family.cazy_class}\t"
                f"Failed to connect to this page of proteins for {family.name}, and "
                f"raised the following error message:\n{errors['url']}"
            )
            # do not attempt to scrape again becuase ran out of attempts
            del family.failed_pages[first_pagination_url]

    if errors["format"] is not None:  # URL incorrectly formatted, so cannot scrape
        format_failures.append(f"{first_pagination_url}\t{errors['format']}")

    return family, failed_scrapes, sql_failures, format_failures


def parse_proteins(
    protein_page_urls,
    total_proteins,
    family,
    taxonomy_filters,
    ec_filters,
    session,
    args,
    failed_scrapes,
    sql_failures,
    format_failures,
):
    """Parse proteins in HTML tables from CAZy.

    :param protein_page_urls: list of URLs to pages containing tables of proteins
    :param total_proteins: int, number of proteins in the CAZy family
    :param family: Family class object, represents a CAZy family
    :param taxonomy_filters: taxonomic filters to restrict the retrieval of CAZymes to
    :param ec_filters: list of EC numbers to restrict the scrape to
    :param session: open SQL db session
    :param args: cmd-line args parser
    :param failed_scrapes: list of URLs that could not be scraped
    :param sql_failures: list of proteins that raised errors when being added to the db
    :param format_failures: list of URLs with incorrect formating for scraping

    Return Family class instance, Bool to rescrape or not, list of failed scrapes, list of SQL
    errors raised during adding CAZymes to the db, and URLs with incorrect formats.
    """
    for protein in tqdm(
        (y for x in (
            parse_protein_table(
                url,
                family.name,
                taxonomy_filters,
                ec_filters,
                session,
                args,
            ) for url in protein_page_urls
        ) for y in x),
        total=total_proteins,
        desc=f"Parsing protein pages for {family.name}",
    ):
        if protein["url"] is not None:  # Protein not retrieved because couldn't connect to CAZy
            # protein["url"][0] = URL
            # protein["url"][0] = Error message
            try:
                family.failed_pages[protein["url"][0]] += 1
            except KeyError:
                family.failed_pages[protein["url"][0]] = 1  # First failed attempt to connect

            if family.failed_pages[protein["url"][0]] == args.retries:
                # maximum attempts to connect have been reached no more attempts made, write to file
                failed_scrapes.append(
                    f"{protein['url'][0]}\t{family.cazy_class}\t"
                    f"Failed to connect to this page of proteins for {family.name}, "
                    f"and raised the following error message:\n{protein['url'][1]}"
                )
                # ... and do no attempt to scrape again
                del family.failed_pages[protein["url"][0]]

        if protein["sql"] is not None:  # Error occured when adding Protein to SQL database
            sql_failures.append(
                f"{protein['sql']} was not added to the database\t"
                f"and raised the following error when atempting to do so:\n{protein['error']}"
            )

        try:
            if protein["format"] is not None:  # Inconsistency in formating from CAZy
                format_failures.append(
                    f"Inconsistently formated. Rhe following error:\n{protein['format']}"
                )
        except KeyError:
            pass  # need to come back and see where the format key is not being added

    return family, failed_scrapes, sql_failures, format_failures


def parse_protein_table(
    protein_page_url,
    family_name,
    taxonomy_filters,
    ec_filters,
    session,
    args,
):
    """Parse proteins from the paginiation page.

    Returns generator of dicts containing error messages if any errors occured during parisng.
    Dict keys: 'url' - [url, error when connecting to CAZy], 'sql' - error when adding protein
    to SQL db, 'format' - unexpected format of data in the CAZy HTML data.

    :param protein_page_url, str, URL to the CAZy family page containing protein records
    :param family_name: str, name of CAZy family
    :param taxonomy_filters: set of genera, species and strains to restrict the scrape to
    :param ec_filters: set of EC numbers to limit scrape to
    :param session: open SQL database session
    :param args: cmd-line args parser

    Return generator object.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Retrieving proteins from {protein_page_url}")

    # connect to page
    protein_page, error = get_page(protein_page_url, args, max_tries=args.retries)
    if protein_page is None:
        logger.warning(
                f"Could not connect to {protein_page_url} after {args.retries} attempts\n"
                f"The following error was raised:\n{error}\n"
                "No protein records from this page will be retried."
        )
        return {"url": [protein_page_url, error], "sql": None, "format": None}, session

    # retrive the table of proteins
    cazyme_table = protein_page.select("table")[1]

    # check page was returned correctly. Sometimes CAZy returns the page but the table of CAZymes
    # is not populated with proteins. This often occurs after multiple timedout conenctions to CAZy
    if len(cazyme_table.select("tr")) == 1:
        # double check it is not a deleted family
        fam_activities = protein_page.select("table")[0].select("tr")[0].select("td")[0].contents
        if type(fam_activities) is list:
            if fam_activities[0] == "Deleted family!":
                return {"url": None, "sql": None, "format": "Deleted family in CAZy"}, session
        else:
            if fam_activities == "Deleted family!":
                return {"url": None, "sql": None, "format": "Deleted family in CAZy"}, session

        # check protein total to be sure there should be proteins
        data = protein_page.find_all("div", {"class": "pos_choix"})

        try:
            protein_total = int(re.findall(
                r"all \(\d+\)", data[0].text, flags=re.IGNORECASE,
            )[0].split("(")[1][:-1])
        except IndexError:
            warning = (
                f"Did not retrieve protein table populated with proteins for {protein_page_url}\n"
                "and can't retrieve protein total for family, therefore, can'y check if an"
                "incomplete page was returned."
            )
            logger.warning(warning)
            return {"url": None, "sql": None, "format": warning}, session

        if protein_total == 0:
            logger.warning(f"{protein_page_url} is an empty family in CAZy")
            return (
                {"url": None, "sql": None, "format": "Empty family in CAZy, family added to db"},
                session,
            )

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
                        f"No protein records from this page will be retried."
                    )
                )
                return {"url": [protein_page_url, error], "sql": None, "format": None}, session

    cazyme_table = protein_page.select("table")[1]

    tax_kingdom = ''  # Archaea, Bacteria, Eukaryota, Viruses, unclassified
    for row in cazyme_table.select("tr"):
        try:
            if (row.attrs["class"] == ['royaume']) and (row.text.strip() != 'Top'):
                # Row defines the taxonomy Kingdom
                tax_kingdom = row.text.strip()
                continue  # row does not contain protein
            else:  # result when row containing 'Top', becuase reached end of the page
                continue  # row does not contain protein
        except KeyError:
            pass

        if ('class' not in row.attrs) and ('id' not in row.attrs):  # row contains protein data
            yield row_to_protein(
                row,
                family_name,
                taxonomy_filters,
                tax_kingdom,
                ec_filters,
                session,
                args,
            )


def parse_protein_table_dict(
    family,
    protein_page_urls,
    taxonomy_filters,
    ec_filters,
    failed_scrapes,
    session,
    args,
):
    """Parse proteins from the paginiation page when building a CAZy dict.

    Returns generator of dicts containing error messages if any errors occured during parisng.
    Dict keys: 'url' - [url, error when connecting to CAZy], 'sql' - error when adding protein
    to SQL db, 'format' - unexpected format of data in the CAZy HTML data.

    :param protein_page_url, str, URL to the CAZy family page containing protein records
    :param family: family class instance
    :param taxonomy_filters: set of genera, species and strains to restrict the scrape to
    :param ec_filters: set of EC numbers to limit scrape to
    :param session: open SQL database session
    :param args: cmd-line args parser

    Return generator object.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Retrieving proteins from {family.name}")

    for url in tqdm(protein_page_urls, desc=f"Parsing {family.name} HTML pages"):
        protein_page, error = get_page(url, args, max_tries=args.retries)

        if protein_page is None:
            warning = (
                f"Could not connect to {url} after 10 attempts\n"
                "The following error was raised:\n"
                f"{error}\n"
                f"No protein records from this page will be retrieved."
            )
            logger.warning(warning)
            family, failed_scrapes = parse_url_error_all(url, family, failed_scrapes, warning, args)
            continue

        # Get the table on the page corresponding to CAZymes
        try:
            cazyme_table = protein_page.select("table")[1]
        except (AttributeError, IndexError):
            warning = f"No protein table found on page {url}"
            logger.warning(warning)
            family, failed_scrapes = parse_url_error_all(url, family, failed_scrapes, warning, args)
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

            try:
                protein_total = int(re.findall(
                    r"all \(\d+\)", data[0].text, flags=re.IGNORECASE,
                )[0].split("(")[1][:-1])
            except IndexError:
                warning = (
                    f"Did not retrieve protein table populated with proteins for {url}\n"
                    "and can't retrieve protein total for family, therefore, cannot check if an"
                    "incomplete page was returned."
                )
                logger.warning(warning)
                return {"url": None, "format": warning}, session

            if protein_total == 0:
                logger.warning(f"{url} is an empty in CAZy")
                failed_scrapes.append(f"{family.name} is an empty family in CAZy")
                continue

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
                family, failed_scrapes = parse_url_error_all(
                    url, family, failed_scrapes, warning, args,
                )
                continue

                cazyme_table = protein_page.select("table")[1]
                if len(cazyme_table.select("tr")) == 1:
                    warning = (
                        f"Previously connected to {url} but complete page was not returned.\nAdded "
                        f"{family.name} to be rescraped later"
                    )
                    logger.warning(warning)
                    family, failed_scrapes = parse_url_error_all(
                        url, family, failed_scrapes, warning, args,
                    )
                    continue

        for row in cazyme_table.select("tr"):
            # check if row contains protein data
            if ('class' not in row.attrs) and ('id' not in row.attrs):  # row contains protein data
                new_report_dict, session = row_to_protein_in_dict(
                    row, family.name, taxonomy_filters, ec_filters, session
                )

                if new_report_dict["error"] is not None:
                    failed_scrapes.append(
                        f"Protein {new_report_dict['protein']} raised the following error "
                        f"when being scrapped:\n{new_report_dict['error']}"
                    )

    return family, failed_scrapes, session


def parse_url_error_all(url, family, failed_scrapes, error_message, args):
    """Parse connection errors when trying to parse protein table pages.

    :param url: str, url of HTML webpage containing the current working protein table page
    :param family: Family class instance, represents a CAZy family
    :param failed_scrapes: list of erros raised when trying to scrape CAZy
    :param error_message: str, error raised when trying to scrape CAZy
    :param args: cmd-line args parser

    Return failed_scrapes and Family class instance"""
    try:
        family.failed_pages[url]

    except KeyError:  # first failed attempt for the family:kingdom
        family.failed_pages[url] = 1

    if family.failed_pages[url] >= (args.retries + 1):
        # Reached maximum attempts number of attempted connections ...
        failed_scrapes.append(
            f"{url}\t{family.cazy_class}\t"
            f"Failed to connect to this page of proteins for {family.name}\t{error_message}"
        )
        # ... and do no attempt to scrape again
        del family.failed_pages[url]

    return family, failed_scrapes
