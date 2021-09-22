#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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
"""Module for crawling through the CAZy website and parsing HTML from the CAZy website."""


import logging
import re
import sys
import time

from socket import timeout
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from tqdm import tqdm
from requests.exceptions import ConnectionError, MissingSchema
from urllib3.exceptions import HTTPError, RequestError

import mechanicalsoup

from scraper.sql import sql_interface


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


class Family:
    """A single CAZy family.

    Used to keep track if family needs to be scraped again.
    """

    def __init__(self, name, cazy_class, url, download_url, path, population_size=None):
        self.name = name
        self.cazy_class = cazy_class
        self.url = url  # path to the cazy family page
        self.download_url = download_url  # url to download the family txt file
        self.path = path  # Path to downloaded family file
        if population_size is None:
            self.population_size = None
        else:
            self.population_size = population_size

    def __str__(self):
        return f"CAZy family {self.name}"

    def __repr__(self):
        return f"<Family: {id(self)}: {self.name}>"


def get_cazy_classes(cazy_home_url, excluded_classes, cazy_dict, args):
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
        sys.exit(1)

    # retrieve the h3 elements with class spip
    h3_spip_elements = homepage.find_all("h3", {"class": "spip"})

    # retrieve the div section containing the h3 element for Enzyme classes catalgoued by CAZy
    try:
        enzyme_classes_div = [
            _ for _ in h3_spip_elements if (
                str(_.contents[0].strip()).replace(u'\xa0', ' ')
                ) == 'Enzyme Classes currently covered'][0].parent

        # Retreive the enzyme class page URLs suffixs
        enzyme_class_urls  = [
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

    except AttributeError:
        logger.error(
            (
                "Attribute error raised during retrieving of CAZy class URLs.\n"
                "Therefore, cannot scrape CAZy classes, or families\n"
                "Terminating program."
            ),
            exc_info=1,
        )
        sys.exit(1)

    # compile the full CAZy class URLs from the homepage url and class suffixes

    if len(enzyme_class_urls) == 0 and len(associated_module_urls) == 0:
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
    cazy_class_urls = enzyme_class_urls + associated_module_urls
    cazy_classes = []

    for url in cazy_class_urls:
        # retrieve class name and standardise it
        class_name = url[20:-5]
        for key in cazy_dict:
            if class_name in cazy_dict[key]:
                class_name = key

        cazy_class = CazyClass(class_name, url, 0)
        cazy_classes.append(cazy_class)

    logger.info(
        "Retrieved URLs for:"
        f"{len(enzyme_class_urls)} Enzyme Classes and\n"
        f"{len(associated_module_urls)} Associated Modules classes"
    )

    return cazy_classes


def get_cazy_family_urls(class_name, class_url, cazy_home_url, cache_dir, args):
    """Retrieve all protein members of each CAZy family within the given CAZy class.

    :param class_name: str, name of CAZy class
    :param class_url: str, URL to CAZy class webpage
    :param cazy_home_url: str, URL to CAZy home page
    :param cache_dir: str representing Path to dir to write out downloaded family file to
    :param args: args parser object

    Returns:
    A CazyClass isntance
    An error message when connecting to CAZy
    A List of incorrectly formated URLs
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Retrieving URLs to families under {class_name}")

    # get the html code of the class page
    class_page, error = get_page(class_url, args, max_tries=(args.retries + 1))

    if class_page is None:
        logger.error(
                f"Couldn't connect to {class_url} after {(args.retries + 1)} attempts.\n"
                f"The following error was raised:\n{error}"
        )
        return None, error, None

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
        logger.warning(f"Failed to retrieve URLs to CAZy families for {class_name}\n")
        return None, f"Failed to retrieve URLs to CAZy families for {class_name}\n", None

    # retrieve URLs to subfamilies
    if args.subfamilies is True:
        subfam_urls = get_subfamily_links(family_h3_element, cazy_home_url)

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
        family_name = url.replace(cazy_home_url, '')
        family_name = family_name.replace('.html', '')
        family_name = family_name.replace('/', '')
        
        # generate link for download family file
        # files are gzipped to reduce storage requirements
        family_download_url = f'IMG/cazy_data/{family_name}.txt'

        family_name = url[(len(cazy_home_url) + 1): -5]

        output_path = f"{cache_dir}/{family_name}.txt"

        family = Family(family_name, class_name, url, family_download_url, output_path )
        family.members = set()  # later used to store Protein members
        cazy_families.append(family)

    if len(incorrect_urls) == 0:
        incorrect_urls = None
    
    logger.info(f"Retrieved URLs for {len(cazy_families)} from {class_name} class page")

    return cazy_families, None, incorrect_urls


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


def parse_family(family, args):
    """Coordinate parsing a CAZy family.
    
    The overal process involved:
    Download text file containing family members. Retrieve number of members listed in family on
    the CAZy website (used for download validation). Parse downloaded text file and add proteins
    to the database and/or dict.

    :param family: Family class instance
    :param args: cmd-line args parser
    
    Return
    - Family class instance
    - Bool indicating successful or unsuccesful retrieval of CAZy family txt file
    """
    logger = logging.getLogger(__name__)

    # retrieve the number of proteins listed in the family by CAZy
    family_population = get_family_population(family, args)

    if family_population == 0 or \
        family_population == 'Deleted family!' or \
            family_population == 'Empty family' or \
                family_population == 'Failed Retrieval':
                logger.warning(
                    f"No proteins found for {family.name}"
                )
    
    family.population_size = family_population
    
    success = download_family_file(family, args)
    if success is False:  # failed to download file, is attempts remain retry later
        return family, success

    # tally number of proteins added to the local CAZyme db
    protein_count = 0

    with open(family.path, 'r') as fh:
        lines = fh.read().splitlines()

    for line in lines:
        print(line)

    # Scrape proteins from text file
    # Tally number of proteins
    # How many could not be added to CAZy
    # track against family population


def get_family_population(family, args):
    """Retrieve the CAZy family population from the CAZy family page.
    
    :param family: Family class instance
    :param args: cmd-line args parser
    
    Return int: number of proteins in the CAZy family
    """
    logger = logging.getLogger(__name__)

    family_page = get_page(family.url, args, max_tries=(args.retries + 1))

    # retrieve the table containing the Family data
    family_data = family_page.find_all("div", {"class": "pos_choix"})

    try:
        family_population = int(re.findall(
            r"Download \D{2,3}\d+? \(\d+?\)",
            family_data[0].text,
            flags=re.IGNORECASE,
        )[0].split("(")[1].replace(")", ""))

    except (IndexError, AttributeError):
        family_population = 0
    
    if family_population == 0:
        # check if the family has been deleted, and labeled as such under 'Activities in Family'
        try:
            family_activities_cell = family_page.select("table")[
                0].select("tr")[0].select("td")[0].contents[0].strip()

            if family_activities_cell == 'Deleted family!':
                family_population = 'Deleted family!'
                logger.warning(f"{family.name} is a deleted family in CAZy")
            else:
                logger.warning(f"{family.name} is an empty CAZy family")
                family_population = 'Empty family'
        except Exception:
            logger.warning(f"Could not retrieve family population for {family.name}")
            family_population = 'Failed Retrieval'
    
    return family_population


def download_family_file(family, args):
    """Download text file with CAZy family members from CAZy
    
    :param family: Family class instance
    :param args: cmd-line args parser
    
    Return
    - True if download is successful
    - False if tried to connected but failed, and ran out of attempts
    - False if failed to download the file successfully after making URL connection
            This indicate to try again if attempts remain (args.retries)
    - 'empty' if an empty family (population_size = 0 and HTTP 404 error)
    - 404 if not labelled as an empty family but HTTP 404 error raised
    """
    logger = logging.getLogger(__name__)


    # Try connection
    tries, response = 0, None
    while (tries < (args.retries + 1)) and (response is None):

        try:
            response = urlopen(family.download_url, timeout=args.timeout)

        except (HTTPError, URLError, timeout) as err:

            if err.code == 404:
                if family.population_size == 0 or \
                    family.population_size == 'Deleted family!' or \
                        family.population_size == 'Empty family':
                        logger.warning(
                            f"{family.name} annotated with family population: "
                            f"{family.population_size}\n"
                            f"Retrieving {family.name} text file raised: HTTP Error 404: Not Found\n"
                            f"{family.name} is an empty family. {family.name} added to the database\n"
                            f"but NO proteins from {family.name} added to the local database."
                        )
                        return 'empty'
                
                else:
                    logger.error(
                        f"{family.name} not marked as empty family but connecting to\n"
                        f"{family.download_url} raised: HTTP Error 404: Not Found\n"
                        f"No proteins from {family.name} added to the database"
                    )
                    return '404'

            logger.warning(
                f"Failed to connected to {family.url}. The following error was raised\n"
                f"{err}"
            )
            tries += 1
    
    if response is None:
        logger.warning(
            f"Failed to download {family.name} file after {args.retries} tries.\n"
            "If family tries are remaining, will retry later."
        )
        return False
    
    # Connection was successfully made
    
    bsize = 1_048_576
    try:
        with open(family.path, 'wb') as fh:
            while True:
                buffer = response.read(bsize)
                if not buffer:
                    break
                fh.write(buffer)
    except IOError as err:
        logger.warning(
            f"Failed to download {family.name}. The following error was raised:\n"
            f"{err}"
        )
        return False
    
    # check if the file was downloaded ok
    if family.path.stat().st_size == 0:
        logger.warning(
            f"File for {family.name} was not downloaded successfully.\n"
            "File size is 0.\n"
            "If retries exist, will retry later"
        )
        return False
    
    return True






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
