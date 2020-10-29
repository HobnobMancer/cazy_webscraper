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

"""

from collections import defaultdict

import mechanicalsoup

from tqdm.notebook import tqdm


class Protein:
    """A single protein.

    Each protein has a name, source organism (source), and links to external databases. The links to
    external databases are stored in a dictionary, keyed by the external database name ('str') with
    'list' values becuase there may be multiple links per database.
    """

    def __init__(self, name, source, links=None):
        self.name = name
        self.source = source
        if links is None:
            self.links = defaultdict(list)
        else:
            self.links = links

    def __str__(self):
        """Create representative string of class object"""
        return f"{self.name} ({self.source}): links to {self.links.keys()}"

    def __repr__(self):
        """Create representative object"""
        return(
            (
                f"<Protein: {id(self)}: {self.name}, "
                f"({self.source}), {len(self.links)} to external databases>"
            )
        )


class Family:
    """A single CAZy family.

    Each family has a name and a set of proteins that are members of the family.
    """

    members = set()  # holds Protein instances

    def __init__(self, name):
        self.name = name

    def __str__(self):
        return f"CAZy family {self.name}: {len(self.members)} protein members"

    def __repr__(self):
        return f"<Family: {id(self)}: {self.name}, {len(self.members)} protein members"


def main():
    """Set up parser, logger and coordinate overal scrapping of CAZy."""
    cazy_home = "http://www.cazy.org"  # the CAZy homepage URL

    # retrieve links to CAZy class pages
    class_pages = get_cazy_class_pages(cazy_home)

    # retrieve links to CAZy family pages
    for class_url in class_pages:
        family_urls = get_cazy_family_pages(class_url, cazy_home, subfam_retrieval)
        for family_url in family_urls:


def get_cazy_class_pages(cazy_home):
    """Returns a list of CAZy class main/home page URLs for each specified class as the CAZy site.

    :param cazy_url: str, URL to the CAZy home page.

    Return list of URLs.
    """
    # define items to be excluded from returned class list, ALWAYS exlide links to genomes
    exclusions = ("<strong>Genomes</strong>")    
    # exclusions.append(config file input)

    # scrape the home page
    home_page = get_page(cazy_home)

    return [f"{cazy_home}/{_['href']}" for _ in home_page.soup.find_all("a", {"class": "spip_out"})
            if (not _["href"].startswith("http")) and (str(_.contents[0]) not in exclusions)]


def get_cazy_family_pages(class_url, cazy_home, subfam_retrieval):
    """Retrieve all protein members of each CAZy family within the given CAZy class.

    :param class_url: str, URL to the CAZy class
    :param cazy_home: str, URL to CAZy home page
    :param subfam_retrieval: bool, enables subfamily retrieval if true

    Returns list of URLs to family pages.
    """
    # scrape the class page
    class_page = get_page(class_url)

    # retrieve the <h3> element that titles the div section containing the tables of family links
    family_h3_element = [_ for _ in class_page.soup.find_all("h3"), {"class": "spip"} if
                         str(_.contents[0]) == "Tables for Direct Access"][0]

    # retrieve all tables within the parent div section of the <h3> element
    tables = family_h3_element.parent.find_all("table")

    # tables[0] is the table containing links to CAZy families
    # tables[1] is the table containing the link to unclassified proteins

    family_urls = family_urls = [f"{cazy_home}/{_[;'href']}" for _ in tables[0].find_all("a")]
    family_urls.append(f"{cazy_home}/{tables[1].a['href']}")
    if subfam_retrieval:
        family_urls.append(get_subfamily_links(family_h3_element, cazy_home))

    return family_urls


def get_subfamily_links(family_h3_element, cazy_home):
    """Retrieve URL links to CAZy subfamilies.

    :param family_h3_element: bs4.element.Tag, h3 element titling the page div
    :param cazy_home: str, URL to CAZy home_page

    Return list of URLs to subfamilies.
    """
    all_links = family_h3_element.find_all("a")

    pattern = re.compile(rf"\D+?\d+?_\d+?\.html")

    urls = []  # empty list to store subfamily URLs

    for link in all_links:
        try:
            search_result = re.search(pattern, link["href"])
            urls.append(f"{cazy_home}/{link['href']}")
        except KeyError, AttributeError as error:
            # KeyError raised if link does not have ['href']
            # AttributeError error raised search_result is None becuase not subfam link
            pass

    return urls


def parse_family_page(family_url, cazy_home):
    """Retrieve all protein records for given CAZy family.

    Protein records are listed in a pagination method, with 1000 proteins per page.

    :param family_url: str, URL to CAZy family main page
    :param cazy_home: str, URL to CAZy home page

    Return list of protein records.
    """
    # compile URL to first family page of protein records
    first_pagniation_url = family_url.replace(".html", "_all.html")
    first_pagination_page = get_page(first_pagniation_url)

    # retrieve the URL to the final page of protein records in the pagination listing
    last_pagination_url = first_pagination_page.find_all("a", {"class": "lien_pagination", "rel": "nofollow"})[-1]

    url_prefix = last_pagination_url["href"].split("PRINC=")[0] + "PRINC="
    last_princ_num = int(last_pagination_url["href"].split("PRINC=")[-1].split("#pagination")[0])
    url_suffix = "#pagination" + last_pagination_url["href"].split("#pagination")[-1]

    # Build list of urls to all pages in the pagination listing, increasing the PRINC increment
    protein_page_urls = [f"{cazy_home}/{prefix}{_}{suffix}" for _ in range(1000, last_princ_num + 1000, 1000)"]
    protein_page_urls.insert(0, first_pagniation_url)

    proteins = []
    for url in protein_page_urls:
        proteins.append(parse_proteins(url))
    
    # MISSING RETURN -- DON'T FORGET TO ADD LATER!!


def parse_proteins(protein_page_url):
    """Returns generator of Protein objects for all protein rows on a single CAZy family page.

    :param protein_page_url, str, URL to the CAZy family page containing protein records

    Return generator object.
    """
    protein_page = get_page(url)

    # retrieve protein record table
    protein_table = protein_page.soup.find_all("table", {"class": "listing"})[0]
    protein_rows = [_ for _ in protein_table.descendants if (_.name == "tr") and
                    ("id" not in _.attrs) and ("class" not in _.attrs)]

    for row in rows:
        yield row_to_protein(row)


def row_to_protein(row):
    """Returns a Protein object representing a single protein row from a CAZy family protein page.

    Each row, in order, contains the protein name, EC number, source organism, GenBank ID(s), 
    UniProt ID(s), and PDB accession(s).

    :param row: tr element from CAZy family protein page

    Return Protein instance.
    """
    # retrieve list of cells ('td' elements) in row
    tds = list(row.find_all("td"))

    protein_name = tds[0].contents[0].strip()
    ec_number = tds[1].contents[0].strip()
    source_organism = tds[2].a.get_text()
    links = {}

    # test for len(tds[x].contents) in case there is no link,
    # the check of .name then ensures link is captured
    if len(tds[3].contents) and tds[3].contents[0].name == "a":
        links["GenBank"] = [f"{_.get_text()}: {_['href']}" for _ in tds[3].contents if _.name == "a"]
    if len(tds[4].contents) and tds[4].contents[0].name == "a":
        links["UniProt"] = [f"{_.get_text()}: {_['href']}" for _ in tds[4].contents if _.name == "a"]
    if len(tds[5].contents) and tds[5].contents[0].name == "a":
        links["PDB"] = [f"{_.get_text()}: {_['href']}" for _ in tds[5].contents if _.name == "a"]

    return Protein


def browser_decorator(func):
    """Decorator to retry the wrapped function up to 'retries' times."""

    def wrapper(*args, retries=10, **kwargs):
        tries, success = 0, False
        while not success and (tries < retries):
            response = func(*args, **kwargs)
            if str(response) == "<Response [200]>":  # response was successful
                success = True
            # if response from webpage was not successful
            tries += 1
        if not success:
            print("Ran out of connection retries, and was unable to connect")
            return None
        else:
            return response

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

    return page


if __name__ == "__main__":
    main()
