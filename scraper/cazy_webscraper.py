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

import re
import sys

import json
import mechanicalsoup

from tqdm import tqdm

from scraper.utilities import build_logger

class Site:
    """A single website, parent to multiple Pages."""

    pages = set()  # enable storing (web)pages within this site

    def __init__(self, base_url):
        self.base_url = base_url
        # base url is the url of the homepage of the website

    def add_page(self, page):
        """Add (web)page to this (web)Site.

        Add all pages, including specific page types such as class, family,
        and protein pages.

        :param page: Page instance to add site
        """
        if page.url not in self.page_urls:
            # do not add duplicate links
            self.pages.add(page)

    def get_class(self, cazy_class):
        """Return ClassPage for passes CAZy class from this Site.

        ClassPage is the summary page of a CAZy class.

        :param cazy_class: str, name of CAZy class.
        """
        for page in self.class_pages:
            if page.cazy_class == cazy_class:
                # once the correct Class page is found
                return page

    def get_class_url(self, cazy_class):
        """Return the URL to the passed CAZy class summary page."""
        for page in self.class_pages:
            if page.cazy_class == cazy_class:
                # once the correct Class page is found
                return page.url

    def get_family_urls(self, cazy_class):
        """Return a list of URLs to all family summary pages for the passed CAZy class."""
        return [_.url for _ in self.family_pages if (_.cazy_class == cazy_class)]

    def get_table_urls(self, cazy_family):
        """Return a list of URLs to all protein table pages for given CAZy family."""
        return [_.url for _ in self.table_pages if (_.cazy_family == cazy_family)]

    def __repr__(self):
        """Create represent instance."""
        return f"<Site: {id(self)}, base_url: {self.base_url}, page_count: {self.page_count}>"

    @property
    def class_pages(self):
        """Return ClassPages for this site."""
        return [_ for _ in self.pages if isinstance(_, ClassPage)]

    @property
    def family_pages(self):
        """Return FamilyPages for this site."""
        return [_ for _ in self.pages if isinstance(_, FamilyPage)]

    @property
    def table_pages(self):
        """Return ProteinTablePages for this site."""
        return [_ for _ in self.pages if isinstance(_, ProteinTablePage)]

    @property
    def page_urls(self):
        """Return list of all URLs for each page in this Site."""
        return [_.url for _ in self.pages]

    @property
    def page_count(self):
        """Return the total number of (web)pages in this (web)site."""
        return len(self.pages)


class Page:
    """Describes the relevant features of a (web)page belonging to the parent CAZy (web)Site."""

    links = set()
    url: str

    def __init__(self, url):
        self.url = url

    def add_link(self, url):
        """Add the passed url to a collection of links retrieved from the current page.

        :param url: str, URL of the link to be added to the collection.
        """
        self.links.add(url)

    def __repr__(self):
        """Create represent instance."""
        return f"<Page: {id(self)}, URl={self.url}>"


class ClassPage(Page):
    """CAZy website page describing and summarising a CAZy class."""

    cazy_class: str  # full name of the CAZy class

    def __init__(self, url, cazy_class):
        self.cazy_class = cazy_class
        Page.__init__(self, url)

    def __repr__(self):
        """Create represent instance."""
        return f"<ClassPage: {id(self)}, CLASS={self.cazy_class}, URL={self.url}>"


class FamilyPage(ClassPage):
    """CAZy website page describing and summarising a CAZy family."""

    cazy_family: str  # abbreviated-CAZy-class _ family-number

    def __init__(self, url, cazy_class, cazy_family):
        self.cazy_family = cazy_family
        ClassPage.__init__(self, url, cazy_class)

    def __repr__(self):
        """Create represent instance."""
        return (
            f"<FamilyPage: {id(self)}, FAMILY={self.cazy_family}, "
            f"CLASS={self.cazy_class}, URL={self.url}>"
        )


class ProteinTablePage(FamilyPage):
    """CAZy website page containing a protein table, for a specific CAZy family."""

    table_number: str  # the number of the table page

    def __init__(self, url, cazy_class, cazy_family, table_number):
        self.table_number = table_number
        FamilyPage.__init__(self, url, cazy_class, cazy_family)

    def __repr__(self):
        """Create respresent instance."""
        return (
            f"<ProteinTablePage: {id(self)}, TABLE#={self.table_number}, "
            f"FAMILY={self.cazy_family}, CLASS={self.cazy_class}, URL={self.url}>"
        )


def main():
    """Coordinate scraping of CAZy website."""
    # build logger
    logger = build_logger()

    # page to start browsing from: the CAZy homepage
    base_url = "http://www.cazy.org"

    # create site class object
    site = Site(base_url)

    # list of CAZy classes to be scraped from CAZy website will eventually
    # be configured by curation file
    classes = [
        "Glycoside Hydrolases",
        "GlycosylTransferases",
        "Polysaccharide Lyases",
        "Carbohydrate Esterases",
        "Auxiliary Activities",
        "Carbohydrate-Binding Modules",
    ]

    cazy_class_definitions = {
        "Glycoside Hydrolases": "GH",
        "GlycosylTransferases": "GT",
        "Polysaccharide Lyases": "PL",
        "Carbohydrate Esterases": "CE",
        "Auxiliary Activities": "AA",
        "Carbohydrate-Binding Modules": "CBM",
    }

    # populate site with webpages from CAZy homepage
    get_links_from_hompage(base_url, classes, site)

    # populate site with webpages from each of the CAZy class pages
    for cazy_class in classes:
        # compile url to CAZy class web page
        url = site.get_class_url(cazy_class)
        class_url = base_url + "/" + str(url)

        # retrieve CAZy abbreviation of class name (e.g. "GH")
        class_abbrev = cazy_class_definitions[cazy_class]

        # populate site with family pages for CAZy class
        get_links_from_classpage(class_url, cazy_class, class_abbrev, site)

        # retrieve list of urls to all family pages for given class
        families = site.get_family_urls(cazy_class)

        for family in families:  # make tqdm progress bar
            # compile full url for CAZy family page
            family_url = base_url + "/" + str(family)

            # populate site with the protein pages for each CAZy family
            get_links_from_familypage(family_url, cazy_class, base_url, site)
            cazy_family = family_url[20:-5]
            print(
                f"**number of protein table pages for {family}:",
                len(site.get_table_urls(cazy_family)),
            )

            # parse protein tables for CAZy family
            # ...


def get_links_from_hompage(base_url, classes, site):
    """Populate site with links to other pages from the homepage.

    Differentiate between links/urls to pages summarising and describing
    CAZy classes (ClassPage, a subclass of Page) and pages that do not (Page).

    :param base_url: str, url of site homepage
    :param classes: list, list of CAZy classes to be scraped
    :param site: Site class object

    return nothing
    """
    # scrape CAZy homepage
    homepage = get_page(base_url)

    # check connection was successful, if not terminate
    if homepage is None:
        print("could not connect to homepage, terminating")
        sys.exit(1)

    # retrieve all links from homepage
    all_links = homepage.soup.select("a")

    # create empty list to store all links to other pages in
    pages = []
    for link in all_links:
        try:
            if link.text in classes:
                # if link is for a CAZy class summary page, create a ClassPage class object
                pages.append(ClassPage(link["href"], link.text))
            else:
                # if links is not for a CAZy class page, create general Page class object
                pages.append(Page(link["href"]))
        except KeyError:
            pass

    # populate site object with all (web)page objects retrieved from the homepage
    for page in pages:
        site.add_page(page)

    return


def get_links_from_classpage(class_url, cazy_class, class_abbrev, site):
    """Add url links from CAZy class page to the Site.

    Differentiate between link/url to a Family page and non-family
    pages.

    :param class_url: str, url to the page of the given cazy_class
    :param cazy_class: str, current working CAZy class
    :param class_abbrev: str, CAZy-standard abbreviation of CAZy class
    :param site: Site class object

    Return nothing.
    """
    # scrape CAZy class page
    class_page = get_page(class_url)

    # check connection was successful, if return nothing
    if class_page is None:
        print(
            (
                f"could not connect to {cazy_class} page, therefore, "
                "retrieving not family pages for f{cazy_class}"
            )
        )
        return

    # obtain links on CAZy class page
    all_links = class_page.soup.select("a")

    # search patter to determine if a link is for a CAZy family page, describing a CAZy famoly
    pattern = re.compile(rf"{class_abbrev}\d+?.*?\.html")

    # create empty list to store all links to other pages in
    pages = []

    for link in all_links:
        # check if link is to the Family page describing the CAZy family
        try:
            search_result = re.search(pattern, link["href"])

            try:
                pages.append(
                    FamilyPage(link["href"], cazy_class, search_result.group()[:-5])
                )
                # if subfam retrieval is false, check if link is for subfamily
                # subfam_index = link["href"].find("_")
                # if subfam_index != -1:  # link not for a sub-family page
                # pages.append(FamilyPage(link["href"], cazy_class, search_result[:-5]))
                # else: # if retrieving family and subfamily pages

            except AttributeError:
                # Raised if search_result is None because CAZy family not found in link.text
                pages.append(Page(link["href"]))

        except KeyError:
            pass

    # populate site object with all (web)page objects retrieved from the working CAZy class page
    for page in pages:
        site.add_page(page)

    return


def get_links_from_familypage(family_url, cazy_class, base_url, site):
    """Add url links from CAZy family page to the Site.

    Differentiate between link/url to a Protein page and non-protein
    pages.

    :param family_url: str, url to the page of the given cazy_class
    :param cazy_class: str, current working CAZy class
    :param base_url: str, base url of CAZy website
    :param site: Site class object

    Return nothing.
    """
    # retrieve name of family from url
    cazy_family = family_url[20:-5]

    # scrape CAZy family page
    family_page = get_page(family_url)

    # obtain all links on family page
    all_links = family_page.soup.select("a")

    # create empty lists to store all links to other pages
    pages = []

    # populate site with pages
    for link in all_links:
        try:
            url = link["href"]

            if link.text == "All":
                # This is the link to the first protein table page for family
                # compile URL for first protein table page
                page_url = url
                # populate site with all protein table pages for given CAZy family
                get_protein_table_pages(page_url, cazy_family, cazy_class, site)

            else:
                # link is to a page not containing a protein table
                pages.append(Page(link["href"]))

        except KeyError:
            pass

    # populate site object with all (web)page objects retrieved from the working CAZy family page
    for page in pages:
        site.add_page(page)

    return


def get_protein_table_pages(first_page_url, cazy_family, cazy_class, site):
    """Add url links to pages on all protein table pages for given CAZy family to the Site.

    Differentiate between links to protein table pages and those not to protein table pages. Finds
    the highest protein table page number. Use maximum number of pages to direct retrieval of pages
    containing protein tables.

    Protein table pages are numbers in the same manner that google search result pages are. Each
    protein table page is called simply after the number of the page in the series of protein table
    pages.

    :param first_page_url: str, url to the first protein table page for given family
    :param cazy_family: str, name of current working CAZy family
    :param cazy_class: str, name of current working CAZy class
    :param site: Site class object

    Return nothing.
    """
    print(f"getting protein table pages for {cazy_family}")
    # scrape the first protein table page
    first_table_page = get_page(first_page_url)
    if first_table_page is None:
        print(
            (
                f"failed to connect to first protein table page for {cazy_family}.\n"
                f"Retrieving not protein table pages for {cazy_family}"
            )
        )
        return

    # obtain all links on protein table page
    all_links = first_table_page.soup.select("a")

    # create empty list to store links from pages in
    pages = []  # to store all pages
    protein_table_numbers = []  # store number of protein page

    # Determine the total number of protein table pages for the given CAZy family
    # Protein table pages are identifiable by being called the number of the protein table, e.g. '2'
    for link in all_links:
        try:
            url = link["href"]

            # test if link name is a interger, indicating protein table page
            try:
                page_name = int(link.text)
                protein_table_numbers.append(page_name)

            except ValueError:
                # Raised if link.text contains non-digit characters in name, thus not table page
                pages.append(Page(url))

        except KeyError:
            pass

    print("Number of pages before adding protein table pages:", len(pages))
    if len(protein_table_numbers) == 0:
        # no links to other protein table pages retrieved, only one protein table page
        pages.append(ProteinTablePage(first_page_url, cazy_class, cazy_family, "1"))
        print(
            "Only 1 protein table page\nNumber of pages after adding protein table pages:",
            len(pages),
        )

    else:
        # Retrieve the highest protein table page number, finding the total number of protein table
        # pages for given protein family
        protein_table_numbers.sort(reverse=True)
        total = protein_table_numbers[0]
        print("total number of protein table pages=", total)
        protein_table_total = int(total) + 1  # ensure capturing final table page

        page_count = 1  # the number of the protein table page of interest

        page_url = first_page_url

        while page_count < protein_table_total:
            if page_count == 1:
                # add first protein table page to list
                pages.append(
                    ProteinTablePage(page_url, cazy_class, cazy_family, page_count)
                )
                page_count = 2  # search for second protein table page

                # add second protein table page from 'all_links' from first protein table page
                for link in all_links:
                    try:
                        url = link["href"]
                        if link.text == str(page_count):
                            page_url = (
                                first_page_url[:20] + url
                            )  # compile ful URL for page
                            pages.append(
                                ProteinTablePage(
                                    page_url, cazy_class, cazy_family, page_count
                                )
                            )
                    except KeyError:
                        pass

                page_count = 3  # first and second pages have been added to [pages]

            else:
                # scrape the next protein table page
                new_page = get_page(page_url)

                # retrieve all links from protein table page
                all_links = new_page.soup.select("a")

                for link in all_links:
                    try:
                        url = link["href"]

                        try:
                            if link.text == str(page_count):
                                page_url = (
                                    first_page_url[:20] + url
                                )  # compile full URL for page
                                pages.append(
                                    ProteinTablePage(
                                        page_url,
                                        cazy_class,
                                        cazy_family,
                                        str(page_count),
                                    )
                                )
                        except ValueError:
                            pages.append(Page(link["href"]))

                    except KeyError:
                        pass

                page_count += 1

        print(f"page count= {page_count}, number of pages=", len(pages))

    for page in pages:
        site.add_page(page)

    return


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
