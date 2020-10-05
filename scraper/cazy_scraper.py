#!/usr/bin/env python
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

import pandas as pd

from urllib.request import urlopen

import mechanicalsoup


def main():
    # page to start browser at: the CAZy homepage
    base_url = "http://www.cazy.org"

    # create browser object
    browser = mechanicalsoup.Browser()

    # Retrieve all links from the CAZy homepage
    links_dict = get_links(browser, base_url)
    print("got links dictionary")

    # list of CAZy classes
    cazy_classes = [
        "Glycoside Hydrolases",
        "GlycosylTransferases",
        "Polysaccharide Lyases",
        "Carbohydrate Esterases",
        "Auxiliary Activities",
        "Carbohydrate-Binding Modules",
    ]

    # Navigate through each CAZy class pages
    # iterating through all pages listing all CAZymes for each CAZy class
    # for cazy_class in cazy_classes:
    # create url for CAZy class main page
    # class_url = base_url + "/" + links_dict[cazy_class]
    # navigate through pages associated for each CAZy class
    class_url = "http://www.cazy.org/Glycoside-Hydrolases.html"
    navigate_class_page(browser, class_url)


def get_links(browser, base_url):
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


def navigate_class_page(browser, class_url):
    """Navigate CAZy class page, iterating through pages listing CAZymes."""
    class_page = browser.get(class_url)

    # obtain links on class main page
    all_links = class_page.soup.select("a")
    # empty to store all links from CAZy class main page
    all_links = []

    # retieve all links from CAZy class main page
    for link in all_links:
        try:
            family_links.append(link)
        except KeyError:
            pass


if __name__ == "__main__":
    main()
