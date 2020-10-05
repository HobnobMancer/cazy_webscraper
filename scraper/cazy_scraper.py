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

import re

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

    # tuple of CAZy classes and abbrievations
    cazy_classes = [
        ["Glycoside Hydrolases", "GH"],
        ["GlycosylTransferases", "GT"],
        ["Polysaccharide Lyases", "PL"],
        ["Carbohydrate Esterases", "CE"],
        ["Auxiliary Activities", "AA"],
        ["Carbohydrate-Binding Modules", "CBM"],
    ]

    # Navigate through each CAZy class pages
    # iterating through all pages listing all CAZymes for each CAZy class
    index = 0

    for index in range(len(cazy_classes)):
        # retrieve full CAZy class name
        cazy_class = cazy_classes[index][0]
        # create url for CAZy class main page
        class_url = base_url + "/" + links_dict[cazy_class]

        # retrieve all links to CAZy familes for given CAZy class
        get_family_links(browser, class_url, cazy_classes[index][1])
        index += 1
    # url = "http://www.cazy.org/Glycoside-Hydrolases.html"
    # get_family_links(browser, url)


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


def get_family_links(browser, class_url, class_abbreviation):
    """Navigate CAZy class page, iterating through pages listing CAZymes.
    
    Return list of links to CAZy family pages for given class.
    """
    class_page = browser.get(class_url)

    # obtain links on class main page
    all_links = class_page.soup.select("a")
    # empty to store all links from CAZy class main page
    family_links = []
    # search pattern to determine if link is for CAZy family or not
    pattern = re.compile(rf"{class_abbreviation}\d+?.*?\.html")

    # retieve all links from CAZy class main page
    for link in all_links:
        try:
            address = link["href"]
        except KeyError:
            pass
        # check if link is for CAZy family page
        search_result = re.match(pattern, address)
        if search_result:
            family_links.append(address)

    print(len(family_links))

    return family_links


if __name__ == "__main__":
    main()
