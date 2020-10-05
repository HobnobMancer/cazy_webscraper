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

import mechanicalsoup

# create browser object
browser = mechanicalsoup.Browser()

# page to start browser at: the CAZy homepage
base_url = "http://www.cazy.org"

# create response object
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

# Note during development: The keys of interest
# Glycoside Hydrolases
# GlycosylTransferases
# Polysaccharide Lyases
# Carbohydrate Esterases
# Auxiliary Activities
# Carbohydrate-Binding Modules
