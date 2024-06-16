#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
# (c) James Hutton Institute 2024
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


import argparse
import logging

from cazy_webscraper.sql import sql_interface, sql_orm
from cazy_webscraper.utilities import termcolour


logger = logging.getLogger(__name__)


def add_main_scrape_message(
    kingdom_filters: set[str],
    taxonomy_filters: set[str],
    taxonomy_filter_dict: dict,
    time_stamp: str,
    config_dict: dict,
    args: argparse.Namespace,
    connection
):
    """add information of scraping CAZy to the local CAZyme database"""
    scrape_config_message = (
        "Configuration:\n"
        f"Classes to scrape: {config_dict['classes']}\n"
        f"GH fams to scrape: {config_dict['Glycoside Hydrolases (GHs)']}\n"
        f"GT fams to scrape: {config_dict['GlycosylTransferases (GTs)']}\n"
        f"PL fams to scrape: {config_dict['Polysaccharide Lyases (PLs)']}\n"
        f"CE fams to scrape: {config_dict['Carbohydrate Esterases (CEs)']}\n"
        f"AA fams to scrape: {config_dict['Auxiliary Activities (AAs)']}\n"
        f"CBM fams to scrape: {config_dict['Carbohydrate-Binding Modules (CBMs)']}\n"
        f"Scraping subfamilies: {args.subfamilies}"
    )
    scrape_config_message += "\nTaxonomy filters applied." if len(taxonomy_filters) != 0 else ""
    scrape_config_message += f"\nScraping only tax kingdoms: {kingdom_filters}" if len(kingdom_filters) < 5 else ""
    
    logger.info(termcolour(scrape_config_message, "cyan"))

    logger.info("Adding log of scrape to the local CAZyme database")
    with sql_orm.Session(bind=connection) as session:
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            kingdom_filters,
            taxonomy_filter_dict,
            set(),  # ec_filters not applied when scraping CAZy
            'CAZy',
            'CAZy annotations',
            session,
            args,
        )
