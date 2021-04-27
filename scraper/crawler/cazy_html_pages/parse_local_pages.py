#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
"""Scrape data from local HTML files"""


import json
import logging
import os
import re
import sys

import pandas as pd

from datetime import datetime
from pathlib import Path
from tqdm import tqdm

from bs4 import BeautifulSoup

from scraper.crawler import row_to_protein
from scraper.utilities import build_logger


def parse_local_pages(
    args,
    cazy_home,
    start_time,
    time_stamp,
    session,
    taxonomy_filters,
    ec_filters,
):
    """Scrape local HTML pages to build CAZy database.

    :param args: cmd-line args parser
    :param start_time: pd object, time program was initiated
    :param time_stamp: str, time program was initiated
    :param cazy_home: str, cazy homepage URL address
    :param session: open SQL database session

    Return nothing.
    """
    if args.output is not sys.stdout:
        out_log_path = args.output
    else:
        out_log_path = None

    sql_failures_logger = build_logger(out_log_path, f"SQL_errors_CW_{time_stamp}.log")
    logger = logging.getLogger(__name__)

    html_files = get_html_files(args)

    for html_file_path in tqdm(html_files, desc="Parsing HTML files"):
        # open the file
        with open(html_file_path) as fp:
            page = BeautifulSoup(fp, features="lxml")

        # get the table of CAZymes
        try:
            cazyme_table = page.select("table")[1]
        except IndexError:
            logger.warning(
                f"No CAZyme table found in {html_file_path}.\n"
                "Retrieving not proteins from this file"
            )
            continue

        # get family name from file path
        try:
            family_name = re.search(
                (
                    r"GH(\d+_all|\d+_\d+_all)|GT(\d+_all|\d+_\d+_all)|PL(\d+_all|\d+_\d+_all)"
                    r"|CE(\d+_all|\d+_\d+_all)|AA(\d+_all|\d+_\d+_all)|CBM(\d+_all|\d+_\d+_all)"
                ),
                str(html_file_path),
            ).group()
            family_name = family_name.replace("_all", "")
        except (AttributeError, TypeError):
            logger.warning(
                f"Incorrect formating of file path {html_file_path}\n"
                "Not retrieing CAZyme from file"
            )
            continue

        # check if a deleted and/or empty family
        try:
            if len(cazyme_table.selected("tr")) == 1:
                activities = page.select("table")[0].select("tr")[0].select("td")[0].contents
                if type(activities) is list:

                    if activities[0] == 'Deleted family!':
                        logger.warning(
                            f"{family_name} from {html_file_path} is a deleted Family in CAZy."
                        )
                        continue
                    else:  # check if got any proteins
                        data = page.find_all("div", {"class": "pos_choix"})
                        try:
                            protein_total = int(re.findall(
                                r"all \(\d+\)", data[0].text, flags=re.IGNORECASE,
                            )[0].split("(")[1][:-1])
                        except IndexError:
                            logger.warning(
                                f"Could not retrieve the number of proteins for {html_file_path}.\n"
                                "Therefore, could not verify page was properly returned.\n"
                                "Will scrape to see if there are proteins present."
                            )
                            protein_total = 1000  # max num that can appear

                        if protein_total == 0:
                            logger.warning(
                                f"{family_name} from {html_file_path} is not listed as a deleted "
                                "Family in CAZy but is an empty."
                            )
                            continue
                        else:
                            logger.error(
                                f"{html_file_path} protein table does not include any CAZymes"
                            )
                            continue

                else:
                    if activities == 'Deleted family!':
                        logger.warning(
                            f"{family_name} from {html_file_path} is a deleted Family in CAZy."
                        )
                        continue
                    else:  # check if a protein total is included and not 0
                        data = page.find_all("div", {"class": "pos_choix"})
                        try:
                            protein_total = int(re.findall(
                                r"all \(\d+\)", data[0].text, flags=re.IGNORECASE,
                            )[0].split("(")[1][:-1])
                        except IndexError:
                            logger.warning(
                                f"Could not retrieve the number of proteins for {html_file_path}.\n"
                                "Therefore, could not verify page was properly returned.\n"
                                "Will scrape to see if there are proteins present."
                            )
                            protein_total = 1000  # max num that can appear

                        if protein_total == 0:
                            logger.warning(
                                f"{family_name} from {html_file_path} is not listed as a deleted "
                                "Family in CAZy but is an empty."
                            )
                            continue
                        else:
                            logger.error(
                                f"{html_file_path} protein table does not include any CAZymes"
                            )
                            continue
        except TypeError:
            logger.warning(f"no CAZyme table in {html_file_path}\nNot retrieving CAZyme from {family_name}")
            continue

        # check if an 'all' or 'kingdom' page:
        path_ = str(args.scrape_files)
        filename = html_file_path.replace(path_, "")

        if str(filename).find("_all") != -1:  # scraping an 'all' page
            for row in tqdm(
                cazyme_table.select("tr"), desc=f"Parsing proteins in {html_file_path}",
            ):
                try:
                    if (row.attrs["class"] == ['royaume']) and (row.text.strip() != 'Top'):
                        # Row defines the taxonomy Kingdom
                        tax_kingdom = row.text.strip()
                        continue  # row does not contain protein
                    else:  # result when row containing 'Top', becuase reached end of the page
                        continue  # row does not contain protein
                except KeyError:
                    pass

                if ('class' not in row.attrs) and ('id' not in row.attrs):  # row contains protein
                    report, session = row_to_protein(
                        row,
                        family_name,
                        taxonomy_filters,
                        tax_kingdom,
                        ec_filters,
                        session,
                        args,
                    )
                    if report["sql"] is not None:
                        sql_failures_logger.warning(report["sql"])

        else:
            # scraping a kingdom page so retrieve the kingdom from the file name
            if filename.find("bacteria") != -1:
                kingdom = "Bacteria"
            elif filename.find("archaea") != -1:
                kingdom = "Archaea"
            elif filename.find("eukaryota") != -1:
                kingdom = "Eukaryota"
            elif filename.find("viruses") != -1:
                kingdom = "Viruses"
            elif filename.find("unclassified") != -1:
                kingdom = "Unclassified"
            else:
                logger.error(
                    f"Could not retrieve taxonomy kingdom from {html_file_path}.\n"
                    "Not scraping this file and adding proteins to the local databsae."
                )
                continue
                report, session = row_to_protein(
                    row,
                    family_name,
                    taxonomy_filters,
                    kingdom,
                    ec_filters,
                    session,
                    args,
                )
                if report["sql"] is not None:
                    sql_failures_logger.warning(report["sql"])

    if type(session) is dict:
        if args.output is not sys.stdout:
            output_path = args.output / f"cazy_dict_{time_stamp}.json"
            json.dump(session, output_path)
        else:
            json.dump(session, sys.stdout)

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    logger.info(
        "Finished scraping CAZy. Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}"
    )
    print(
        "=====================cazy_webscraper=====================\n"
        "Finished scraping CAZy\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}\n"
    )
    sys.exit(1)


def get_html_files(args):
    """Get paths to the HTML files.

    :param args: cmd-line args parser

    Return list of paths to HTML files.
    """
    logger = logging.getLogger(__name__)
    html_files = []

    # retrieve all files in dir
    files_in_dir = (entry for entry in Path(args.scrape_files).iterdir() if entry.is_file())

    for item in files_in_dir:
        if item.name.endswith(".html"):
            html_files.append(item)

    if len(html_files) == 0:
        logger.error(f"Not HTML files found in {args.scrape_files}.\nTerminating program")
        sys.exit(1)

    logger.warning(f"Retrieved {len(html_files)} from {args.scrape_files}")

    return html_files
