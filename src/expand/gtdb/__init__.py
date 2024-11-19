#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) Jame Hutton Institute 2022
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
#
# Bio.PDB reference:
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
# implemented in Python. Bioinformatics 19: 2308–2310
"""Retrieve taxonomic classifications from GTDB"""


import logging
import sys
import time

from requests.exceptions import ConnectionError, MissingSchema
from socket import timeout
from urllib3.exceptions import HTTPError, RequestError
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import mechanicalsoup

from tqdm import tqdm


GTDB_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"


def get_gtdb_data(args, cache_dir, arch, bact):
    """Download taxonomic data files from the GTDB website

    :param args: cmd-line args parser
    :param cache_dir: path to cache directory
    :param arch: bool, whether to download arch datafile
    :param bact: bool, whether to download the bact datafile

    Return paths to downloaded archea and bacteria data files. None if kingdom not selected for DL
    """
    logger = logging.getLogger(__name__)

    gtdb_release_page, error_message = get_page(
        GTDB_URL,
        max_tries=args.retries,
    )

    if gtdb_release_page is None:
        logger.error(f"Failed to get GTDB data release page:\n{error_message}\nTerminating program")
        sys.exit(1)

    archaea_link, bacteria_link = None, None

    if args.archaea_file is None or args.bacteria_file is None:
        # need to retrieve links from the website
        for i in gtdb_release_page.select("table")[0].select("tr"):
            for j in i.select("td"):
                if j.contents[0] is not None:
                    try:
                        if j.contents[0]['href'].endswith('_taxonomy.tsv.gz'):
                            if j.contents[0]['href'].split("/")[-1].startswith('ar'):
                                archaea_link = f"{GTDB_URL}{j.contents[0]['href']}"
                            else:
                                bacteria_link = f"{GTDB_URL}{j.contents[0]['href']}"
                    except (KeyError, TypeError):
                        continue

        if archaea_link is None or bacteria_link is None:
            if archaea_link is None and 'archaea' in args.taxs:
                logger.error(
                    "Failed to get archeae GTDB data release page\n"
                    "Retrieved datafile urls:\n"
                    f"Archaea: {archaea_link}"
                )
            if bacteria_link is None and 'bacteria' in args.taxs:
                logger.error(
                    "Failed to get bacteria GTDB data release page\n"
                    "Retrieved datafile urls:\n"
                    f"Archaea: {bacteria_link}"
                )
            logger.error("Failed to retrieve download URLs\nTerminating program")
            sys.exit(1)

    archaea_file, bacteria_file = None, None

    if 'archaea' in args.taxs:
        if args.archaea_file is not None:
            archaea_file = args.archaea_file
        else:
            arch_release = archaea_link.split("/")[-1].split("_")[0]
            archaea_file = cache_dir / f"archaea_data-{arch_release}.gz"

            downloaded = download_gtdb_data(archaea_link, archaea_file, 'Archaea')

            if downloaded is False:
                logger.error("Failed to download archaea GTDB data file\nTerminating program")
                sys.exit(1)

    if 'bacteria' in args.taxs:
        if args.bacteria_file is not None:
            bacteria_file = args.bacteria_file
        else:
            bact_release = bacteria_link.split("/")[-1].split("_")[0]
            bacteria_file = cache_dir / f"bacteria_data-{bact_release}.gz"

            downloaded = download_gtdb_data(bacteria_link, bacteria_file, 'Bacteria')

            if downloaded is False:
                logger.error("Failed to download bacteria GTDB data file\nTerminating program")
                sys.exit(1)

    return archaea_file, bacteria_file


def download_gtdb_data(url, out_file_path, gtdb_group):
    logger = logging.getLogger(__name__)

    # Try URL connection
    try:
        response = urlopen(url, timeout=45)
    except (HTTPError, URLError, timeout) as e:
        logger.error(
            f"Failed to download {gtdb_group} gtdb data file from: {url}", exc_info=1,
        )
        return False

    logger.info(f"Downloading data file from:\n{url}")

    file_size = int(response.info().get("Content-length"))
    bsize = 1_048_576

    try:
        with open(out_file_path, "wb") as out_handle:
            with tqdm(
                total=file_size,
                leave=False,
                desc=f"Downloading gtdb {gtdb_group} data file",
            ) as pbar:
                while True:
                    buffer = response.read(bsize)
                    if not buffer:
                        break
                    pbar.update(len(buffer))
                    out_handle.write(buffer)
    except IOError:
        logger.error(f"Download failed GTDB {gtdb_group} data file", exc_info=1)
        return False

    return True


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
def get_page(url, **kwargs):
    """Create browser and use browser to retrieve page for given URL.

    :param url: str, url to webpage
    :param args: cmd-line args parser
    :kwargs max_tries: max number of times connection to CAZy can be attempted

    Return browser response object (the page).
    """
    # create browser object
    browser = mechanicalsoup.Browser()
    # create response object
    page = browser.get(url, timeout=10)
    page = page.soup

    return page
