#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Web scraper to scrape the CAZy website."""


import logging

import pandas as pd

from datetime import datetime


__version__ = "2.0.0-beta"

VERSION_INFO = f"cazy_webscraper version: {__version__}"

CITATION_INFO = (
    "If you use cazy_webscraper in your work, please cite the following publication:\n"
    "\tHobbs, E. E. M., Pritchard, L., Chapman, S., Gloster, T. M.,\n"
    "\t(2021) cazy_webscraper Microbiology Society Annual Conference 2021 poster.\n"
    "\tFigShare. Poster.\n"
    "\thttps://doi.org/10.6084/m9.figshare.14370860.v7"
)


def closing_message(job, start_time, args):
    """Write closing messsage to terminal"""
    logger = logging.getLogger(__name__)

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    if args.verbose:
        logger.info(
            "Finished scraping CAZy. Terminating program.\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {VERSION_INFO}\n"
            f"Citation: {CITATION_INFO}"
        )
    else:
        print(
            f"====================={job}=====================\n"
            "Finished scraping CAZy\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {VERSION_INFO}\n"
            f"Citation: {CITATION_INFO}"
        )

    return
