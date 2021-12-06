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
"""Explore the number of GenBank genomes annotated by CAZy."""


import logging

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio import Entrez
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from cazy_webscraper import cazy_scraper


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # Program preparation
    if argv is None:
        parser = genbank_parser.build_parser()
        args = parser.parse_args()
    else:
        parser = genbank_parser.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)

    Entrez.email = args.email

    # connect to the local CAZyme database
    connection, logger_name, cache_dir = cazy_scraper.connect_existing_db(args, time_stamp)

    # check if need to build output dir
    ???



    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    if args.verbose:
        logger.info(
            "Finished getting data from UniProt\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_scraper.VERSION_INFO}\n"
            f"Citation: {cazy_scraper.CITATION_INFO}"
        )
    else:
        print(
            "=====================cazy_webscraper=====================\n"
            "Finished getting data from UniProt\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_scraper.VERSION_INFO}\n"
            f"Citation: {cazy_scraper.CITATION_INFO}"
        )


def get_genomic_accessions():
    """Retrieve genomic accessions for the protein accessions in the local db.
    
    :param connection: open sqlalchemy connection to a SQLite db
    :param args: cmd-line args parser
    
    Return dict {kingdom: {genomic_acc: {'proteins': [p_acc], 'count':#ofProteins}}}
    """
    # load the Genbank and Kingdom records from the local CAZyme db
    genbank_kingdom_records = []

    # extract the GenBank accessions and Kingdoms from the records
    genbank_kingdom_dict = {}  # kingdom: {genbank_accessions}
    for record in genbank_kingdom_records:
        kingdom = record[1].kingdom
        gbk_acc = record[0].genbank_accession

        try:
            genbank_kingdom_dict[kingdom].add(gbk_acc)
        except KeyError:
            genbank_kingdom_dict[kingdom] = gbk_acc

    genomic_accession_dict = {}  #  {kingdom: {genomic_acc: {'proteins': [p_acc], 'count':#ofProteins}}}
    
    for kingdom in tqdm(genbank_kingdom_dict, desc="Retrieving genomic accessions per kingdom"):
        gbk_accessions = genbank_kingdom_dict[kingdom]

        # break up the list into a series of smaller lists that can be batched querried
        batch_queries = []

        for batch_query in batch_queries:

            # 

if __name__ == "__main__":
    main()
