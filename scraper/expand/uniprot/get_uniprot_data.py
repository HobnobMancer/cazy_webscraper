#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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
"""Retrieve data from UniProtKB and adding it to a local CAZyme db"""


import io
import logging
import urllib.parse
import urllib.request

import pandas as pd

from datetime import datetime
from typing import List, Optional
from urllib.error import HTTPError

from bioservices import UniProt
from tqdm import tqdm

from scraper import cazy_webscraper
from scraper.expand import get_chunks
from scraper.sql import sql_interface
from scraper.sql.sql_interface.add_uniprot_db import (
    add_ec_numbers,
    add_pdb_accessions,
    add_uniprot_accessions,
)
from scraper.sql import sql_orm
from scraper.utilities import config_logger, file_io
from scraper.utilities.parsers import uniprot_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # Program preparation
    if argv is None:
        parser = uniprot_parser.build_parser()
        args = parser.parse_args()
    else:
        parser = uniprot_parser.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)
    
    # parse the configuration data (cache the uniprot data as .csv files)
    connection, logger_name, cache_dir = cazy_webscraper.connect_existing_db(args, time_stamp)

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "uniprot_data_retrieval"
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)

    config_dict, kingdom_filters, taxonomy_filter_dict = parse_uniprot_config(args)

    # add log to the local CAZyme database
    logger.info("Adding log of scrape to the local CAZyme database")
    with sql_orm.Session(bind=connection) as session:
        retrieved_annotations = "UniProt accessions"
        if len(config_dict['ec']) != 0:
            retrieved_annotations = f"{retrieved_annotations}, EC numbers"
        if len(config_dict['pdb']) != 0:
            retrieved_annotations = f"{retrieved_annotations}, PDB accessions"
        if len(config_dict['seq']) != 0:
            retrieved_annotations = f"{retrieved_annotations}, Protein sequence"
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            kingdom_filters,
            taxonomy_filter_dict,
            'UniProt',
            retrieved_annotations,
            session,
            args,
        )

    # retrieve dict of genbank accession and genbank accession ids from the local CAZyme db
    gbk_dict = sql_interface.get_gbk_table_dict(connection)
    genbank_accessions = list(gbk_dict.keys())

    # retrieve the uniprot accessions for the genbank accessions
    uniprot_gkb_dict = get_uniprot_accessions(genbank_accessions)  # {uniprot_acc: gbk_acc}
    uniprot_dict, all_ecs = get_uniprot_data(uniprot_gkb_dict, cache_dir)

    # add uniprot accessions (and sequences if seq retrieval is enabled)
    add_uniprot_accessions(uniprot_dict, gbk_dict, connection)

    # add ec numbers
    if args.ec:
        add_ec_numbers(uniprot_dict, all_ecs, gbk_dict, connection)

    # add pdb accessions
    if args.pdb:
        add_pdb_accessions(uniprot_dict, gbk_dict, connection)

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    if args.verbose:
        logger.info(
            "Finished scraping CAZy. Terminating program.\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_webscraper.V}\n"
            f"Citation: {cazy_webscraper.CITATION_INFO}"
        )
    else:
        print(
            "=====================cazy_webscraper=====================\n"
            "Finished scraping CAZy\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_webscraper.VERSION_INFO}\n"
            f"Citation: {cazy_webscraper.CITATION_INFO}"
        )


def get_uniprot_accessions(genbank_accessions, args):
    """Retrieve UniProt accessions for the GenBank accessions from UniProtKB.
    
    UniProt requests batch queries of no larger than 20,000

    :param genbank_accessions: set, GenBank accessions to retrieve UniProt accessions for
    :param args: cmd-line args parser

    Return dict of {uniprot_accession: genbank_accession}
    """
    uniprot_url = 'https://www.uniprot.org/uploadlists/'

    gbk_accs_chunks = get_chunks(genbank_accessions, args.uniprot_batch_size)

    accession_dict = {}  # {genbank_accession: uniprot_accession}

    for query_chunk in tqdm(gbk_accs_chunks, desc='Batch retrieving UniProt accessions'):
        # convert the set of gbk accessions into str format
        query = ' '.join(query_chunk)

        params = {
            'from': 'EMBL',
            'to': 'ACC',
            'format': 'tab',
            'query': query
        }

        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(uniprot_url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        uniprot_batch_response = response.decode('utf-8')

        uniprot_batch_response = uniprot_batch_response.split('\n')

        for line in uniprot_batch_response[1:-1]:  # the first line includes the titles, last line is an empty str
            accession_dict[line.split('\t')[1]] = line.split('\t')[0]

    return accession_dict


def get_uniprot_data(uniprot_gkb_dict, cache_dir, args):
    """Batch query UniProt to retrieve protein data. Save data to cache directory.
    
    Bioservices requests batch queries no larger than 200.

    :param gkb_uniprot_dict: dict, keyed by GenBank accession and valued by UniProt accession
    :param cache_dir: path to directory to write out cache
    :param args: cmd-line args parser
    
    Return
    Dict of data retrieved from UniProt and to be added to the db 
        {uniprot_acccession: {gbk_acccession: str, uniprot_name: str, pdb: set, ec: set}}
    Set of all retrieved EC numbers
    """
    logger = logging.getLogger(__name__)

    uniprot_dict = {}  # {uniprot: {genbank_accession: str, uniprot_name: str, pdb: set, ec: set}}

    # add PDB column to columns to be retrieved
    UniProt()._valid_columns.append('database(PDB)')
    
    query_chunks = get_chunks(list(uniprot_gkb_dict.keys()), args.bioservices_batch_size)

    all_ecs = set()  # store all retrieved EC numbers as one tuple per unique EC number

    for query_list in tqdm(query_chunks, "Batch quering UniProt"):
        uniprot_df = UniProt().get_df(entries=query_list)

        _time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        _path = cache_dir / f'uniprot_query_response_{_time_stamp}.csv'
        uniprot_df.to_csv(_path)

        index = 0
        uniprot_df = uniprot_df[[
            'Entry',
            'Protein names',
            'EC number',
            'Sequence',
            'Cross-reference (PDB)',
        ]]

        for index in tqdm(range(len(uniprot_df['Entry'])), desc="Parsing UniProt response"):
            row = uniprot_df.iloc[index]
            uniprot_acc = row['Entry']
            uniprot_name = row['Protein names']

            # checked if parsed before incase bioservices returned duplicate proteins
            try:
                uniprot_dict[uniprot_acc]
                logger.warning(
                    f'Multiple entries for UniProt:{uniprot_acc}, '
                    f'GenBank:{uniprot_gkb_dict[uniprot_acc]} retrieved from UniProt,\n'
                    'compiling data into a single record'
                )
            except KeyError:
               uniprot_dict[uniprot_acc] = {
                   "genbank_accession": uniprot_gkb_dict[uniprot_acc],
                   "name": uniprot_name,
                }
            
            if args.ec:
                ec_numbers = row['EC number']
                try:
                    ec_numbers = ec_numbers.split('; ')
                except AttributeError:
                    # no EC numbers listed
                    continue

                try:
                    uniprot_dict[uniprot_acc]["ec"]
                except KeyError:
                    uniprot_dict[uniprot_acc]["ec"] = set()

                for ec in ec_numbers:
                    all_ecs.add( (ec,) )
                    uniprot_dict[uniprot_acc]["ec"].add(ec)

            if args.pdb:
                pdb_accessions = row['Cross-reference (PDB)']
                try:
                    pdb_accessions = pdb_accessions.split('; ')
                except AttributeError:
                    # no EC numbers listed
                    continue

                try:
                    uniprot_dict[uniprot_acc]["pdb"]
                except KeyError:
                    uniprot_dict[uniprot_acc]["pdb"] = set()

                for pdb in pdb_accessions:
                    uniprot_dict[uniprot_acc]["pdb"].add(pdb)
            
            if args.sequences:
                sequence = row['Sequence']

                try:
                    uniprot_dict[uniprot_acc]["sequence"]
                    logger.warning(
                        f'Multiple sequences retrieved for {uniprot_acc}\n'
                        'Adding first retrieved sequence to the local CAZyme db'
                    )
                except KeyError:
                    uniprot_dict[uniprot_acc]["sequence"] = sequence

    return uniprot_dict, all_ecs


if __name__ == "__main__":
    main()
