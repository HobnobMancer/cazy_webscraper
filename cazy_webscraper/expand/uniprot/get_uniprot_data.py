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


import logging
import urllib.parse
import urllib.request

import pandas as pd

from datetime import datetime
from typing import List, Optional
from urllib.error import HTTPError

from bioservices import UniProt
from saintBioutils.uniprot import get_uniprot_accessions
from saintBioutils.utilities.file_io import make_output_directory
from tqdm import tqdm

from cazy_webscraper import cazy_scraper, closing_message
from cazy_webscraper.expand import get_chunks_list
from cazy_webscraper.sql import sql_interface
from cazy_webscraper.sql.sql_interface import get_selected_gbks
from cazy_webscraper.sql.sql_interface.add_uniprot_data import (
    add_ec_numbers,
    add_pdb_accessions,
    add_uniprot_accessions,
)
from cazy_webscraper.sql import sql_orm
from cazy_webscraper.utilities import config_logger
from cazy_webscraper.utilities.parsers import uniprot_parser
from cazy_webscraper.utilities.parse_configuration import get_expansion_configuration


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
    connection, logger_name, cache_dir = cazy_scraper.connect_existing_db(args, time_stamp, start_time)

    # build cache directory
    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "uniprot_data_retrieval"
        make_output_directory(cache_dir, args.force, args.nodelete_cache)

    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = get_expansion_configuration(args)

    # add log to the local CAZyme database
    logger.info("Adding log of scrape to the local CAZyme database")

    retrieved_annotations = "UniProt accessions"
    if args.ec:
        retrieved_annotations += ", EC numbers"
    if args.pdb:
        retrieved_annotations = ", PDB accessions"
    if args.sequence:
        retrieved_annotations = ", Protein sequence"
    if args.seq_update:
        retrieved_annotations = ", Updated UniProt protein sequences"

    with sql_orm.Session(bind=connection) as session:
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            'UniProt',
            retrieved_annotations,
            session,
            args,
        )

    # retrieve dict of genbank accession and genbank db ids from the local CAZyme db
    gbk_dict = get_selected_gbks.get_genbank_accessions(
        class_filters,
        family_filters,
        taxonomy_filter_dict,
        kingdom_filters,
        ec_filters,
        connection,
    )

    logger.warning(f"Retrieving UniProt data for {len(gbk_dict.keys())}")

    # retrieve the uniprot accessions for the genbank accessions
    uniprot_gkb_dict = get_uniprot_accessions(gbk_dict, args)  # {uniprot_acc: {'gbk_acc': str, 'db_id': int}}

    uniprot_dict, all_ecs = get_uniprot_data(uniprot_gkb_dict, cache_dir)

    # add uniprot accessions (and sequences if seq retrieval is enabled)
    add_uniprot_accessions(uniprot_dict, gbk_dict, connection, args)

    # add ec numbers
    if (args.ec) and (len(all_ecs) != 0):
        add_ec_numbers(uniprot_dict, all_ecs, gbk_dict, connection, args)

    # add pdb accessions
    if args.pdb:
        add_pdb_accessions(uniprot_dict, gbk_dict, connection, args)

    closing_message("get_uniprot_data", start_time, args)


def get_uniprot_data(uniprot_gbk_dict, cache_dir, args):
    """Batch query UniProt to retrieve protein data. Save data to cache directory.
    
    Bioservices requests batch queries no larger than 200.

    :param uniprot_gbk_dict: dict, keyed by UniProt accession and valued by dict of a GenBank accession
        and the local CAZyme db record ID {uniprot_acc: {'gbk_acc': str, 'db_id': int}}
    :param cache_dir: path to directory to write out cache
    :param args: cmd-line args parser
    
    Return
    Dict of data retrieved from UniProt and to be added to the db 
        {uniprot_acccession: {gbk_acccession: str, uniprot_name: str, pdb: set, ec: set}}
    Set of all retrieved EC numbers
    """
    logger = logging.getLogger(__name__)

    uniprot_dict = {}  # {uniprot: {genbank_accession: str, uniprot_name: str, pdb: set, ec: set}}
    all_ecs = set()  # store all retrieved EC numbers as one tuple per unique EC number

    # add PDB column to columns to be retrieved
    UniProt()._valid_columns.append('database(PDB)')
    
    bioservices_queries = get_chunks_list(
        list(uniprot_gbk_dict.keys()),
        args.bioservices_batch_size,
    )

    for query in tqdm(bioservices_queries, "Batch retrieving protein data from UniProt"):
        uniprot_df = UniProt().get_df(entries=query, limit=None)

        # cache UniProt response
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
            uniprot_acc = row['Entry'].strip()
            uniprot_name = row['Protein names'].strip()

            # remove quotation marks from the protein name, else an SQL error will be raised on insert
            uniprot_name = uniprot_name.replace("'", "")
            uniprot_name = uniprot_name.replace('"', '')
            uniprot_name = uniprot_name.replace("`", "")

            # checked if parsed before incase bioservices returned duplicate proteins
            try:
                uniprot_gbk_dict[uniprot_acc]
            except KeyError:
                logger.warning(
                    f"UniProt ID {uniprot_acc} was retrieved from UniProt\n"
                    "But no corresponding record was found in the local CAZyme database\n"
                    "Skipping this UniProt ID"
                )

            try:
                uniprot_dict[uniprot_acc]
                logger.warning(
                    f'Multiple entries for UniProt:{uniprot_acc}, '
                    f'GenBank:{uniprot_gbk_dict[uniprot_acc]} retrieved from UniProt,\n'
                    'compiling data into a single record'
                )
            except KeyError:
                try:
                    uniprot_dict[uniprot_acc] = {
                        "genbank_accession": uniprot_gbk_dict[uniprot_acc],
                        "name": uniprot_name,
                    }
                except KeyError:
                    logger.warning(
                        f"Retrieved record with UniProt accession {uniprot_acc} but this "
                        "accession was not\nretrieved from the UniProt REST API"
                    )
                    continue
            
            if args.ec:
                # retrieve EC numbers
                ec_numbers = row['EC number']
                try:
                    ec_numbers = ec_numbers.split('; ')
                except AttributeError:
                    # no EC numbers listed
                    ec_numbers = set()
                
                try:
                    uniprot_dict[uniprot_acc]["ec"]
                except KeyError:
                    uniprot_dict[uniprot_acc]["ec"] = set()

                # add EC numbers to dict
                for ec in ec_numbers:
                    all_ecs.add( (ec,) )
                    uniprot_dict[uniprot_acc]["ec"].add(ec.strip())

            if args.pdb:
                # retrieve PDB accessions
                pdb_accessions = row['Cross-reference (PDB)']
                print('PDBS:', pdb_accessions)
                try:
                    pdb_accessions = pdb_accessions.split('; ')
                except AttributeError:
                    pdb_accessions = set()

                try:
                    uniprot_dict[uniprot_acc]["pdb"]
                except KeyError:
                    uniprot_dict[uniprot_acc]["pdb"] = set()

                # add PDB accessions to dict
                for pdb in pdb_accessions:
                    uniprot_dict[uniprot_acc]["pdb"].add(pdb.strip())
            
            if args.sequence:
                sequence = row['Sequence']

                try:
                    uniprot_dict[uniprot_acc]["sequence"]
                    existing_date = uniprot_dict[uniprot_acc]["seq_date"]
                    new_date = row['Date of last sequence modification']
                    
                    # check which sequence is newer
                    existing_date.split('-')
                    existing_date = datetime(existing_date[0], existing_date[1], existing_date[2])
                    new_date.split('-')
                    new_date = datetime(existing_date[0], existing_date[1], existing_date[2])

                    if new_date > existing_date:  # past < present is True
                        uniprot_dict[uniprot_acc]["sequence"] = sequence
                        uniprot_dict[uniprot_acc]["seq_date"] = row['Date of last sequence modification']
                    # else keep the existing sequence
                    logger.warning(
                        f'Multiple sequences retrieved for {uniprot_acc}\n'
                        'Using most recently updated sequence'
                    )
                    
                except KeyError:
                    uniprot_dict[uniprot_acc]["sequence"] = sequence
                    uniprot_dict[uniprot_acc]["seq_date"] = row['Date of last sequence modification']

    return uniprot_dict, all_ecs


if __name__ == "__main__":
    main()
