#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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


import json
import logging

import pandas as pd

from datetime import datetime
from typing import List, Optional

from bioservices import UniProt
from saintBioutils.misc import get_chunks_list
from saintBioutils.uniprot import get_uniprot_accessions
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from cazy_webscraper import closing_message, connect_existing_db
from cazy_webscraper.sql import sql_interface
from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks
from cazy_webscraper.sql.sql_interface.add_data.add_uniprot_data import (
    add_ec_numbers,
    add_pdb_accessions,
    add_uniprot_accessions,
)
from cazy_webscraper.sql import sql_orm
from cazy_webscraper.utilities.parsers.uniprot_parser import build_parser
from cazy_webscraper.utilities.parse_configuration import get_expansion_configuration

import sys


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # Program preparation
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)
    
    # parse the configuration data (cache the uniprot data as .csv files)
    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)

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
    add_db_log(
        config_dict,
        taxonomy_filter_dict,
        kingdom_filters,
        ec_filters,
        time_stamp,
        connection,
        args,
    )

    # get dick to GenBank accessions and local db IDs
    gbk_dict = get_db_gbk_accs(
        class_filters,
        family_filters,
        taxonomy_filter_dict,
        kingdom_filters,
        ec_filters,
        connection,
        args,
    )

    # get cached data, or return an empty dict and set
    uniprot_dict, all_ecs, gbk_data_to_download = get_uniprot_cache(gbk_dict, args)

    if args.skip_download is False:
        logger.warning(f"Retrieving data for {len(gbk_data_to_download)} proteins")

        downloaded_uniprot_data, all_ecs = get_uniprot_data(uniprot_dict, gbk_data_to_download, cache_dir, args)

        uniprot_dict.update(downloaded_uniprot_data)

    else:
        logger.warning(f"Using only data from cache:\n{args.use_uniprot_cache}")

    if len(list(downloaded_uniprot_data.keys())) != 0:
        cache_uniprot_data(uniprot_dict, cache_dir, time_stamp)

    # add uniprot accessions (and sequences if seq retrieval is enabled)
    logger.warning("Adding data to the local CAZyme database")
    add_uniprot_accessions(uniprot_dict, gbk_dict, connection, args)

    # add ec numbers
    if (args.ec) and (len(all_ecs) != 0):
        logger.warning("Adding EC numbers to the local CAZyme database")
        add_ec_numbers(uniprot_dict, all_ecs, gbk_dict, connection, args)

    # add pdb accessions
    if args.pdb:
        logger.warning("Adding RSCB PDB IDs to the local CAZyme database")
        add_pdb_accessions(uniprot_dict, gbk_dict, connection, args)

    closing_message("get_uniprot_data", start_time, args)


def add_db_log(
    config_dict,
    taxonomy_filter_dict,
    kingdom_filters,
    ec_filters,
    time_stamp,
    connection,
    args,
):
    """Add log to local db

    :param connection: open connection to SQLite db engine
    :param args: CLI args parser
    
    Return nothing
    """
    retrieved_annotations = "UniProt accessions, Protein names"
    if args.ec:
        retrieved_annotations += ", EC numbers"
    if args.pdb:
        retrieved_annotations += ", PDB accessions"
    if args.sequence:
        retrieved_annotations += ", Protein sequence"
    if args.seq_update:
        retrieved_annotations += ", Updated UniProt protein sequences"

    with sql_orm.Session(bind=connection) as session:
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            taxonomy_filter_dict,
            kingdom_filters,
            ec_filters,
            'UniProt',
            retrieved_annotations,
            session,
            args,
        )


def get_db_gbk_accs(
    class_filters,
    family_filters,
    taxonomy_filter_dict,
    kingdom_filters,
    ec_filters,
    connection,
    args
):
    """Get the GenBank accessions of proteins to retrieve UniProt data for.

    :param args: CLI args parser

    Return dict {gbk_acc: local db id}
    """
    # retrieve dict of genbank accession and genbank db ids from the local CAZyme db
    if args.genbank_accessions is not None:
        logger.warning(f"Getting GenBank accessions from file: {args.genbank_accessions}")

        with open(args.genbank_accessions, "r") as fh:
            lines = fh.read().splitlines()
        
        accessions = [line.strip() for line in lines]
        accessions = set(accessions)

        gbk_dict = get_selected_gbks.get_ids(accessions, connection)

    else:
        gbk_dict = get_selected_gbks.get_genbank_accessions(
            class_filters,
            family_filters,
            taxonomy_filter_dict,
            kingdom_filters,
            ec_filters,
            connection,
        )

    logger.warning(f"Retrieving UniProt data for {len(gbk_dict.keys())}")

    return gbk_dict


def get_uniprot_cache(gbk_dict, args):
    """Get cached UniProt data, or return empty dict and set.
    
    :param gbk_dict: {gbk acc: local db id}
    :param args: CLI args parser

    Return dict {uniprot: {genbank_accession: str, uniprot_name: str, pdb: set, ec: set}}
    and set of EC numbers
    and list of GenBank accessions to download UniProt data for
    """
    # if using cachce skip accession retrieval
    uniprot_dict = {}  # {uniprot: {genbank_accession: str, uniprot_name: str, pdb: set, ec: set}}
    all_ecs = set()
    gbk_data_to_download = []

    if args.use_uniprot_cache is not None:
        logger.warning(f"Getting UniProt data from cache: {args.use_uniprot_cache}")

        with open(args.use_uniprot_cache, "r") as fh:
            uniprot_dict = json.load(fh)

        if args.ec:
            all_ecs = get_ecs_from_cache(uniprot_dict)
    
    if args.skip_download:  # only use cached data
        return uniprot_dict, all_ecs, gbk_data_to_download

    # else: check for which GenBank accessions data still needs be retrieved from UniProt
    # if some of the data is used from a cache, if no data is provided from a cache
    # retrieve data for all GenBank accesisons matching the provided criteria
    if len(list(uniprot_dict.keys())) != 0: 
        for uniprot_acc in tqdm(uniprot_dict):
            gbk_data_to_download.append(uniprot_dict[uniprot_acc]['genbank_accession'])
    else:  # get data for all GenBank accessions from the local db matching the user criteria
        gbk_data_to_download = list(gbk_dict.keys())
    
    return uniprot_dict, all_ecs, gbk_data_to_download


def get_uniprot_data(uniprot_dict, gbk_data_to_download, cache_dir, args):
    """Batch query UniProt to retrieve protein data. Save data to cache directory.
    
    Bioservices requests batch queries no larger than 200.

    :param uniprot_dict: dict, keyed by UniProt accession and valued by dict of a GenBank accession
        and the local CAZyme db record ID {uniprot_acc: {'gbk_acc': str, 'db_id': int}}
    :param gbk_data_to_download: list of NCBI protein accessions to query UniProt with
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
            'Date of last sequence modification',
            'Cross-reference (PDB)',
        ]]


        for index in tqdm(range(len(uniprot_df['Entry'])), desc="Parsing UniProt response"):
            row = uniprot_df.iloc[index]
            uniprot_acc = row['Entry'].strip()
            try:
                uniprot_name = row['Protein names'].strip()

            except AttributeError:
                logger.warning(
                    f"Protein name {row['Protein names']} was returned as float not string. Converting to string"
                )
                uniprot_name = str(row['Protein names'])

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
                continue  # uniprot ID does not match any retrieved previously

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
                try:
                    pdb_accessions = pdb_accessions.split(';')
                    pdb_accessions = [pdb.strip() for pdb in pdb_accessions if len(pdb.strip()) > 0]
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


def get_ecs_from_cache(uniprot_dict):
    """Extract all unique EC numbers from the UniProt data cache.
    
    :param uniprot_dict: dict of data retrieved from UniProt.
    
    Return set of EC numbers.
    """
    all_ecs = set()

    for uniprot_acc in tqdm(uniprot_dict, desc="Getting EC numbers from cached data"):
        try:
            ecs = uniprot_dict[uniprot_acc]["ec"]
            for ec in ecs:
                all_ecs.add( (ec,) )
        except (ValueError, TypeError, KeyError):
            pass

    return all_ecs


def cache_uniprot_data(uniprot_dict, cache_dir, time_stamp):
    """Cache data retrieved from UniProt.

    :param

    Return nothing
    """
    # cache updated UniProt data
    for uniprot_accession in uniprot_dict:
        try:
            uniprot_dict[uniprot_accession]['ec'] = list(uniprot_dict[uniprot_accession]['ec'])
        except KeyError:
            pass
        try:
            uniprot_dict[uniprot_accession]['pdb'] = list(uniprot_dict[uniprot_accession]['pdb'])
        except KeyError:
            pass

    uniprot_acc_cache = cache_dir / f"uniprot_data_{time_stamp}.json"
    with open(uniprot_acc_cache, "w") as fh:
        json.dump(uniprot_dict, fh) 


if __name__ == "__main__":
    main()
