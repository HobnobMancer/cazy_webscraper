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
import time

from datetime import datetime
from typing import List, Optional

import pandas as pd

from Bio import Entrez
from bioservices import UniProt
from saintBioutils.misc import get_chunks_list
from saintBioutils.uniprot import get_uniprot_accessions
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from cazy_webscraper import closing_message, connect_existing_db
from cazy_webscraper.ncbi.gene_names import get_linked_ncbi_accessions
from cazy_webscraper.expand.uniprot.uniprot_cache import get_uniprot_cache, cache_uniprot_data
from cazy_webscraper.sql import sql_interface
from cazy_webscraper.sql.sql_interface.add_data.add_uniprot_data import (
    add_ec_numbers,
    add_pdb_accessions,
    add_uniprot_accessions,
)
from cazy_webscraper.sql.sql_interface.delete_data import delete_old_relationships, delete_old_annotations
from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks, get_table_dicts
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

    Entrez.email = args.email
    
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

    # get cached data, or return an empty dict
    # and return a list of GenBank/NCBI protein version accessions to query UniProt with
    # to download protein data
    uniprot_dict, gbk_data_to_download = get_uniprot_cache(gbk_dict, args)
    # uniprot_data[ncbi_acc] = {
    #     'uniprot_acc': uniprot_acc,
    #     'uniprot_entry_id': uniprot_entry_id,
    #     'protein_name': protein_name,
    #     'ec_numbers': ec_numbers,
    #     'sequence': sequence,
    #     'pdbs': all_pdbs,
    # }

    # gbk_data_to_download = list of GenBank accs to download data for

    if args.skip_download is False:
        logger.warning(f"Retrieving data for {len(gbk_data_to_download)} proteins")

        downloaded_uniprot_data, all_ecs = get_uniprot_data(gbk_data_to_download, cache_dir, args)

        uniprot_dict.update(downloaded_uniprot_data)

    else:
        logger.warning(f"Using only data from cache:\n{args.use_uniprot_cache}")

    if len(list(downloaded_uniprot_data.keys())) != 0:
        cache_uniprot_data(uniprot_dict, cache_dir, time_stamp)

    # add uniprot accessions (and sequences if seq retrieval is enabled)
    logger.warning("Adding data to the local CAZyme database")
    add_uniprot_accessions(uniprot_dict, connection, args)

    # add ec numbers
    if args.ec:
        logger.warning("Adding EC numbers to the local CAZyme database")
        add_ec_numbers(uniprot_dict, connection, args)
        logger.warning("Adding Genbanks-ECnumber relationships to local CAZyme db")
        add_genbank_ec_relationships(uniprot_dict, gbk_dict, connection, args)

        if args.delete_old_ec_relationships:
            logger.warning(
                "Deleting Genbanks-EC number annotations in the local CAZyme database\n"
                "that were not included for the protein whose additional data was just\n"
                "downloaded from UniProt"
            )
            # load ec numbers and relationships with Genbanks records from the local db
            ec_table_dict = get_table_dicts.get_ec_table_dict(connection)
            ec_gbk_table_dict = get_table_dicts.get_ec_gbk_table_dict(connection)

            delete_old_relationships(
                uniprot_dict,
                gbk_dict,
                ec_table_dict,
                ec_gbk_table_dict,
                'ec_numbers',
                'Genbanks_Ecs',
                connection,
                args,
            )

        if args.delete_old_ecs:
            logger.warning(
                "Deleting EC numbers in local db that are not linked to any Genbanks table records"
            )
            # load ec numbers and relationships with Genbanks records from the local db
            ec_table_dict = get_table_dicts.get_ec_table_dict(connection)
            ec_gbk_table_dict = get_table_dicts.get_ec_gbk_table_dict(connection)

            delete_old_annotations(ec_table_dict, ec_gbk_table_dict, 'Ecs', connection, args)

    # add pdb accessions
    if args.pdb:
        logger.warning("Adding RSCB PDB IDs to the local CAZyme database")
        add_pdb_accessions(uniprot_dict, gbk_dict, connection, args)
        add_pdb_gbk_relationships(uniprot_dict, gbk_dict, connection, args)

        if args.delete_old_pdb_relationships:
            logger.warning(
                "Deleting Genbanks-PDB annotations in the local CAZyme database\n"
                "that were not included for the protein whose additional data was just\n"
                "downloaded from UniProt"
            )
            # load ec numbers and relationships with Genbanks records from the local db
            pdb_table_dict = get_table_dicts.get_pdb_table_dict(connection)
            gbk_pdb_rel_table_dict = get_table_dicts.get_gbk_pdb_table_dict(connection)

            delete_old_relationships(
                uniprot_dict,
                gbk_dict,
                pdb_table_dict,
                gbk_pdb_rel_table_dict,
                'pdbs',
                'Genbanks_Pdbs',
                connection,
                args,
            )

        if args.delete_old_pdbs:
            logger.warning(
                "Deleting PDB accessions in local db that are not linked to any Genbanks table records"
            )
            # load ec numbers and relationships with Genbanks records from the local db
            pdb_table_dict = get_table_dicts.get_pdb_table_dict(connection)
            gbk_pdb_rel_table_dict = get_table_dicts.get_gbk_pdb_table_dict(connection)

            delete_old_annotations(pdb_table_dict, gbk_pdb_rel_table_dict, 'Pdbs', connection, args)

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
    logger = logging.getLogger(__name__)

    # retrieve dict of genbank accession and genbank db ids from the local CAZyme db
    if args.genbank_accessions is not None:
        logger.warning(f"Getting GenBank accessions from file: {args.genbank_accessions}")

        with open(args.genbank_accessions, "r") as fh:
            lines = fh.read().splitlines()
        
        accessions = [line.strip() for line in lines if len(line) != 0]
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


def get_uniprot_data(ncbi_accessions, cache_dir, args):
    """Batch query UniProt to retrieve protein data. Save data to cache directory.
    
    Bioservices requests batch queries no larger than 200.

    Note that according to Uniprot (June 2022), there are various limits on ID Mapping Job Submission:

    ========= =====================================================================================
    Limit	  Details
    ========= =====================================================================================
    100,000	  Total number of ids allowed in comma separated param ids in /idmapping/run api
    500,000	  Total number of "mapped to" ids allowed
    100,000	  Total number of "mapped to" ids allowed to be enriched by UniProt data
    10,000	  Total number of "mapped to" ids allowed with filtering
    ========= =====================================================================================

    :param ncbi_accessions: list of NCBI protein accessions to query UniProt with
    :param cache_dir: path to directory to write out cache
    :param args: cmd-line args parser
    
    Return
    Dict of data retrieved from UniProt and to be added to the db 
        uniprot_data[ncbi_acc] = {
            'uniprot_acc': uniprot_acc,
            'uniprot_entry_id': uniprot_entry_id,
            'protein_name': protein_name,
            'ec_numbers': ec_numbers,
            'sequence': sequence,
            'pdbs': all_pdbs,
        }
    """
    logger = logging.getLogger(__name__)

    failed_ids_cache = cache_dir / "ncbi_acc_not_in_uniprot"
    failed_connections_cache = cache_dir / "failed_connections_ncbi_acc"

    uniprot_dict = {}  # see doc string
    all_batches = {}  # {batch: int(tries)}

    bioservices_queries = get_chunks_list(
        ncbi_accessions,
        args.uniprot_batch_size,
    )

    for batch in bioservices_queries:
        all_batches[",".join(batch)] = 0  # num of attempts at connecting to UniProt

    while len(list(all_batches.keys())) > 0:
        batches_to_process = copy(all_batches)

        for batch in tqdm(batches_to_process, "Batch retrieving protein data from UniProt"):
            # batch is a string of comma separated NCBI protein version accessions
            mappings = map_to_uniprot(batch)

            if mappings is None:  # could not connect to UniProt
                all_batches[batch] += 1

                if all_batches[batch] > args.retries:  # run out of attempts to retry the connection
                    del all_batches[batch]
                    logger.warning(
                        f"Failed to retrieve data for batch after {args.retries} attemps\n"
                        "Out of retries.\n"
                        "NCBI accessions will be written to cache"
                    )
                    failed_acc = batch.split(",")
                    with open("failed_connections_cache", a) as fh:
                        for acc in failed_acc:
                            fh.write(f"{acc}\n")
                        
                else:  # still attempts remaining
                    logger.warning(
                        f"Failled to retrieve data from UniProt for batch after {all_batches[batch]}"
                        f"/{args.retries} attemtps\n"
                        "Will retry later."
                    )

                continue

            try:  # some NCBI accessions could not be mapped to a record in UniProt
                with open(failed_ids_cache, "a") as fh:
                    for not_catalogued_acc in mappings['failedIds']:
                        fh.write(f"{not_catalogued_acc}\n")
            except KeyError:  # may not be any failed Ids
                pass
            
            try:  # Mapped UniProt records
                mapping_results = mappings['results']  # used mostly for try/except

                for mapping_result in mapping_results:
                    ncbi_acc = mapping_result['from']
                    mapped_record = mapping_result['to']

                    (
                        uniprot_acc,
                        uniprot_entry_id,
                        protein_name,
                        ec_numbers,
                        sequence,
                        all_pdbs,
                        matching_record,
                    ) = extract_protein_data(mapped_record, ncbi_acc)

                    if matching_record is False:  # bool, if mapped record contains ncbi accession
                        continue
                    
                    uniprot_data[ncbi_acc] = {
                        'uniprot_acc': uniprot_acc,
                        'uniprot_entry_id': uniprot_entry_id,
                        'protein_name': protein_name,
                        'ec_numbers': ec_numbers,
                        'sequence': sequence,
                        'pdbs': all_pdbs,
                    }

            except KeyError: # may not be any ncbi acc that mapped to a UniProt record
                pass

    return uniprot_data


def mapping_decorator(func):
    """Decorator to retry the wrapped function up, up to 'retries' times"""

    def wrapper(*args, retries=10, **kwards):
        tries, success, response = 0, False, None
        
        while not success and (tries < retries):
            response = func(*args, **kwards)
            
            if response is not None:
                success = True
            else:
                print(
                    f"Could not connect to UniProt on attempt no.{tries} of {retries} tries.\n"
                    "Retrying in 10 seconds"
                )
                tries += 1
                time.sleep(10)

        return response
    
    return wrapper


@mapping_decorator
def map_to_uniprot(accessions):
    """Map accessions to records in UniProt
    
    :param accessions: list, list of NCBI protein version accessions
    
    Return dict. Successful mappings under 'results', and failed mappings under 'failedIds' 
    Failed mappings are NCBI accessions that are not listed in UniProt
    
    Or returns None if connection could not be made
    """
    mapping_results = UniProt().mapping(
        fr="EMBL-GenBank-DDBJ_CDS",
        to="UniProtKB",
        query=",".join(accessions),  # str of ids, separated by commas
    )
    
    return mapping_results


if __name__ == "__main__":
    main()
