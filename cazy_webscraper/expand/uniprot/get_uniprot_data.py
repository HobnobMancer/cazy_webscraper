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

from Bio import Entrez
from bioservices import UniProt
from saintBioutils.misc import get_chunks_list
from saintBioutils.uniprot import get_uniprot_accessions
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from cazy_webscraper import closing_message, connect_existing_db
from cazy_webscraper.ncbi.gene_names import get_linked_ncbi_accessions
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

    # get cached data, or return an empty dict and set
    #####
    ####
    ####
    ### NEEDS UPDATING
    ####
    ####
    #####
    uniprot_dict, all_ecs, gbk_data_to_download = get_uniprot_cache(gbk_dict, args)

    # uniprot_dict = {uniprot_acc: {gene_name: str, protein_name: str, pdb: set, ec: set, sequence:str, seq_data:str}}
    # all_ecs = set of EC numers
    # gbk_data_to_download = list of GenBank accs to download data for

    if args.skip_download is False:
        logger.warning(f"Retrieving data for {len(gbk_data_to_download)} proteins")

        downloaded_uniprot_data, all_ecs = get_uniprot_data(gbk_data_to_download, cache_dir, args)

        uniprot_dict.update(downloaded_uniprot_data)

    else:
        logger.warning(f"Using only data from cache:\n{args.use_uniprot_cache}")

    if len(list(downloaded_uniprot_data.keys())) != 0:
        cache_uniprot_data(uniprot_dict, cache_dir, time_stamp)

    # get genbank accessions by mapping uniprot accession to ncbi
    # if can't get genbank accession and did not retrieve a gene name from uniprot
    # map the uniprot acc to the gene name and then retrieve the genbank accession from ncbi
    uniprot_dict = get_mapped_genbank_accessions(uniprot_dict, cache_dir, args)

    acc_to_remove = set()
    for uniprot_acc in uniprot_dict:
        try:
            uniprot_dict[uniprot_acc]['genbank_accession']
        except KeyError:
            logger.error(
                f"Could not map the UniProt accession '{uniprot_acc}' to a GenBank accession\n"
                "directly via the UniProt mapping service or via its gene name.\n"
                f"Not adding protein data for the UniProt accession '{uniprot_acc}' to the\n"
                "local CAZyme database."
            )
            acc_to_remove.add(uniprot_acc)
    for uniprot_acc in acc_to_remove:
        del uniprot_dict[uniprot_acc]

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

    logger = logging.getLogger(__name__)

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

                for mapping result in mapping_results:
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


def get_mapped_genbank_accessions(uniprot_dict, cache_dir, args):
    """Map uniprot accessions to GenBank protein version accessions.
    
    :param uniprot_dict: {uniprot_acc: {gene_name: str, protein_name: str, pdb: set, ec: set, sequence:str, seq_data:str}}
    :param cache_dir: path to cache directory
    :param args: CLI parser

    Return uniprot_dict with GenBank accessions
    """
    logger = logging.getLogger(__name__)
    mapping_batch_size = 25
    cache_path = cache_dir / "failed_mapped_uniprot_ids"

    bioservices_queries = get_chunks_list(
        list(uniprot_dict.keys()),
        mapping_batch_size,
    )

    failed_ids = set()

    for batch in tqdm(bioservices_queries, desc="Mapping UniProt acc to GenBank acc"):
        mapping_dict = UniProt().mapping(
            fr="UniProtKB_AC-ID",
            to="EMBL-GenBank-DDBJ_CDS",
            query=batch,
        )
        try:
            for mapped_pair in mapping_dict['results']:
                uniprot_acc = mapped_pair['from']
                gbk_acc = mapped_pair['to']

                try:
                    uniprot_dict[uniprot_acc]['genbank_accession'] = gbk_acc
                except KeyError:
                    logger.error(
                        f"Retrieved UniProt accessions {uniprot_acc} from UniProt but accession was\n"
                        "not retrieved when quering by GenBank accession"
                    )
                    pass
        except KeyError:
            pass

        try:
            for acc in mapping_dict['failedIds']:
                if acc in (list(uniprot_dict.keys())):
                    failed_ids.add(acc)
        except KeyError:
            pass

    if len(failed_ids) != 0:
        logger.warning(f"Could not map {len(failed_ids)} UniProt accessions to GenBank accessions")

        failed_ids_to_parse = set()

        for failed_id in failed_ids:
            if uniprot_dict[failed_id]['gene_name'] == 'nan':
                failed_ids_to_parse.add(failed_id)

        if len(failed_ids_to_parse) != 0:
            failed_ids = set()

            logger.warning(
                f"Could not map {len(failed_ids_to_parse)} UniProt accessions to GenBank accessions\n"
                "and could did not retrieve a gene name from UniProt.\n"
                "Will try mapping the UniProt acc to the gene name"
            )

            # retry getting gene names if did not have gene names before
            bioservices_queries = get_chunks_list(
                list(failed_ids_to_parse),
                mapping_batch_size,
            )

            for batch in tqdm(bioservices_queries, desc="Getting gene names"):
                mapping_dict = UniProt().mapping(
                    fr="UniProtKB_AC-ID",
                    to="EMBL-GenBank-DDBJ",
                    query=batch,
                )

                try:
                    for mapped_pair in mapping_dict['results']:
                        uniprot_acc = mapped_pair['from']
                        gene_name = mapped_pair['to']

                        try:
                            uniprot_dict[uniprot_acc]['gene_name'] = gene_name
                        except KeyError:
                            logger.error(
                                f"Retrieved UniProt accessions {uniprot_acc} from UniProt but accession was\n"
                                "not retrieved when quering by GenBank accession"
                            )
                            pass
                except KeyError:
                    pass

                try:
                    for acc in mapping_dict['failedIds']:
                        if acc in list(uniprot_dict.keys()):
                            failed_ids.add(acc)
                except KeyError:
                    pass

    if len(failed_ids) != 0:
        with open(cache_path, "w") as fh:
            for acc in failed_ids:
                if acc in list(uniprot_dict.keys()):
                    if uniprot_dict[acc]['gene_name'] == 'nan':
                        logger.error(
                            f"Could not map the UniProt accession '{acc}' to a NCBI GenBank protein\n"
                            "accession or gene name, so can't map the UniProt accession to a GenBank\n"
                            "accession in the local CAZyme database.\n"
                            f"Not adding protein data for UniProt accession '{acc}' to the local\n"
                            "CAZyme database"
                        )
                        del uniprot_dict[acc]
                        fh.write(f"Could not map UniProt accession '{acc}' to a GenBank accession or gene name\n")

    uniprot_dict = get_linked_ncbi_accessions(uniprot_dict, args)
    
    return uniprot_dict
    


def get_gene_names(uniprot_dict, args):
    """Map UniProt accessions to Ensemble-GenBank gene name
    
    :param uniprot_dict: {uniprot_acc: {gene_name: str, protein_name: str, pdb: set, ec: set, sequence:str, seq_data:str}}
    :param args: CLI parser
    
    Return UniProt dict
    """
    uniprot_acc_to_parse = [acc for acc in uniprot_dict if uniprot_dict[acc]['gene_name'] == 'nan']

    bioservices_queries = get_chunks_list(
        uniprot_acc_to_parse,
        args.uniprot_batch_size,
    )

    for batch in tqdm(bioservices_queries, "Getting gene names from UniProt"):
        mapping_dict = UniProt().mapping(
            fr="UniProtKB_AC-ID",
            to="EMBL-GenBank-DDBJ",
            query=batch,
        )
        for result in mapping_dict['results']:
            uniprot_acc = result['from']
            gene_name = result['to']
            uniprot_dict[uniprot_acc]['gene_name'] = gene_name
    
    return uniprot_dict


if __name__ == "__main__":
    main()
