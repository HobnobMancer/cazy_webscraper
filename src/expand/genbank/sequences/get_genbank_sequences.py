#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
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
"""Retrieve proteins sequences from GenBank and populate the local database"""


from http.client import IncompleteRead
import json
import logging
import pandas as pd
import re
import shutil

from datetime import datetime
from typing import List, Optional

from Bio import Entrez, SeqIO
from Bio.Entrez.Parser import NotXMLError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from saintBioutils.genbank import entrez_retry
from saintBioutils.misc import get_chunks_list
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from cazy_webscraper import closing_message, connect_existing_db
from cazy_webscraper.cache.ncbi import get_cache_seqs
from cazy_webscraper.ncbi.sequences import post_accessions_to_entrez, fetch_ncbi_seqs, get_protein_accession
from cazy_webscraper.sql import sql_orm, sql_interface
from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks
from cazy_webscraper.sql.sql_interface.add_data import add_genbank_data
from cazy_webscraper.utilities.parse_configuration import get_expansion_configuration
from cazy_webscraper.utilities.parsers.gbk_seq_parser import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)
    date_today = datetime.now().strftime("%Y-%m-%d")  # used as seq_update_date in the db

    # parse cmd-line arguments
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        args = build_parser(argv).parse_args()

    if logger is None:
        logger = logging.getLogger(__package__)
        config_logger(args)

    logger.info("Providing user email address to NCBI.Entrez")
    Entrez.email = args.email

    if args.seq_update:
        logger.warning("Enabled updating sequences")

    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "genbank_seq_retrieval"
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
    with sql_orm.Session(bind=connection) as session:
        retrieved_data = "GenBank protein sequences"
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            'GenBank',
            retrieved_data,
            session,
            args,
        )

    seq_dict, seq_records = get_cache_seqs(start_time, args)
    # seq_dict = {genbank_acc: Bio.Seq}
    # seq_records = list of BioPython SeqRecords

    # Get GenBank accessions from a file or records matching the provided criteria, and get list
    # of genbank accesions for proteins for whom seqs need to be downloaded from NCBI
    # gbk_dict = {gbk_acc: db id}
    # accs_seqs_to_retrieve = [list of gbk acc to query ncbi with]
    gbk_dict, accs_seqs_to_retrieve = get_records_to_retrieve(
        class_filters,
        family_filters,
        taxonomy_filter_dict,
        kingdom_filters,
        ec_filters,
        seq_dict,
        start_time,
        connection,
        cache_dir,
        args,
    )

    if args.file_only:
        logger.warning(
            "Only adding Seqs in JSON and/or FASTA file.\n"
            "Not retrieving seqs from NCBI\n"
        )

    else:  # retrieve data from NCBI for seqs in the local db
        seq_dict = get_seqs_from_ncbi(
            seq_records,
            accs_seqs_to_retrieve,
            seq_dict,
            start_time,
            cache_dir,
            args,
        )

    logger.warning(f"Adding {len(list(seq_dict.keys()))} to the local CAZyme db")

    # get acc and db ids for seqs to be added to the db from the cache
    acc_with_unknown_db_id = list(set(seq_dict.keys()) - set(gbk_dict.keys()))
    if len(acc_with_unknown_db_id) > 0:
        gbk_dict.update(get_selected_gbks.get_ids(acc_with_unknown_db_id, connection, cache_dir))

    try:
        shutil.make_archive(cache_dir, 'tar',  cache_dir)
        shutil.rmtree(cache_dir)
    except FileNotFoundError:
        logger.error(
            "Failed to compress cache directory"
        )

    if len(list(seq_dict.keys())) == 0:
        logger.warning("Adding 0 sequences to the local CAZyme database")
        closing_message("Get GenBank Sequences", start_time, args)

    # make subdir in cache dir for logging failed seqs
    make_output_directory(cache_dir, force=True, nodelete=True)
    add_genbank_data.add_gbk_seqs_to_db(seq_dict, date_today, gbk_dict, connection, cache_dir, args)

    closing_message("Get GenBank Sequences", start_time, args)


def get_records_to_retrieve(
    class_filters,
    family_filters,
    taxonomy_filter_dict,
    kingdom_filters,
    ec_filters,
    seq_dict,
    start_time,
    connection,
    cache_dir,
    args,
):
    """Get Genbank accessions to retrieved data for from NCBI and genbank accessions and local 
    database IDs for the GenBank accessions. Remove seqs that are not in the local db.

    :param seq_dict: dict {id: seq} of seqs retrieved from cache json/fasta files
    :param start_time: str: time program was called
    :param connection: open connection to a SQLite db engine
    :param cache_dir: Path to cache directory
    :param args: CLI args parser

    Return a dict {gbk_acc: gbk db id} and list of genbank accs to retrieve seqs for.
    """
    logger = logging.getLogger(__name__)

    all_accessions = list(seq_dict.keys())  # cache + (file or db)

    if args.file_only:
        gbk_dict = get_selected_gbks.get_ids(all_accessions, connection)
        accs_to_retrieve = []
        return gbk_dict, accs_to_retrieve

    accessions_of_interest = []  # acc from file or from db that match provided criteria

    # retrieve dict of genbank accession and genbank db ids from the local CAZyme db
    if args.genbank_accessions is not None:
        logger.warning(f"Getting GenBank accessions from file: {args.genbank_accessions}")

        try:
            with open(args.genbank_accessions, "r") as fh:
                lines = fh.read().splitlines()
        except FileNotFoundError:
            logging.error(
                "Could not find list of GenBank accessions at:\n"
                f"{args.genbank_accessions}"
                "Check the path is correct\n"
                "Terminating program"
            )
            closing_message("Get GenBank seqs", start_time, args, early_term=True)

        # don't add accessions that are already in the cache
        accessions_of_interest = [line.strip() for line in lines if line not in list(seq_dict.keys())]
        all_accessions += list(set(accessions_of_interest))

    else:
        logger.warning("Getting GenBank accessions of proteins matching the provided criteria from the local CAZyme db")
        gbk_dict = get_selected_gbks.get_genbank_accessions(
            class_filters,
            family_filters,
            taxonomy_filter_dict,
            kingdom_filters,
            ec_filters,
            connection,
        )

        # don't get seq for accession retrieved from cache
        accessions_of_interest = [acc for acc in gbk_dict if acc not in list(seq_dict.keys())]
        all_accessions += accessions_of_interest

    gbk_dict = get_selected_gbks.get_ids(set(all_accessions), connection, cache_dir)

    return gbk_dict, accessions_of_interest


def get_seqs_from_ncbi(
    seq_records,
    accs_seqs_to_retrieve,
    seq_dict,
    start_time,
    cache_dir,
    args,
):
    """Coordinate retrieving sequence data from NCBI for proteins not retrieved from cache files

    :param seq_records: list of BioPython SeqRecords
    :param accs_seqs_to_retrieve: list of gbk accs of protein seqs to retrieve from ncbi
    :param seq_dict: dict {id: seq} of seqs retrieved from cache json/fasta files
    :param start_time: str: time program was called
    :param cache_dir: path to cache directory
    :param args: CLI args parser

    Return updated seqs_dict, {acc: Bio.Seq}, containing seqs retrieved from NCBI
    """
    logger = logging.getLogger(__name__)

    logger.warning(f"Retrieving GenBank sequences from NCBI for {len(accs_seqs_to_retrieve)}")

    if len(accs_seqs_to_retrieve) == 0:
        return seq_dict

    cache_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    cache_path = cache_dir / f"cw_genbank_seqs_{cache_time}.fasta"
    invalid_ids_cache = cache_dir / "invalid_ids"

    # break up long list into managable chunks
    batches = get_chunks_list(accs_seqs_to_retrieve, args.batch_size)

    (
        seq_records,
        batches_with_invalid_ids,
        failed_connections_batches,
        accs_still_to_fetch,
    ) = get_seqs(seq_records, batches, cache_path, invalid_ids_cache, args)

    if len(list(failed_connections_batches.keys())) > 0:
        # Retry failed batches before parsing invalid IDs as the failed connections
        # may contain invalid IDs
        (
            seq_records,
            batches_with_invalid_ids,
        ) = parse_failed_connections(
            failed_connections_batches,
            seq_records,
            accs_still_to_fetch,
            cache_dir,
            cache_path,
            batches_with_invalid_ids,
            invalid_ids_cache,
            args,
        )

    if len(batches_with_invalid_ids) > 0:
        ## retry IDs individually and identify invalid IDs
        # break up batches into batches of len 1 (1 acc per list/batch)
        individual_batches = []
        for batch in batches_with_invalid_ids:
            for acc in batch:
                individual_batches.append([acc])

        seq_records = parse_invalid_id_batches(
            seq_records,
            batches_with_invalid_ids,
            cache_path,
            invalid_ids_cache,
            args,
        )

    # add seq records to seq dict
    for record in seq_records:
        try:
            seq = seq_dict[record.id]
            if seq != record.seq:
                logger.warning(
                    f"Sequence downloaded from NCBI for protein '{record.id}' does not match the seq\n"
                    "retrieved from the cache:\n"
                    f"Downloaded seq:\n{record.seq}\nSeq from cache:\n{seq}\n"
                    "Overwriting seq from cache and using downloaded seq"
                )
            seq_dict[record.id] = record.seq
        except KeyError:
            seq_dict[record.id] = record.seq

    accs_still_to_fetch = set()
    for acc in accs_seqs_to_retrieve:
        try:
            seq_dict[acc]
        except KeyError:
            accs_still_to_fetch.add(acc)

    if len(accs_still_to_fetch) > 0:
        logger.error(
            f"Could not fetch seqs for {len(accs_still_to_fetch)} NCBI protein accessions.\n"
            "Caching these protein accessions"
        )
        failed_ids_cache = cache_dir / "failed_retrieve_ids"
        with open(failed_ids_cache, "a") as fh:
            for acc in accs_still_to_fetch:
                fh.write(f"{acc}\n")

    logger.warning(f"Retrieved {len(seq_records)} from NCBI")

    return seq_dict


def get_seqs(seq_records, batches, cache_path, invalid_ids_cache, args, disable_prg_bar=False):
    """Retrieve protein sequences from Entrez.

    :param seq_records: list of Bio SeqRecords
    :param batches: list of lists, one list be batch of gbk acc to query against NCBI
    :param cache_path: Path, to cache fasta file or protein seqs
    :param invalid_id_cache: Path to file to cache invalid IDs
    :param args: cmb-line args parser
    failed query

    Return
    * list of SeqRecords
    * list of nested lists of batches containing one or more invalid IDs (i.e. no longer listd in NCBI)
    * dict of batches which suffered failed connections when querying NCBI
    * list of accessions for whom a seq has not been retrieved yet
    """
    logger = logging.getLogger(__name__)

    gbk_acc_to_retrieve = []
    for batch in batches:
        gbk_acc_to_retrieve += [acc for acc in batch]

    invalid_ids = set()
    batches_with_invalid_ids = []  # nested lists of batches with invalid IDs
    failed_connections_batches = {}  # dict of batches during the querying of which the connection failed

    successful_accessions = []  # accessions for which seqs were retrieved

    for batch in tqdm(batches, desc="Batch querying NCBI.Entrez", disable=disable_prg_bar):
        epost_webenv, epost_query_key, success = post_accessions_to_entrez(batch, args)

        if success == "Invalid ID":
            # invalid ID, cache, and do not query again
            invalid_ids.add(batch[0])
            with open(invalid_ids_cache, "a") as fh:
                fh.write(f"{batch[0]}\n")
            continue

        elif success == "Contains invalid ID":
            # batch contains one or more invalid IDs.
            # Don't know which are invalid so will try later
            for acc in batch:
                batches_with_invalid_ids.append(acc)
            continue

        elif success == "Failed connection":
            # Connection failed so retry later
            failed_connections_batches["_".join(batch)] = {
                'batch': batch,
                '#_of_attempts': 1 
            }
            continue
        
        # else was a success, so proceed to retrieving protein seqs
        seq_records, success, temp_successful_accessions = fetch_ncbi_seqs(
            seq_records,
            epost_webenv,
            epost_query_key,
            gbk_acc_to_retrieve,
            args,
        )
        
        if success == "Invalid ID":
            invalid_ids.add(batch[0])
            with open(invalid_ids_cache, "a") as fh:
                fh.write(f"{batch[0]}\n")
            continue

        elif success == "Contains invalid ID":
            for acc in batch:
                batches_with_invalid_ids.append(acc)
            continue

        elif success == "Failed connection":
            failed_connections_batches["_".join(batch)] = {
                'batch': batch,
                '#_of_attempts': 1 
            }
            continue

        SeqIO.write(seq_records, cache_path, "fasta")
        successful_accessions += temp_successful_accessions

    accs_still_to_fetch = [acc for acc in gbk_acc_to_retrieve if ((acc not in invalid_ids) or (acc not in successful_accessions))]

    return seq_records, batches_with_invalid_ids, failed_connections_batches, accs_still_to_fetch


def parse_failed_connections(
    failed_connections_batches,
    seq_records,
    accs_to_fetch,
    cache_dir,
    cache_path,
    batches_with_invalid_ids,
    invalid_ids_cache, 
    args,
):
    """Parse batches that suffered failed connections on the first attempt.

    :param failed_connections_batches: dict of batches, key str(batch): {'batch': [], '#_of_attempts': int}
    :param seq_records: list of Bio.SeqRecords
    :param accs_to_fetch: list of NCBI protein accs to retrieve seqs for
    :param cache_path: path to cache downloaded seqs
    :param batches_with_invalid_ids: list of batches containing invalid IDs
    :param invalid_ids_cache: path to cache invalid IDs
    :param args: CLI args parser

    Return 
    * updated list of seq_records
    * batches with invalid IDs
    """
    logger = logging.getLogger(__name__)
    failed_connection_cache = cache_dir / "failed_entrez_connection_accessions"
    failed_batches = []
    for batch_name in failed_connections_batches:
        failed_batches.append(failed_connections_batches[batch_name]['batch'])

    while len(list(failed_connections_batches.keys())) != 0:
        # remove processed batches from failed_batches
        for batch in failed_batches:
            try:
                failed_connections_batches["_".join(batch)]
            except KeyError:
                # batch has been processed and no longer in failed_connections
                # delete batch from list
                failed_batches.remove(batch)

        if len(failed_batches) > 0:
            logger.warning(
                f"Retrying connection for {len(failed_batches)} remaining batches suffering failed connections"
            )
        
            # batch query failed batches
            (
                new_seq_records,
                new_batches_with_invalid_ids,
                new_failed_connections_batches,
                new_accs_still_to_fetch,
            ) = get_seqs(
                seq_records,
                failed_batches,
                cache_path,
                invalid_ids_cache,
                args,
            )

            seq_records += new_seq_records
            # remove successfully processed batches
            for batch in failed_batches:
                if (batch not in new_batches_with_invalid_ids) and (batch not in list(new_failed_connections_batches.keys())):
                    # not in invalid IDs or new failed batches, presume batch was processed
                    failed_batches.remove(batch)

            batches_with_invalid_ids += new_batches_with_invalid_ids
            # remove batches with invalid IDs, so don't retry connection
            for batch in batches_with_invalid_ids:
                failed_connections_batches["_".join(batch)]
                failed_batches.remove(batch)

            # remove batches that were processed successfully and 
            # increate attempt count for batches whose connections failed
            for batch_name in new_failed_connections_batches:
                if failed_connections_batches[batch_name]['#_of_attempts'] >= args.retries:
                    logger.error(
                        "Ran out of connection attempts for the following batch:\n"
                        f"{batch}\n"
                        "Will not add seqs for these proteins to the local CAZyme database."
                    )
                    with open(failed_connection_cache, "a") as fh:
                        for protein_acc in batch:
                            fh.write(f"{protein_acc}\n")
                    del failed_connections_batches[batch_name]

                else:
                    failed_connections_batches[batch_name]['#_of_attempts'] += 1

    return seq_records, batches_with_invalid_ids


def parse_invalid_id_batches(
    seq_records,
    batches_with_invalid_ids,
    cache_path,
    invalid_ids_cache,
    args,
):
    """Separate accessions in batches containing invalid IDs. Retrieve seqs for valid IDs.

    :param seq_records: list of Bio.SeqRecords
    :param batches_with_invalid_ids: list of accessions from batches containing invalid IDs
        NOT NESTED LISTS
    :param cache_path: path to cache downloaded seqs
    :param invalid_ids_cache: path to cache invalid IDs
    :param args: CLI args parser

    Return updated list of seq_records
    """
    logger = logging.getLogger(__name__)

    # turn into list of lists of len 1
    batches = [acc for acc in batches_with_invalid_ids]

    logger.warning(
        "Separted accessions in batches containing invalid IDs.\n"
        f"Processing {len(batches)} batches"
    )

    # inital parse
    (
        new_seq_records,
        batches_with_invalid_ids,
        failed_connections_batches,
        accs_still_to_fetch,
    ) = get_seqs(
        seq_records,
        batches,
        cache_path,
        invalid_ids_cache,
        args,
    )

    seq_records += new_seq_records

    if len(list(failed_connections_batches.keys())) > 0:
        # Retry failed batches before parsing invalid IDs as the failed connections
        # may contain invalid IDs
        (
            new_seq_records,
            batches_with_invalid_ids,
            accs_still_to_fetch,
        ) = parse_failed_connections(
            failed_connections_batches,
            seq_records,
            accs_still_to_fetch,
            cache_dir,
            cache_path,
            batches_with_invalid_ids.
            invalid_ids_cache, 
            args,
        )
        seq_records += new_seq_records

    return seq_records


if __name__ == "__main__":
    main()
