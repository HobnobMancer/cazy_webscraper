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

from datetime import datetime
from typing import List, Optional

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from saintBioutils.genbank import entrez_retry
from saintBioutils.misc import get_chunks_list
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from cazy_webscraper import closing_message, connect_existing_db
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
        cache_dir = cache_dir / "genbank_data_retrieval"
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

    genbank_accessions = list(gbk_dict.keys())
    logger.warning(f"Retrieving GenBank sequences for {len(gbk_dict.keys())}")  

    if args.seq_dict:
        logger.warning(f"Getting sequences from cache: {args.seq_dict}")
        with open(args.seq_dict, "r") as fh:
            cache_dict = json.load(fh)

        # convert strs to SeqRecords
        seq_dict = {}
        for key in cache_dict:
            seq_dict[key] = Seq(cache_dict[key])

    else:
        seq_dict, no_seq = get_sequences(genbank_accessions, args)  # {gbk_accession: seq}

        # only cache the sequence. Seq obj is not JSON serializable
        cache_dict = {}
        for key in seq_dict:
            cache_dict[key] = str(seq_dict[key])

        # cache the retrieved sequences
        cache_path = cache_dir / f"genbank_seqs_{time_stamp}.json"
        with open(cache_path, "w") as fh:
            json.dump(cache_dict, fh)

        if len(no_seq) != 0:
            no_seq_cache = cache_dir / f"no_seq_retrieved_{time_stamp}.txt"
            logger.warning(
                f"No protein sequence retrieved for {len(no_seq)}\n"
                f"The GenBank accessions for these proteins have been written to: {no_seq_cache}"
            )
            with open(no_seq_cache, "a") as fh:
                for acc in no_seq:
                    fh.write(f"{acc}\n")

    logger.warning(f"Adding {len(list(seq_dict.keys()))} protein seqs to the db")

    seq_records = []
    for acc in seq_dict:
        sequence = Seq(seq_dict[acc])
        record = SeqRecord(sequence, id=acc)
        seq_records.append(record)
    
    SeqIO.write(seq_records, "ce_gbk_protein_seqs.fasta", "fasta")

    add_genbank_data.add_gbk_seqs_to_db(seq_dict, date_today, gbk_dict, connection, args)

    closing_message("get_genbank_sequences", start_time, args)


def get_sequences(genbank_accessions, args, retry=False):
    """Retrieve protein sequences from Entrez.

    :param genbank_accessions: list, GenBank accessions
    :param args: cmb-line args parser
    :param retry: bool, default False, if get_sequences is being called for retrying a previously failed query

    Return dict keyed by GenBank accession and valued by Seq instance, and a list of all GenBank accessions 
    for which no record from NCBI was retrieved.
    """
    logger = logging.getLogger(__name__)

    seq_dict = {}  # {gbk_accession: SeqRecord}

    # the list of accessions is to long, break down into smaller chunks for batch querying
    all_queries = get_chunks_list(genbank_accessions, args.batch_size)

    failed_queries = []  # lists which raised an error, likely because contain an accession not in NCBI

    irregular_accessions = []

    success_accessions = set()  # accessions for which seqs were retrieved

    for query_list in tqdm(all_queries, desc="Batch querying NCBI.Entrez"):

        try:
            epost_webenv, epost_query_key = bulk_query_ncbi(query_list, args)
        except RuntimeError:
            logger.warning(
                "Runtime error raised when batch quering\n"
                "Possible result of a accessions not being in NCBI\n"
                "Attempt identification of the causal accession later\n"
            )

            if retry:
                return None, None

            failed_queries.append(query_list)
            continue

        try:
            # retrieve the protein sequences
            with entrez_retry(
                args.retries,
                Entrez.efetch,
                db="Protein",
                query_key=epost_query_key,
                WebEnv=epost_webenv,
                rettype="fasta",
                retmode="text",
            ) as seq_handle:
                for record in SeqIO.parse(seq_handle, "fasta"):
                    temp_accession = record.id

                    # check if multiple items returned in ID
                    temp_accession = temp_accession.split("|")
                    retrieved_accession = None

                    for acc in temp_accession:
                        if acc.strip() in genbank_accessions:
                            retrieved_accession = acc

                    if retrieved_accession is None:  # if could not retrieve GenBank accession from the record
                        logger.error(
                            "Could not retrieve a GenBank protein accession matching an accession from the local database from:\n"
                            f"{record.id}\n"
                            "The sequence from this record will not be added to the db"    
                        )
                        irregular_accessions.append(temp_accession)
                        continue
                        
                    seq_dict[retrieved_accession] = record.seq

                    success_accessions.add(retrieved_accession)
        
        except IncompleteRead as err:
            logger.warning(
                "IncompleteRead error raised:\n"
                f"{err}\n"
                "Will reattempt NCBI query later"
            )

            if retry:
                return None, None

            failed_queries.append(all_queries)
            continue

    # list of GenBank accessions for which no protein sequence was retrieved

    no_seq = [acc for acc in genbank_accessions if acc not in success_accessions]
    no_seq += irregular_accessions

    if retry:
        return seq_dict, no_seq

    if len(failed_queries) != 0:
        for failed_query in tqdm(failed_queries, desc="Reparsing failed queries"):
            first_half = failed_query[:int((len(failed_query)/2))]

            seq_dict, success_accessions, failed_accessions = retry_failed_queries(
                first_half,
                seq_dict,
                success_accessions,
                args,
            )

            no_seq += failed_accessions

            second_half = failed_query[int((len(failed_query)/2)):]

            seq_dict, success_accessions, failed_accessions = retry_failed_queries(
                second_half,
                seq_dict,
                success_accessions,
                args,
            )

            no_seq += failed_accessions

    logger.warning(f"Retrieved sequences for {len(success_accessions)} proteins")

    return seq_dict, no_seq


def bulk_query_ncbi(accessions, args):
    """Bulk query NCBI and retrieve webenvironment and history tags

    :param accessions: list of GenBank protein accessions
    :param args: cmd-line args parser

    Return webenv and query key
    """
    # perform batch query of Entrez
    try:
        accessions_string = ",".join(accessions)
    except TypeError:
        accessions_string = accessions

    # Runtime error captured by try/except function call
    epost_result = Entrez.read(
        entrez_retry(
            args.retries,
            Entrez.epost,
            db="Protein",
            id=accessions_string,
        )
    )

    # retrieve the web environment and query key from the Entrez post
    epost_webenv = epost_result["WebEnv"]
    epost_query_key = epost_result["QueryKey"]

    return epost_webenv, epost_query_key


def retry_failed_queries(query, seq_dict, success_accessions, args):
    """Parse queries which previously raised a Runtime Error, likey the result of containing accessions 
    that are not present in NCBI.
    
    :param query: list GenBank accessions
    :param seq_dict: dict, keyed by GenBank accession, valued by BioPython Seq obj
    :param success_acessions: list of GenBank accessions for which a protein sequence was retrieved
    :param args: cmd-line args parser
    
    Return seq_dict and success_accessions with data from the reparsed failed queries
    """
    failed_accessions = []  # accessions of proteins for which no record was retrieved from NCBI.Entrez

    new_seq_dict, no_seq = get_sequences(query, args, retry=True)
    
    if new_seq_dict is None:  # failed retrieval of data from NCBI
        for acc in query:
            new_seq_dict, no_seq = get_sequences([acc], args, retry=True)
        
            if new_seq_dict is not None:  # retrieve data for this accession
                seq_dict[acc] = new_seq_dict[acc]

            else:  # failed to retrieve data for this accession
                failed_accessions.append(acc)
    
    else:  # successful retrieval of data from NCBI
        for acc in new_seq_dict:
            seq_dict[acc] = new_seq_dict[acc]
            success_accessions.add(acc)

        failed_accessions += no_seq

    return seq_dict, success_accessions, failed_accessions


if __name__ == "__main__":
    main()
