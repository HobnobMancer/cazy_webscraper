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

    seq_dict = get_cache_seqs(start_time, args)
    gbk_dict = {}

    if args.file_only:
        logger.warning(
            "Only adding Seqs in JSON and/or FASTA file.\n"
            "Not retrieving seqs from NCBI\n"
        )

    else:  # retrieve data from NCBI for seqs in the local db
        seq_dict, gbk_dict = get_seqs_from_ncbi(
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

    logger.warning(f"Adding {len(list(seq_dict.keys()))} to the local CAZyme db")

    # get acc and db ids for seqs to be added to the db from the cache
    acc_with_unknown_db_id = list(set(seq_dict.keys()) - set(gbk_dict.keys()))
    if len(acc_with_unknown_db_id) > 0:
        gbk_dict.update(get_selected_gbks.get_ids(acc_with_unknown_db_id, connection))

    try:
        shutil.make_archive(cache_dir, 'tar',  cache_dir)
        shutil.rmtree(cache_dir)
    except FileNotFoundError:
        logger.error(
            "Failed to compress cache directory"
        )

    if len(list(seq_dict.keys())) != 0:
        closing_message()

    add_genbank_data.add_gbk_seqs_to_db(seq_dict, date_today, gbk_dict, connection, args)

    closing_message("Get GenBank Sequences", start_time, args)


def get_cache_seqs(start_time, args):
    """Extract protein sequences from FASTA and/or JSON file, which will be added to the
    local CAZyme database

    :param seq_dict: dict, {genbank_acc: Bio.Seq}

    Return update seq_dict
    """
    logger = logging.getLogger(__name__)

    seq_dict = {}

    if args.seq_dict:
        logger.warning(f"Getting sequences from JSON cache:\n{args.seq_dict}")

        try:
            with open(args.seq_dict, "r") as fh:
                cache_dict = json.load(fh)

        except FileNotFoundError:
            logger.error(
                f"Could not find JSON file of protein sequences at:\n"
                f"{args.seq_dict}\n"
                "Check the path is correct"
                "Terminating program"
            )
            closing_message("Get GenBank seqs", start_time, args, early_term=True)

        # convert strs to SeqRecords
        for key in cache_dict:
            seq_dict[key] = Seq(cache_dict[key])

    if args.seq_file:
        logger.warning(f"Getting sequences from FASTA cache:\n{args.seq_file}")

        try:
            for record in SeqIO.parse(args.seq_file, "fasta"):
                try:
                    seq_dict[record.id]
                    if seq_dict[record.id] != record.seq:
                        logger.warning(
                            f"Retrieved seq for {record.id} from JSON file which does NOT match "
                            "the seq in the FASTA file.\n"
                            "Adding seq from the FASTA file to the local CAZyme database\n"
                            f"JSON seq: {seq_dict[record.id]}\n"
                            f"FASTA seq: {record.seq}"
                        )
                        seq_dict[record.id] = record.seq
                except KeyError:
                    seq_dict[record.id] = record.seq

        except FileNotFoundError:
            logger.error(
                f"Could not find FASTA file of protein sequences at:\n"
                f"{args.seq_file}\n"
                "Check the path is correct"
                "Terminating program"
            )
            closing_message("Get GenBank seqs", start_time, args, early_term=True)

    return seq_dict


def get_records_to_retrieve(
    class_filters,
    family_filters,
    taxonomy_filter_dict,
    kingdom_filters,
    ec_filters,
    seq_dict,
    start_time,
    connection,
    args,
):
    """Get Genbank accessions to retrieved data from NCBI for

    :param seq_dict: dict {id: seq} of seqs retrieved from cache json/fasta files
    :param start_time: str: time program was called
    :param connection: open connection to a SQLite db engine
    :param args: CLI args parser

    Return a dict {gbk_acc: gbk db id
    """
    logger = logging.getLogger(__name__)

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

        accessions = [line.strip() for line in lines if line not in list(seq_dict.keys())]
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

        # don't retrieve data for cached seqs
        for gbk_acc in gbk_dict:
            if gbk_acc in list(seq_dict.keys()):
                del gbk_dict[gbk_acc]

    return gbk_dict


def get_seqs_from_ncbi(
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
    """Coordinate retrieving sequence data from NCBI for proteins not retrieved from cache files

    :param seq_dict: dict {id: seq} of seqs retrieved from cache json/fasta files
    :param start_time: str: time program was called
    :param connection: open connection to a SQLite db engine
    :param args: CLI args parser

    Return dicts of {acc: Bio.Seq} and {gbk acc: db id}
    """
    logger = logging.getLogger(__name__)

    gbk_dict = get_records_to_retrieve(
        class_filters,
        family_filters,
        taxonomy_filter_dict,
        kingdom_filters,
        ec_filters,
        seq_dict,
        start_time,
        connection,
        args,
    )

    genbank_accessions = list(gbk_dict.keys())
    logger.warning(f"Retrieving GenBank sequences from NCBI for {len(gbk_dict.keys())}")

    if len(genbank_accessions) == 0:
        return seq_dict

    cache_path = cache_dir / f"genbank_seqs_{start_time}.fasta"

    # break up long list into managable chunks
    all_queries = get_chunks_list(genbank_accessions, args.batch_size)

    # list of downloaded SeqRecords and list of gbk acc for which no seq was retrieved from NCBI
    downloaded_seqs, failed_queries = get_sequences(all_queries, cache_dir, args)

    # retry failed accs
    no_seq_acc = []
    if len(failed_queries) != 0:
        # break up and query individually
        retrying_acc = {}  # {acc: # of tries}
        for batch in failed_queries:
            for acc in batch:
                retrying_acc[acc] = 0

        finished_retry = False

        acc_list = list(retrying_acc)

        no_seq_acc = set()  # set of accessions for which no seq could be retrieved

        while finished_retry is False:
            acc_list = list(retrying_acc)

            for accession in acc_list:
                new_seq, failed_seq = get_sequences([[accession]], cache_dir, args)

                if len(new_seq) == 0:
                    retrying_acc[accession] += 1

                    if retrying_acc[accession] >= args.retries:
                        del retrying_acc[accession]
                        no_seq_acc.add(accession)

            acc_list = list(retrying_acc)

            if len(acc_list) > 0:
                finished_retry = True

    # cache accs of proteins for which not seq could be retrieved from NCBI
    if len(no_seq_acc) != 0:
        no_seq_cache = cache_dir / f"no_seq_retrieved_{start_time}.txt"

        logger.warning(
            f"No protein sequence retrieved for {len(no_seq_acc)} proteins\n"
            f"The GenBank accessions for these proteins have been written to: {no_seq_cache}"
        )

        try:
            with open(no_seq_cache, "a") as fh:
                for acc in no_seq_acc:
                    fh.write(f"{acc}\n")
        except FileNotFoundError:
            logger.error(
                "Could not cache acc of proteins for which not seqs were retrieved from NCBI to:\n"
                f"{no_seq_cache}"
            )

    # only cache the sequence. Seq obj is not JSON serializable
    cache_dict = {}
    for key in seq_dict:
        cache_dict[key] = str(seq_dict[key])

    # cache the retrieved sequences
    cache_path = cache_dir / f"genbank_seqs_{start_time}.json"
    with open(cache_path, "w") as fh:
        json.dump(cache_dict, fh)

    for record in downloaded_seqs:
        seq_dict[record.id] = record.seq

    return seq_dict, gbk_dict


def get_sequences(batches, cache_dir, args):
    """Retrieve protein sequences from Entrez.

    :param batches: list of lists, one list be batch of gbk acc to query against NCBI
    :param cache_dir: Path, to cache directory
    :param args: cmb-line args parser
    failed query

    Return
    * list of SeqRecords
    * list of failed batches
    """
    logger = logging.getLogger(__name__)

    cache_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    cache_path = cache_dir / f"genbank_seqs_{cache_time}.json"

    gbk_acc_to_retrieve = []
    for batch in batches:
        gbk_acc_to_retrieve += [acc for acc in batch]

    seq_records = []  # SeqRecords, can't be set, SeqRecords are unhasable

    failed_queries = []  # lists which raised an error,
    # likely because contain an accession not in NCBI

    irregular_accessions = []

    success_accessions = set()  # accessions for which seqs were retrieved

    for batch in tqdm(batches, desc="Batch querying NCBI.Entrez"):

        # POST IDS
        try:
            epost_webenv, epost_query_key = bulk_query_ncbi(batch, args)
        except RuntimeError:
            logger.warning(
                "Runtime error raised when batch quering\n"
                "Possible result of a accessions not being in NCBI\n"
                "Attempt identification of the causal accession later\n"
            )

            failed_queries.append(batch)
            continue

        if epost_webenv is None:
            failed_queries.append(batch)
            continue

        # Retrieve seqs
        try:
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
                    try:
                        retrieved_accession = [_ for _ in record.id.split("|") if _.strip() in gbk_acc_to_retrieve][0]
                    except IndexError:
                        logger.error(
                            "Could not retrieve a GenBank protein accession matching "
                            "an accession from the local database from:\n"
                            f"{record.id}\n"
                            "The sequence from this record will not be added to the db"
                        )
                        irregular_accessions.append(temp_accession)
                        continue

                    if retrieved_accession not in success_accessions:
                        seq_records.append(record)
                        success_accessions.add(retrieved_accession)

        except IncompleteRead as err:
            logger.warning(
                "IncompleteRead error raised:\n"
                f"{err}\n"
                "Will reattempt NCBI query later"
            )

            failed_queries.append(batch)
            continue

        except Exception as err:
            logger.warning(
                "Error raised:\n"
                f"{err}\n"
                "Will reattempt NCBI query later"
            )

            failed_queries.append(batch)
            continue

        # cache the currently retrieved seqs
        SeqIO.write(seq_records, cache_path, "fasta")

    # list of GenBank accessions for which no protein sequence was retrieved

    return seq_records, failed_queries


def bulk_query_ncbi(accessions, args):
    """Bulk query NCBI and retrieve webenvironment and history tags

    :param accessions: list of GenBank protein accessions
    :param args: cmd-line args parser

    Return webenv and query key
    """
    logger = logging.getLogger(__name__)

    # perform batch query of Entrez
    try:
        accessions_string = ",".join(accessions)
    except TypeError:
        accessions_string = accessions

    # Runtime error captured by try/except function call
    try:
        epost_result = Entrez.read(
            entrez_retry(
                args.retries,
                Entrez.epost,
                db="Protein",
                id=accessions_string,
            ),
            validate=False,
        )
    except NotXMLError:
        logger.error("Could not parse Entrez output")
        return None, None

    # retrieve the web environment and query key from the Entrez post
    epost_webenv = epost_result["WebEnv"]
    epost_query_key = epost_result["QueryKey"]

    return epost_webenv, epost_query_key


if __name__ == "__main__":
    main()
