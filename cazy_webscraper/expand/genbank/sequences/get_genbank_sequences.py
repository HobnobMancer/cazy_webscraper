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

    seq_dict = get_cache_seqs(start_time, args)  # {genbank_acc: Bio.Seq}
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
        args,
    )

    if args.file_only:
        logger.warning(
            "Only adding Seqs in JSON and/or FASTA file.\n"
            "Not retrieving seqs from NCBI\n"
        )

    else:  # retrieve data from NCBI for seqs in the local db
        seq_dict = get_seqs_from_ncbi(
            accs_seqs_to_retrieve
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
    """Get Genbank accessions to retrieved data for from NCBI and genbank accessions and local 
    database IDs for the GenBank accessions. Remove seqs that are not in the local db.

    :param seq_dict: dict {id: seq} of seqs retrieved from cache json/fasta files
    :param start_time: str: time program was called
    :param connection: open connection to a SQLite db engine
    :param args: CLI args parser

    Return a dict {gbk_acc: gbk db id} and list of genbank accs to retrieve seqs for.
    """
    logger = logging.getLogger(__name__)

    all_accessions = list(seq_dict.keys())

    if args.file_only:
        gbk_dict = get_selected_gbks.get_ids(all_accessions, connection)
        accs_to_retrieve = []
        return gbk_dict, accs_seqs_to_retrieve

    accessions = []

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
        all_accessions += list(set(accessions))

    else:
        logger.warning("Getting GenBank accessions of proteins matching the provided criteria from the local CAZyme db")
        selected_gbk_dict = get_selected_gbks.get_genbank_accessions(
            class_filters,
            family_filters,
            taxonomy_filter_dict,
            kingdom_filters,
            ec_filters,
            connection,
        )

        accessions = list(gbk_dict.keys())
        all_accessions += accessions

        # don't retrieve data for cached seqs
        for gbk_acc in gbk_accs:
            if gbk_acc in list(seq_dict.keys()):
                del gbk_dict[gbk_acc]

    gbk_dict = get_selected_gbks.get_ids(all_accessions, connection)
    accs_to_retrieve = [for acc in all_accessions if acc not in list(seq_dict.keys())]

    return gbk_dict, accs_to_retrieve


def get_seqs_from_ncbi(
    accs_seqs_to_retrieve
    seq_dict,
    start_time,
    cache_dir,
    args,
):
    """Coordinate retrieving sequence data from NCBI for proteins not retrieved from cache files

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
    ) = get_seqs(batches, cache_path, invalid_ids_cache, args)

    if len(failed_connections_batches) > 0:
        # Retry failed batches before parsing invalid IDs as the failed connections
        # may contain invalid IDs
        (
            seq_records,
            batches_with_invalid_ids,
            accs_still_to_fetch,
        ) = parse_failed_connections(
            failed_connections_batches,
            seq_records,
            accs_still_to_fetch,
            cache_path,
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
        
        seq_records, accs_still_to_fetch = parse_invalid_id_batches(
            batches_with_invalid_ids,
            cache_path,
            invalid_ids_cache,
            args,
        )


        

    if len(accs_still_to_fetch) > 0:
        logger.error(
            f"Could not fetch seqs for {len(accs_still_to_fetch)} NCBI protein accessions.\n"
            "Caching these protein accessions"
        )
        failed_ids_cache = cache_dir / "failed_retrieve_ids"
        with open(failed_ids_cache, "a") as fh:
            for acc in accs_still_to_fetch:
                fh.write(f"{acc}\n")

    # add seq records to seq dict
    for record in seq_records:
        try:
            seq = seq_dict[record.id]
            if seq != record.seq:
                logger.warning(
                    f""
                )
            seq_dict[record.id] = record.seq
        except KeyError:
            seq_dict[record.id] = record.seq

    logger.warning(f"Retrieved {len(seq_records)} from NCBI")
    
    return seq_dict

    # list of downloaded SeqRecords and list of gbk acc for which no seq was retrieved from NCBI
    downloaded_seqs, failed_queries = get_sequences(all_queries, cache_dir, args)

    # retry failed accs
    no_seq_acc = []
    if len(failed_queries) != 0:
        logger.warning(
            "Parsing accessions which could not retrieve a seq for the first time\n"
            f"Handling {len(failed_queries)} failed batches"
        )
        # break up and query individually
        retrying_acc = {}  # {acc: # of tries}
        for batch in failed_queries:
            for acc in batch:
                retrying_acc[acc] = 0

        finished_retry = False

        acc_list = list(retrying_acc.keys())

        no_seq_acc = set()  # set of accessions for which no seq could be retrieved

        while finished_retry is False:
            acc_list = list(retrying_acc.keys())

            logger.warning(
                "Handling failed protein accessions\n"
                f"{len(acc_list)} protein accessions remaining"
            )

            for accession in tqdm(acc_list, desc="Handling failed accessions"):
                new_seq, failed_seq = get_sequences([[accession]], cache_dir, args, disable_prg_bar=True)

                if len(new_seq) == 0:
                    logger.warning(
                        f"Failed to retrieve sequence for {accession} on try no. {retrying_acc[accession]}"
                    )
                    retrying_acc[accession] += 1

                    if retrying_acc[accession] >= args.retries:
                        logger.warning(
                            f"Could not retrieve sequence for {accession}\n"
                            "Ran out of attempts"
                        )
                        del retrying_acc[accession]
                        no_seq_acc.add(accession)

            acc_list = list(retrying_acc.keys())

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


def get_seqs(batches, cache_path, invalid_ids_cache, args, disable_prg_bar=False):
    """Retrieve protein sequences from Entrez.

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
    failed_connections_batches = []  # nested list of batches during the querying of which the connection failed

    seq_records = []

    successful_accessions = []  # accessions for which seqs were retrieved

    for batch in tqdm(batches, desc="Batch querying NCBI.Entrez", disable=disable_prg_bar):
        epost_webenv, epost_query_key, success = post_accessions_to_entrez(batch)

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
            failed_connections_batches.append(batch)
            continue
        
        # else was a success, so proceed to retrieving protein seqs
        seq_records, success, temp_successful_accessions = fetch_ncbi_seqs(seq_records, epost_webenv, epost_query_key, gbk_acc_to_retrieve)
        
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
            failed_connections_batches.append(batch)
            continue

        SeqIO.write(seq_records, cache_path, "fasta")
        successful_accessions += temp_successful_accessions

    accs_still_to_fetch = [acc for acc in gbk_acc_to_retrieve if ((acc not in invalid_ids) or (acc not in successful_accessions))]

    return seq_records, batches_with_invalid_ids, failed_connections_batches, accs_still_to_fetch
        

def post_accessions_to_entrez(batch):
    """Post NCBI protein accessions to NCBI via Entrez, and capture error message if one returned.

    :param batch: list of genbank accessions
    :param entrez_function: Entrez function (epost or efetch)

    Return Entrez ePost web environment and query key.
    """
    logger = logging.getLogger(__name__)
    success, epost_result = None, None
    
    try:
        epost_result = Entrez.read(
            entrez_retry(
                args.retries,
                Entrez.epost,
                db="Protein",
                id=",".join(batch),
            ),
            validate=False,
        )

    except RuntimeError as err:
        if repr(err).startswith("RuntimeError('Some IDs have invalid value and were omitted.") or repr(err).startswith("RuntimeError('Empty ID list; Nothing to store')"):
            
            if len(batch) == 1:
                logger.warning(
                    f"Accession '{batch[0]}' not listed in NCBI. Querying NCBI returns the error:\n"
                    f"{repr(err)}\n"
                    f"Not retrieving seq for '{batch[0]}'"
                )
                success = "Invalid ID"
            
            else:
                logger.warning(
                    "Batch contains an accession no longer listed in NCBI\n"
                    f"Querying NCBI returns the error:\n{err}\n"
                    "Will identify invalid ID(s) later"
                )
                success = "Contains invalid ID"

        else:  # unrecognised Runtime error
            if len(batch) == 1:
                logger.warning(
                    f"Runtime error occured when querying NCBI by accession '{batch[0]}'\n"
                    f"Error returned:\n{err}\n"
                    "Interrupting error as recognition of an invalid ID. Therefore,\n"
                    f"not adding seq for '{batch[0]} to the local CAZyme db'"
                )
                success = "Invalid ID"

            else:
                logger.warning(
                    f"Runtime error occured when querying NCBI with batch of gbk accessions\n"
                    f"Error returned:\n{err}\n"
                    "Interrupting error as batch containing invalid ID.\n"
                    "Will identify invalid ID(s) later"
                )
                success = "Contains invalid ID"

    except IncompleteRead as err:
        logger.warning(
            "IncompleteRead error raised when querying NCBI:\n"
            f"{err}\n"
            "Will reattempt NCBI query later"
        )
        success = "Failed connection"

    except NotXMLError as err:
        logger.warning(
            "NotXMLError raised when querying NCBI:\n"
            f"{err}\n"
            "Will reattempt NCBI query later"
        )
        success = "Failed connection"

    except (TypeError, AttributeError) as err:  # if no record is returned from call to Entrez
        logger.warning(
            f"Error occurenced when batch quering NCBI\n"
            "Error retrieved:\n"
            f"{repr(err)}\n"
            "Will retry batch later"
        )
        success = "Failed connection"

    except Exception as err:  # if no record is returned from call to Entrez
        logger.warning(
            f"Error occurenced when batch quering NCBI\n"
            "Error retrieved:\n"
            f"{repr(err)}\n"
            "Will retry batch later"
        )
        success = "Failed connection"

    # retrieve the web environment and query key from the Entrez post
    try:
        epost_webenv = epost_result["WebEnv"]
        epost_query_key = epost_result["QueryKey"]
        success = "Complete"
    except (TypeError, AttributeError):
        epost_webenv, epost_query_key = None, None
        pass  # raised when could not epost failed

    return epost_webenv, epost_query_key, success


def fetch_ncbi_seqs(seq_records, epost_webenv, epost_query_key, acc_to_retrieve):
    """Retrieve protein sequences from NCBI from ePost of protein v.accs.

    :param seq_records: list of Bio.SeqRecords
    :param epost_websenv: Entrez ePost webenvironment
    :param epost_query_key: Entrez ePost query key
    :param acc_to_retrieve: set of NCBI protein version accessions to retrieve seqs for

    Return updated list of SeqRecords and string marking success/failure or seq retrieval.
    """
    logger = logging.getLogger(__name__)
    success, successful_accessions = None, []

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
                retrieved_accession = None

                # check if multiple items returned in ID
                try:
                    retrieved_accession = [_ for _ in record.id.split("|") if _.strip() in gbk_acc_to_retrieve][0]
                except IndexError:
                    # try re search for accession in string
                    try:
                        retrieved_accession = re.match(r"\D{3}\d+\.\d+", record.id).group()
                    except AttributeError:
                        try:
                            retrieved_accession = re.match(r"\D{2}_\d+\.\d+", record.id).group()
                        except AttributeError:
                            logger.warning(
                                f"Could not fetch protein acc from record id '{record.id}'.\n"
                                "Will search all target accessions against record id"
                            )
                            for acc in acc_to_retrieve:
                                if record.id.find(acc) != -1:
                                    retrieved_accession = acc
                
                if retrieved_accession is None:
                    logger.error(
                        "Could not retrieve a NCBI protein version accession matching\n"
                        f"an accession from the local database from the record id '{record.id}'\n"
                        "The sequence from this record will not be added to the db"
                    )
                    continue

                seq_records.append(record)
                successful_accessions.add(retrieved_accession)

    except RuntimeError as err:
        if repr(err).startswith("RuntimeError('Some IDs have invalid value and were omitted.") or repr(err).startswith("RuntimeError('Empty ID list; Nothing to store')"):
            
            if len(batch) == 1:
                logger.warning(
                    f"Accession '{batch[0]}' not listed in NCBI. Querying NCBI returns the error:\n"
                    f"{repr(err)}\n"
                    f"Not retrieving seq for '{batch[0]}'"
                )
                success = "Invalid ID"
            
            else:
                logger.warning(
                    "Batch contains an accession no longer listed in NCBI\n"
                    f"Querying NCBI returns the error:\n{err}\n"
                    "Will identify invalid ID(s) later"
                )
                success = "Contains invalid ID"

        else:  # unrecognised Runtime error
            if len(batch) == 1:
                logger.warning(
                    f"Runtime error occured when fetching record from NCBI for accession '{batch[0]}'\n"
                    f"Error returned:\n{err}\n"
                    "Interrupting error as recognition of an invalid ID. Therefore,\n"
                    f"not adding seq for '{batch[0]} to the local CAZyme db'"
                )
                success = "Invalid ID"

            else:
                logger.warning(
                    f"Runtime error occured when fetching records from NCBI\n"
                    f"Error returned:\n{err}\n"
                    "Interrupting error as batch containing invalid ID.\n"
                    "Will identify invalid ID(s) later"
                )
                success = "Contains invalid ID"

    except IncompleteRead as err:
        logger.warning(
            "IncompleteRead error raised when fetching record from NCBI:\n"
            f"{err}\n"
            "Will reattempt NCBI query later"
        )
        success = "Failed connection"

    except NotXMLError as err:
        logger.warning(
            "NotXMLError raised when fetching record(s) from NCBI:\n"
            f"{err}\n"
            "Will reattempt NCBI query later"
        )
        success = "Failed connection"

    except (TypeError, AttributeError) as err:  # if no record is returned from call to Entrez
        logger.warning(
            f"Error occurenced when fetching records from NCBI\n"
            "Error retrieved:\n"
            f"{repr(err)}\n"
            "Will retry batch later"
        )
        success = "Failed connection"

    except Exception as err:  # if no record is returned from call to Entrez
        logger.warning(
            f"Error occurenced when batch quering NCBI\n"
            "Error retrieved:\n"
            f"{repr(err)}\n"
            "Will retry batch later"
        )
        success = "Failed connection"

    return seq_records, success, success_accessions


def parse_failed_connections(
    failed_connections_batches,
    seq_records,
    accs_to_fetch,
    cache_path,
    invalid_ids_cache, 
    args,
):
    """Parse batches that suffered failed connections on the first attempt.

    :param failed_connections_batches: list of batches
    :param seq_records: list of Bio.SeqRecords
    :param accs_to_fetch: list of NCBI protein accs to retrieve seqs for
    :param cache_path: path to cache downloaded seqs
    :param invalid_ids_cache: path to cache invalid IDs
    :param args: CLI args parser

    Return 
    * updated list of seq_records
    * batches with invalid IDs
    * accs still to fetch seqs for
    """
    


if __name__ == "__main__":
    main()
