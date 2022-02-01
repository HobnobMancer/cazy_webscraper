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
"""Retrieve proteins sequences from GenBank and populate the local database"""


import json
import logging
import re

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio import Entrez, SeqIO
from saintBioutils.utilities import file_io
from saintBioutils.genbank import entrez_retry
from saintBioutils.utilities.logger import config_logger

from cazy_webscraper import cazy_scraper, closing_message
from cazy_webscraper.sql import sql_orm, sql_interface
from cazy_webscraper.sql.sql_interface import get_selected_gbks, add_genbank_data
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

    connection, logger_name, cache_dir = cazy_scraper.connect_existing_db(args, time_stamp, start_time)

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "genbank_data_retrieval"
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)

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
            seq_dict = json.load(fh)

    else:
        seq_dict = get_sequences(genbank_accessions, args)  # {gbk_accession: seq}

        # cache the retrieved sequences
        cache_path = cache_dir / f"genbank_seqs_{time_stamp}.json"
        with open(cache_path, "w") as fh:
            json.dump(seq_dict, fh)

    add_genbank_data.add_gbk_seqs_to_db(seq_dict, date_today, connection, args)

    closing_message("get_genbank_sequences", start_time, args)


def get_sequences(genbank_accessions, args):
    """Retrieve protein sequences from Entrez.

    :param genbank_accessions: list, GenBank accessions
    :param args: cmb-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    seq_dict = {}  # {gbk_accession: seq}
    protein_records = set()

    epost_webenv, epost_query_key = bulk_query_ncbi(genbank_accessions, args)

    success_accessions = set()  # accessions for which seqs were retrieved

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
            index, success = 0, False
            
            while (index < len(temp_accession)) and (success is False):
                acc = temp_accession[index]
                try:
                    re.match(
                        (
                            r"(\D{3}\d{5,7}\.\d+)|"
                            r"(\D\d(\D|\d){3}\d)|"
                            r"((\D\d(\D|\d){3}\d)\.\d)|"
                            r"(\D\d(\D|\d){3}\d\D(\D|\d){2}\d)|"
                            r"((\D\d(\D|\d){3}\d\D(\D|\d){2}\d)\.\d)|"
                            r"(\D\D_\d{9}\.\d)"
                        ),
                        acc,
                    ).group()
                    seq_accession = acc
                    success = True
                except AttributeError:
                    index += 1
                    pass

            if success is False:  # if could not retrieve GenBank accession from the record
                logger.error(
                    "Could not retrieve a GenBank protein accession from the record with the id:\n"
                    f"{record.id}\n"
                    "The sequence from this record will not be added to the db"    
                )
                continue

            if seq_accession not in genbank_accessions:  # if the retrieved Gbk acc does not match an acc in the db
                logger.warning(
                    f"Retrieved the accession {temp_accession} from the record with the id:\n"
                    f"{record.id},\n"
                    "Not adding the protein seq from this record to the db"
                )
                continue
                
            seq_dict[seq_accession] = record.seq

            success_accessions.add(seq_accession)

            protein_records.add(record)

    SeqIO.write(protein_records, 'pl_seq.fasta', "fasta")

    #         # remove the accession from the list
    #         try:
    #             genbank_accessions.remove(temp_accession)
    #         except ValueError:
    #             logger.warning(
    #                 f"Tried to remove {temp_accession} from list of accessions, "
    #                 "but it was not in the list of accessions.\n"
    #                 "The returned accession and the one present in CAZy do not match."
    #             )
    
    # if len(genbank_accessions) != 0:
    #     # logger.warning("Protein sequences for the following CAZymes were not retrieved:")
    #     # for acc in genbank_accessions:
    #     #     logger.warning(f"GenBank accession: {acc}")

    #     cache_path = cache_dir / "no_seq_retrieved.txt"
    #     with open(cache_path, "a") as fh:
    #         for acc in genbank_accessions:
    #             fh.write(f"{acc}\n")

    return seq_dict


def bulk_query_ncbi(accessions, args):
    """Bulk query NCBI and retrieve webenvironment and history tags
    :param accessions: list of GenBank protein accessions
    :param args: cmd-line args parser
    Return webenv and query key
    """
    # perform batch query of Entrez
    accessions_string = ",".join(accessions)

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


if __name__ == "__main__":
    main()
