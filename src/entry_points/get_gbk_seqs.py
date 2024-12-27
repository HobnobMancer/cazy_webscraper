#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
# (c) James Hutton Institute 2024
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# ehobbs@ebi.ac.uk
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
"""Retrieve proteins sequences from GenBank and populate the local database"""


import logging

from argparse import Namespace

from Bio import Entrez
from saintBioutils.utilities.file_io import make_output_directory

from src import closing_message
from src.cache.ncbi import get_cache_seqs
from src.databases.ncbi.sequences import get_seqs_from_ncbi
from src.sql import sql_orm
from src.sql.interface.add_data.scrape_log import log_scrape_in_db
from src.sql.interface.add_data.add_ncbi_seqs import update_ncbi_seqs
from src.sql.interface.connect import connect_existing_db
from src.sql.interface.filter_data.protein import filter_to_db_acc
from src.sql.interface.get_data.get_selected_gbks import get_ncbi_prot_accessions
from src.utilities.parse_configuration import get_expansion_configuration
from src.utilities.sanity_checks.get_gbk_seqs import sanity_check_inputs


logger = logging.getLogger(__name__)


def main(args: Namespace, time_stamp: str, start_time):
    sanity_check_inputs(args)
    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)
    date_today = time_stamp.split("_")[0]
    Entrez.email = args.email

    if args.seq_update:
        logger.warning("Enabled updating sequences")
    if args.cache_dir:
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

    with sql_orm.Session(bind=connection) as session:
        retrieved_data = "GenBank protein sequences"
        log_scrape_in_db(
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

    if args.seq_dict or args.seq_file:
        seq_dict = get_cache_seqs(start_time, args)  # {genbank_acc: Bio.Seq}
        # only keep acc that are in the db
        acc_in_db = filter_to_db_acc(args.database, set(seq_dict.keys()))
        seq_dict = {k: v for k, v in seq_dict.items() if k in acc_in_db}
    else:
        seq_dict = {}

    if args.file_only:
        logger.warning("Only adding Seqs in JSON and/or FASTA file. Not retrieving seqs from NCBI")
    else:
        # get accession of records to retrieve sequences for
        if args.genbank_accessions:
            seq_acc_to_retrieve = get_acc_from_file(
                args,
                start_time
            )
        else:
            seq_acc_to_retrieve = get_ncbi_prot_accessions(
                class_filters,
                family_filters,
                kingdom_filters,
                taxonomy_filter_dict,
                ec_filters,
                args.database
            )

        if seq_acc_to_retrieve:
            new_seqs = get_seqs_from_ncbi(seq_acc_to_retrieve, cache_dir, args)
            seq_dict.update(new_seqs)

    if not seq_dict:
        logger.warning("No seqs to add to db")
        return("get_gbk_seqs")

    update_ncbi_seqs(seq_dict, args.database, args.seq_update)
    return("get_gbk_seqs")


def get_acc_from_file(
    args: Namespace,
    start_time: str
) -> set[str]:
    """Get accession from user file and add to seq_acc_to_retrieve"""
    logger.warning("Getting GenBank accessions from file: %s", args.genbank_accessions)
    seq_acc_to_retrieve = set()
    try:
        with open(args.genbank_accessions, "r") as fh:
            for line in fh:
                seq_acc_to_retrieve.add(line.strip())
    except FileNotFoundError:
        logging.error(
            "Could not find list of GenBank accessions at: %s\n"
            "Check the path is correct\n"
            "Terminating program", args.genbank_accessions
        )
        closing_message("Get GenBank seqs", start_time, args, early_term=True)

    seq_acc_to_retrieve = filter_to_db_acc(args.database, seq_acc_to_retrieve)

    return seq_acc_to_retrieve
