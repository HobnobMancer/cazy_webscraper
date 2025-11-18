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

from src.cache.ncbi import get_cache_seqs
from src.databases.ncbi.sequences import get_seqs_from_ncbi
from src.sql import sql_orm
from src.sql.interface.add_data.scrape_log import log_scrape_in_db
from src.sql.interface.add_data.add_ncbi_seqs import update_ncbi_seqs
from src.sql.interface.connect import connect_existing_db
from src.sql.interface.filter_data.protein import filter_to_db_acc
from src.sql.interface.get_data.get_proteins import get_ncbi_prot_accessions
from src.utilities.parse_configuration import get_expansion_configuration
from src.utilities.parse_configuration.accession_files import get_acc_from_file
from src.utilities.sanity_checks.get_ncbi_seqs import sanity_check_inputs


logger = logging.getLogger(__name__)


def main(args: Namespace, time_stamp: str, start_time):
    sanity_check_inputs(args)
    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)
    Entrez.email = args.email

    if args.update:
        logger.warning("Enabled updating sequences")

    cache_dir = args.cache_dir if args.cache_dir else (cache_dir / "ncbi_seq_retrieval")
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

    # get accession of records to retrieve sequences for
    if args.genbank_accessions:
        seq_acc_to_retrieve = get_acc_from_file(
            args.genbank_accessions,
            args.database,
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
        logger.warning("Retrieving %s sequences from NCBI", len(seq_acc_to_retrieve))
        get_seqs_from_ncbi(seq_acc_to_retrieve, cache_dir, time_stamp, args)

    return "get_ncbi_seqs"
