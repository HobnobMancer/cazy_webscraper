#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2025
# (c) University of Strathclyde 2025
# (c) James Hutton Institute 2025
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

import logging

from argparse import Namespace

from saintBioutils.utilities.file_io import make_output_directory

from src.databases.uniprot import get_uniprot_data
from src.sql import sql_orm
from src.sql.interface.add_data.scrape_log import log_scrape_in_db
from src.sql.interface.get_data.get_proteins import get_ncbi_prot_accessions
from src.sql.interface.connect import connect_existing_db
from src.utilities.parse_configuration import get_expansion_configuration
from src.utilities.parse_configuration.accession_files import get_acc_from_file


logger = logging.getLogger(__name__)

def main(args: Namespace, time_stamp: str, start_time):
    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)

    cache_dir = args.cache_dir if args.cache_dir else (cache_dir / "uniprot_data_retrieval")
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
        retrieved_data = "UniProt protein data"
        log_scrape_in_db(
            time_stamp,
            config_dict,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            'UniProt',
            retrieved_data,
            session,
            args,
        )

    if args.genbank_accessions:
        seq_acc_to_retrieve = get_acc_from_file(
            args.genbank_accessions,
            args.database,
        )
    elif args.update or args.go or args.pdb or args.ec:
        # retreive everything as we can update these fields
        seq_acc_to_retrieve = get_ncbi_prot_accessions(
            class_filters,
            family_filters,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            args.database
        )
    else:
        # only retrieve accessions without uniprot data
        seq_acc_to_retrieve = get_ncbi_prot_accessions(
            class_filters,
            family_filters,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            args.database,
            additional_filter="P.uniprot_id IS NULL"
        )

    if seq_acc_to_retrieve:
        logger.warning("Retrieving UniProt data for %d accessions", len(seq_acc_to_retrieve))
        get_uniprot_data(seq_acc_to_retrieve, cache_dir, time_stamp, args)

    return "get_uniprot_data"
