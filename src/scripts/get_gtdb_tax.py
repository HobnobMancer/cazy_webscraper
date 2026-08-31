#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# Author:
# Emma E. M. Hobbs

# Contact
# ehobbbs@ebi.ac.uk

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
"""Get genome genome data from NCBI Assembly database"""

import gzip
import logging
import re
import urllib.error
import requests
import time
import sqlite3

"""Coordinate retrieving taxonomic classifications from GTDB."""


import logging

from argparse import Namespace

from saintBioutils.utilities.file_io import make_output_directory

from src.databases.gtdb import get_gtdb_taxonomies
from src.sql import sql_orm
from src.sql.interface.add_data.add_gtdb_tax import persist_gtdb_data, rm_temp_gtdb_table
from src.sql.interface.add_data.scrape_log import log_scrape_in_db
from src.sql.interface.connect import connect_existing_db
from src.sql.interface.get_data.get_genomes import get_selected_genomes
from src.sql.interface.get_data.get_proteins import get_ncbi_acc_for_uniprot_acc
from src.utilities.parse_configuration import get_expansion_configuration
from src.utilities.parse_configuration.accession_files import get_acc_from_file


logger = logging.getLogger(__name__)


def main(args: Namespace, time_stamp: str, start_time):
    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)

    cache_dir = args.cache_dir if args.cache_dir else (cache_dir / "gtdb_retrieval")
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
        retrieved_data = "GTDB taxonomic classifications"
        log_scrape_in_db(
            time_stamp,
            config_dict,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            'GTDB',
            retrieved_data,
            session,
            args,
        )

    protein_accessions = None
    if args.genbank_accessions or args.uniprot_accessions:
        protein_accessions = set()
        if args.genbank_accessions:
            protein_accessions.update(
                get_acc_from_file(args.genbank_accessions, args.database)
            )
        if args.uniprot_accessions:
            uniprot_accessions = get_acc_from_file(
                args.uniprot_accessions, args.database, table="UniProt"
            )
            protein_accessions.update(
                get_ncbi_acc_for_uniprot_acc(uniprot_accessions, args.database)
            )

        if not protein_accessions:
            logger.warning("No proteins in the local db matched the accessions provided")
            return "get_gtdb_taxs"

        protein_accessions = sorted(protein_accessions)

    accession_to_genome_id = get_selected_genomes(
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
        args.database,
        protein_accessions=protein_accessions,
        # without --update only genomes that have no lineage yet are worth looking up
        only_missing_gtdb=not args.update,
    )

    if not accession_to_genome_id:
        logger.warning(
            "No genomes matched the criteria provided.\n"
            "GTDB classifications are retrieved for the genomes of proteins in the local db, so "
            "run 'cazy_webscraper get_ncbi_genomes' first if the database has no genomic data yet."
        )
        return "get_gtdb_taxs"

    genome_count = len(set(accession_to_genome_id.values()))
    logger.warning("Retrieving GTDB classifications for %s genomes", genome_count)

    get_gtdb_taxonomies(accession_to_genome_id, cache_dir, args)

    persist_gtdb_data(args)

    rm_temp_gtdb_table(args.database)

    return "get_gtdb_taxs"
