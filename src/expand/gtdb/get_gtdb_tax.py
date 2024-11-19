#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) Jame Hutton Institute 2022
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
#
# Bio.PDB reference:
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
# implemented in Python. Bioinformatics 19: 2308–2310
"""Retrieve taxonomic classifications from GTDB"""


import logging
import pandas as pd
import sys

from datetime import datetime
from typing import List, Optional

from requests.exceptions import ConnectionError, MissingSchema
from socket import timeout
from urllib3.exceptions import HTTPError, RequestError
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import mechanicalsoup

from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from cazy_webscraper import closing_message, connect_existing_db
from cazy_webscraper.expand.gtdb import get_gtdb_data
from cazy_webscraper.utilities.parse_configuration import get_expansion_configuration
from cazy_webscraper.sql import sql_orm
from cazy_webscraper.sql import sql_interface
from cazy_webscraper.sql.sql_interface.get_data.get_table_dicts import (
    get_gbk_table_dict,
    get_uniprot_table_dict,
)
from cazy_webscraper.sql.sql_interface.get_data.get_records import (
    get_user_genbank_sequences,
    get_user_uniprot_sequences
)
from cazy_webscraper.sql.sql_interface.get_data.get_selected_gbks import get_genbank_accessions
from cazy_webscraper.sql.sql_interface.get_data.get_assemblies import get_genomes
from cazy_webscraper.sql.sql_interface.add_data.add_gtdb_tax import (
    add_gtdb_taxs,
    add_genome_gtdb_relations,
)
from cazy_webscraper.utilities.parsers.get_gtdb_parser import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # parse cmd-line arguments
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        args = build_parser(argv).parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args, logger_name=__name__)

    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)

    logger.info(f"Connected to local db: {args.database}")

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "gtdb"
        make_output_directory(cache_dir, args.force, args.nodelete_cache)

    logger.info(f"Using cache dir: {cache_dir}")

    # get configuration data
    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = get_expansion_configuration(args)

    with sql_orm.Session(bind=connection) as session:
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            'Genome Taxonomy DataBase (GTDB)',
            'GTDB taxonomc lineages',
            session,
            args,
        )

    # get the GenBank verion accessions and local db IDs of proteins matching the config criteria
    gbk_dict = get_gbks_of_interest(
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
        connection,
        args,
    )

    # genome_dict = {local db genome id: {
    #     'gbk_genome': str-version acc,
    #     'refseq_genome': str
    # }}
    # selected_genomes = set of genome version accessions
    genome_dict, selected_genomes = get_genomes(gbk_dict, args, connection)

    if len(selected_genomes) == 0:
        logger.error(
            "Retrieved no genomes matching specified criteria from the local database.\n"
            "Terminating program"
        )
        sys.exit(1)

    logger.info(f"Retrieving GTDB tax classification for {len(selected_genomes)} genomes")

    archaea, bacteria = False, False
    if 'archaea' in args.taxs:
        archaea = True
    if 'bacteria' in args.taxs:
        bacteria = True

    archaea_file, bacteria_file = get_gtdb_data(args, cache_dir, arch=archaea, bact=bacteria)

    # get the lineage data from the GTDB data files
    genome_lineage_dict = {}  # {genome_version_accession: {'lineage': lineage. 'release': release-str}
    if 'archaea' in args.taxs:
        genome_lineage_dict.update(get_lineage_data(archaea_file, selected_genomes))

    if 'bacteria' in args.taxs:
        genome_lineage_dict.update(get_lineage_data(bacteria_file, selected_genomes))

    if len(list(genome_lineage_dict.keys())) == 0:
        logger.error(
            "Retrieved no lineage data for genomes retrieved from the local db.\n"
            "Terminating program"
        )
        sys.exit(1)

    # add data to the local CAZyme db
    add_gtdb_taxs(genome_lineage_dict, connection)
    add_genome_gtdb_relations(genome_lineage_dict, args, connection)

    closing_message("Get-GTDB-taxs", start_time, args)


def get_gbks_of_interest(
    class_filters,
    family_filters,
    kingdom_filters,
    taxonomy_filter_dict,
    ec_filters,
    connection,
    args,
):
    """Retrieve the Gbk protein accessions of proteins matching the user criteria

    :param class_filters: set of CAZy classes of interest
    :param family_filters: set of CAZy families and subfamilies of interest
    :param kingdom_filters: set of taxonomic kingdoms of interest
    :param taxonomy_filter_dict: dict of genera, species and strains of interst
    :param ec_filters: set of EC numbers of interest
    :param connection: connection to local SQL db
    :param args: cmd-line args parser

    Return dict {gbk_acc: db_id}
    """
    logger = logging.getLogger(__name__)

    gbk_dict = {}  # {gbk_acc: db_id}

    if args.genbank_accessions is not None or args.uniprot_accessions is not None:
        gbk_table_dict = get_gbk_table_dict(connection)

        if args.genbank_accessions is not None:
            logger.info(f"Retrieving PDB structures for GenBank accessions listed in {args.genbank_accessions}")
            gbk_dict.update(get_user_genbank_sequences(gbk_table_dict, args))

        if args.uniprot_accessions is not None:
            logger.info(f"Extracting protein sequences for UniProt accessions listed in {args.uniprot_accessions}")
            uniprot_table_dict = get_uniprot_table_dict(connection)
            gbk_dict.update(get_user_uniprot_sequences(gbk_table_dict, uniprot_table_dict, args))

    else:
        gbk_dict = get_genbank_accessions(
            class_filters,
            family_filters,
            taxonomy_filter_dict,
            kingdom_filters,
            ec_filters,
            connection,
        )

    return gbk_dict


def get_lineage_data(gtdb_file_path, selected_genomes):
    """Iterate through GTDB datafile, extracting info for genomes of interest

    :param gtdb_file_path: Path, GTDB datafile, csv
    :param selected_genomes: dict, {local db id: {refseq: str, gbk: str}}

    Return dict {genome_version_accession: {'lineage': lineage. 'release': release-str}}
    """
    release = gtdb_file_path.name.split("-")[-1].replace('.gz', '')
    genome_lineage_dict = {}
    i = 0
    for line in tqdm(
        pd.read_csv(gtdb_file_path, sep='\t', chunksize=1, names=['genome', 'lineage']),
        desc=f"Parsing GTDB file {gtdb_file_path.name}"
    ):
        genome = line['genome'][i].replace("RS_", "").replace("GB_", "")
        if genome in selected_genomes:
            genome_lineage_dict[genome] = {
                'lineage': line['lineage'][i],
                'release': release,
            }
        i += 1

    return genome_lineage_dict
