#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
"""Retrieve the accessions for proteins of interest, and store accessions in the local db"""


from asyncio.log import logger
from cmath import log
import logging
import sys
from numpy import False_

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio import Entrez
from tqdm import tqdm
from saintBioutils.genbank import entrez_retry
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger

from cazy_webscraper import closing_message, connect_existing_db
from cazy_webscraper.sql.sql_interface import get_selected_gbks, get_table_dicts
from cazy_webscraper.expand import get_chunks_list
from cazy_webscraper.utilities import parse_configuration
from cazy_webscraper.utilities.parsers.get_genomes_parser import build_parser
from cazy_webscraper.sql.sql_interface.get_records import (
    get_user_genbank_sequences,
    get_user_uniprot_sequences
)


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
        logger = logging.getLogger(__package__)
        config_logger(args)

    Entrez.email = args.email

    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)

    logger.info(f"Connected to local db: {args.database}")

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "genomes"
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    
    logger.info(f"Using cache dir: {cache_dir}")

    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = parse_configuration.get_expansion_configuration(args)

    logger.info("Retrieving Genbank records from the local db")

    genbank_accessions = get_gbks_of_interest(
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
        connection,
        args,
    )
    if len(genbank_accessions) == 0:
        logger.warning(f"No records matching the given criteria found in the local CAZyme database:\n{args.database}")
        closing_message("get_genomic_accessions", start_time, args)
        sys.exit(1)
    
    logger.info(f"Retrieved {len(genbank_accessions)} from the local db")

    # get assembly accessions, URLs and IDs from NCBI

    # download assemblies and associate with protein accessions
    
    # add data to the local db

    closing_message("extract_sequences", start_time, args)


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
    
    Return list of genbank protein accessions
    """
    logger = logging.getLogger(__name__)

    gbk_dict = {}  # {gbk_acc: db_id}

    if args.genbank_accessions is not None or args.uniprot_accessions is not None:
        gbk_table_dict = get_table_dicts.get_gbk_table_dict(connection)

        if args.genbank_accessions is not None:
            logger.info(f"Retrieving PDB structures for GenBank accessions listed in {args.genbank_accessions}")
            gbk_dict.update(get_user_genbank_sequences(gbk_table_dict, args))

        if args.uniprot_accessions is not None:
            logger.info(f"Extracting protein sequences for UniProt accessions listed in {args.uniprot_accessions}")
            uniprot_table_dict = get_table_dicts.get_uniprot_table_dict(connection)
            gbk_dict.update(get_user_uniprot_sequences(gbk_table_dict, uniprot_table_dict, args))

        return list(gbk_dict.keys())
    
    gbk_dict = get_selected_gbks.get_genbank_accessions(
        class_filters,
        family_filters,
        taxonomy_filter_dict,
        kingdom_filters,
        ec_filters,
        connection,
    )

    return list(gbk_dict.keys())


def get_genomic_accessions(protein_accessions, args):
    """Retrieve genomic accessions for a set of proteins from NCBI.
    
    :param protein_accessions: list, GenBank protein verion accessions
    :param args: cmd-line args parser
    
    Return dict {assembly_name: {'genbank': str, 'refseq': str, 'gbk_url': str, 'ref_url': str}}
    """
    logger = logging.getLogger(__name__)

    logger.info("Starting retrieval of genomic accessions")

    failed_batches = []  # store lists that causes issues

    genome_dict = {}  # used for storing retrieved genomic accessions

    # break up long list into smaller batches
    batches = get_chunks_list(protein_accessions, args.batch_size)

    for batch in tqdm(batches, desc="Batch quering NCBI"):
        # post portein accessions

        try:
            query_key, web_env = post_ids(batch, "Protein", args)

            if query_key is None:
                failed_batches.append(batch)
                continue
        except RuntimeError:
            logger.warning("Batch contains invalid NCBI Protein db accessions")
            failed_batches.append(batch)
            continue

        # link protein records to nuccore records
        nuccore_ids = get_linked_nuccore_ids(query_key, web_env, args)

        if nuccore_ids is None:
            failed_batches.append(batch)
            continue

        # post nuccore IDs
        try:
            query_key, web_env = post_ids(nuccore_ids, "Nuccore", args)

            if query_key is None:
                failed_batches.append(batch)
                continue
        except RuntimeError:
            logger.warning("Batch contains invalid NCBI Nuccore IDs")
            failed_batches.append(batch)
            continue

        # fetch nuccore records
        nuccore_id_protein_acc, records_with_unknown_id = link_nuccore_ids_to_gbk_accs(
            query_key,
            web_env,
            nuccore_ids,
            batch,
            args,
        )

        # replace nuccore version acc if could not retrieve from record
        if len(records_with_unknown_id) != 0:
            nuccore_id_protein_acc = get_nuccore_acc(nuccore_id_protein_acc, records_with_unknown_id, args)


def post_ids(ids, database, args):
    """Post protein IDs to Entrez
    
    :param ids: list, GenBank protein accession numbers
    :param database: str, Name of database from which IDs are sourced
    :param args: cmd-line args parser
    
    Return None (x2) if fails
    Else return query_key and web_env
    """
    logger = logging.getLogger(__name__)

    try:
        with entrez_retry(
            args.retries,
            Entrez.epost,
            db=database,
            id=",".join(ids),
        ) as handle:
            posted_records = Entrez.read(handle, validate=False)

    # if no record is returned from call to Entrez
    except (TypeError, AttributeError) as err:
        logger.warning(
            f"Failed to post IDs to Entrez {database} db"
        )
        return None, None

    query_key = posted_records['QueryKey']
    web_env = posted_records['WebEnv']

    return query_key, web_env


def get_linked_nuccore_ids(query_key, web_env, args):
    """Retrieve the IDs of nuccore records linkded to posted protein version accessions
    
    :param query_key: str, query key from Entrez epost
    :param web_env: str, web env from Entrez epost
    :param args: cmd-line args parser
    
    Return set of nuccore db IDs
    """
    logger = logging.getLogger(__name__)

    try:
        with entrez_retry(
            args.retries,
            Entrez.elink,
            query_key=query_key,
            WebEnv=web_env,
            dbfrom="Protein",
            db="Nuccore",
            linkname="protein_nuccore",
        ) as handle:
            nuccore_records = Entrez.read(handle, validate=False)

    except (TypeError, AttributeError, RuntimeError) as error:
        logger.warning(
            f"Entrez failed to link Protein records to nuccore records numbers\n",
            error
        )
        return
    
    nuccore_ids = set()

    for record in tqdm(nuccore_records, desc="Get Nuccore record IDs"):
        if len(record['LinkSetDb']) != 0:
            for nuc_id in record['LinkSetDb'][0]['Link']:
                nuccore_ids.add(nuc_id['Id'])
    
    if len(nuccore_records) == 0:
        return
    
    return nuccore_records


def link_nuccore_ids_to_gbk_accs(query_key, web_env, nuccore_ids, protein_accs, args):
    """Retrieve the IDs of nuccore records and associated with Protein accession versions
    
    :param query_key: str, query key from Entrez epost
    :param web_env: str, web env from Entrez epost
    :param nuccore_ids: list of nuccore db IDs
    :param protein_accs: list of GenBank Protein db version accessions
    :param args: cmd-line args parser
    
    Return dict {nuccore_id: protein_acc}
    """
    logger = logging.getLogger(__name__)

    try:
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            db="Nucleotide",
            query_key=query_key,
            WebEnv=web_env,
            retmode="xml",
        ) as handle:
            nuccore_records = Entrez.read(handle, validate=False)
    except (TypeError, AttributeError, RuntimeError) as err:
        logger.warning(f"Failed to fetch nuccore records:\n{err}")
    
    records_with_unknown_id = set()
    nuccore_id_protein_acc = {}  # {nuccore id: protein acc}

    # retrieve the protein sequence accession
    for record in nuccore_records:
        nuccore_acc = record['GBSeq_accession-version']
        
        # find id in the record
        record_id = None
        for nuccore_id in nuccore_ids:
            if str(record).find(nuccore_id) != -1:
                record_id = nuccore_id
        
        # get the protein accession
        protein_acc = None
        
        for feature in record['GBSeq_feature-table']:
            if feature['GBFeature_key'] == 'CDS':
                for qual in feature['GBFeature_quals']:
                    if qual['GBQualifier_name'] == 'protein_id':
                        if qual['GBQualifier_value'] in protein_accs:
                            protein_acc = qual['GBQualifier_value']
        
        if protein_acc is not None:
            if record_id is None:
                records_with_unknown_id.add(nuccore_acc)
                nuccore_id_protein_acc[nuccore_acc] = protein_acc
            else:
                nuccore_id_protein_acc[record_id] = protein_acc
    
    return nuccore_id_protein_acc, records_with_unknown_id


def get_nuccore_acc(nuccore_id_protein_acc, records_with_unknown_id, args):
    """Replace nuccore version acc with nuccore ID retrieved from record summary
    
    :param nuccore_id_protein_acc: dict, {nuccore_acc or ID: protein_acc}
    :param records_with_unknown_ids: xxx
    :param args: cmd-line args parser
    
    Return dict {nuccore_id: protein_acc}, and bool to represent if process completed successfully
    """
    logger = logging.getLogger(__name__)

    # post accessions
    try:
        query_key, web_env = post_ids(records_with_unknown_id, "Nuccore", args)
    except RuntimeError as err:
        logger.warning(f"Failed to post nuccore accessions:\n{err}")
        return nuccore_id_protein_acc, False
    
    if query_key is None:
        return nuccore_id_protein_acc, False
    
    # get nuccore esummaries
    try:
        with entrez_retry(
            args.retries,
            Entrez.esummary,
            db="Nucleotide",
            query_key=query_key,
            WebEnv=web_env,
            retmode="xml",
        ) as handle:
            nuccore_summaries = Entrez.read(handle, validate=False)
    except (TypeError, AttributeError, RuntimeError) as err:
        logger.warning(f"Failed to retrieve nuccore records\n:{err}")
        return nuccore_id_protein_acc, False
    
    for record in tqdm(nuccore_summaries, desc="Retrieving nuccore IDs"):
        nuccore_id = record['Id']
        nuccore_acc = record['AccessionVersion']
        # replace nuccore acc with nuccore id as key in dict
        try:
            nuccore_id_protein_acc[nuccore_id] = nuccore_id_protein_acc.pop(nuccore_acc)
        except KeyError:
            pass
    
    # check if failed to replace any accessions with IDs
    failed_replacements = set()
    for nuccore_id in nuccore_id_protein_acc:
        try:
            int(nuccore_id)
        except ValueError:
            failed_replacements.add(nuccore_id)
    
    if len(failed_replacements) != 0:
        logger.warning(
            "Failed to retrieve nuccore IDs for the following version accs:\n"
            f"{failed_replacements}\n"
            "Will no retrieve genomic accessions for these nuccore accessions"
        )
        for nuccore_acc in failed_replacements:
            del nuccore_id_protein_acc[nuccore_acc]
    
    return nuccore_id_protein_acc
