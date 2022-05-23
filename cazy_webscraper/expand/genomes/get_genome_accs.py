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


import logging
import sys

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
        logger = logging.getLogger(__name__)
        config_logger(args, logger_name=__name__)

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

    logger.info(f"Retrieving Genbank records from the local db:\n{str(args.database)}")

    # # genbank_accessions = get_gbks_of_interest(
    # #     class_filters,
    # #     family_filters,
    # #     kingdom_filters,
    # #     taxonomy_filter_dict,
    # #     ec_filters,
    # #     connection,
    # #     args,
    # # )
    # # if len(genbank_accessions) == 0:
    # #     logger.warning(f"No records matching the given criteria found in the local CAZyme database:\n{args.database}")
    # #     closing_message("get_genomic_accessions", start_time, args)
    # #     sys.exit(1)

    with open('genbank_accessions.txt', 'r') as fh:
        genbank_accessions = fh.read().splitlines()

    logger.info(f"Retrieved {len(genbank_accessions)} from the local db")

    genome_dict, failed_batches = get_genomic_assembly_data(genbank_accessions, args)

    if len(failed_batches) != 0:
        logger.warning("Retrying failed batches")
        for batch in tqdm(failed_batches, desc="Retrying failed batches"):

            for protein in batch:
                new_assembly_data, failed_batches, parsed_nuccore_ids, parsed_assembly_ids = get_genome_data(
                    [protein],
                    failed_batches,
                    parsed_nuccore_ids,
                    parsed_assembly_ids,
                    args,
                )

                genome_dict.update(new_assembly_data)
    # genome_dict = {assembly_name: {gbk and refseq accessions and uids, and urls to download feature tables}}

    logger.info(f"Identfied {len(genome_dict.keys())} assembly names")

    print(genome_dict)
    sys.exit(1)

    # download assemblies and associate with protein accessions
    genome_protein_relationships = get_relationships(genome_dict, genbank_accessions, cache_dir, args)

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


def get_genomic_assembly_data(genbank_accessions, args):
    """Retrieve the meta data for genomic assemblies in NCBI linked to a list of protein version accessions
    
    :param genbank_accessions: list of Protein GenBank version accessions
    :param args: cmd-line args parser
    
    Return dict of assembly names, ids, version accessions and urls to download feature tables,
        and dict of failed protein, nuccore and assembly IDs -- ids for which no data was retrieved
        from NCBI
    """
    logger = logging.getLogger(__name__)

    failed_batches = {'proteins': [], 'nuccores': [], 'assemblies': []}  # store lists of IDs that cause issues

    genome_dict = {}  # used for storing retrieved genomic accessions

    parsed_nuccore_ids = set()
    parsed_assembly_ids = set()

    # break up long list into smaller batches
    batches = get_chunks_list(genbank_accessions, args.batch_size)

    logger.info("Starting retrieval of genomic accessions")

    for protein_accs in tqdm(batches, desc="Batch quering NCBI"):
        # get nuccore IDs for the set of protein accessions
        nuccore_ids, failed_batches = get_nuccore_ids(protein_accs, failed_batches, args)

        if len(nuccore_ids) == 0:
            failed_batches['proteins'].append(protein_accs)

        # check nuccore IDs have laredy been parsed because multiple proteins are linked to the same nuccore record
        nuccore_ids_to_fetch = [_ for _ in nuccore_ids if _ not in parsed_nuccore_ids]

        if len(nuccore_ids_to_fetch) == 0:
            continue

        # get assembly IDs linked to the nuccore IDs
        aseembly_ids, failed_batches = get_assembly_ids(nuccore_ids_to_fetch)

        assembly_ids_to_fetch = [_ for _ in aseembly_ids if _ not in parsed_assembly_ids]

        if len(assembly_ids_to_fetch == 0):
            continue

        new_assembly_data, failed_batches, parsed_assembly_ids = get_assembly_data(
            assembly_ids_to_fetch,
            failed_batches,
            parsed_assembly_ids,
            args,
        )

        genome_dict.update(new_assembly_data)

        # update parsed nuccore and assembly ids
        for nuccore_id in nuccore_ids:
            parsed_nuccore_ids.add(nuccore_id)

    return genbank_accessions, failed_batches


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
            f"Failed to post IDs to Entrez {database} db:\n{err}"
        )
        return None, None

    query_key = posted_records['QueryKey']
    web_env = posted_records['WebEnv']

    return query_key, web_env


def get_nuccore_ids(batch, failed_batches, args, retry=False):
    """Retrieve the IDs of nuccore records linkded to a list of protein version accessions

    :param batch: list of GenBank protein accessions
    :param failed_batches: dict listing batches for which data was not retrieved
    :param args: cmd-line args parser
    :param retry: bool, is this a retry of previously failed retrieval of genome data?

    Return set() of nuccore db IDs and failed_batches
    """
    logger = logging.getLogger(__name__)

    nuccore_ids = set()

    logger.info("Posting protein accessions to NCBI")

    try:
        query_key, web_env = post_ids(batch, "Protein", args)

        if query_key is None:
            failed_batches['proteins'].append(batch)
            return nuccore_ids, failed_batches

    except RuntimeError:
        if retry:
            logger.warning(f"{batch[0]} is not listed in NCBI")
        else:
            logger.warning("Batch contains invalid NCBI Protein db accessions")
            failed_batches.append(batch)

        return nuccore_ids, failed_batches

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
        return nuccore_ids, failed_batches

    for record in tqdm(nuccore_records, desc="Get Nuccore record IDs"):
        if len(record['LinkSetDb']) != 0:
            for nuc_id in record['LinkSetDb'][0]['Link']:
                nuccore_ids.add(nuc_id['Id'])
    
    return nuccore_ids, failed_batches


def get_assembly_ids(nuccore_ids, failed_batches, args, retry=False):
    """Retrieve the IDs of assembly records linkded to a list of nuccore record Ids

    :param nuccore_ids: list of uccore ncbi db IDs
    :param failed_batches: dict listing batches for which data was not retrieved
    :param args: cmd-line args parser
    :param retry: bool, is this a retry of previously failed retrieval of genome data?

    Return set() of assebly db IDs and failed_batches
    """
    logger = logging.getLogger(__name__)

    assembly_ids = set()

    logger.info("Posting nuccore IDs")
    # post nuccore IDs
    try:
        query_key, web_env = post_ids(nuccore_ids, "Nuccore", args)

        if query_key is None:
            failed_batches['nuccores'].append(nuccore_ids)
            return assembly_ids, failed_batches

    except RuntimeError:
        if retry:
            logger.warning(f"Nuccore ID '{nuccore_ids[0]}' is not listed in NCBI")
        else:
            logger.warning("Batch contains invalid NCBI Nuccore db IDs")
            failed_batches.append(nuccore_ids)

        return nuccore_ids, failed_batches

    logger.info("Getting linked assembly IDs")
    try:
        logger.info("Try")
        with entrez_retry(
            args.retries,
            Entrez.elink,
            query_key=query_key,
            WebEnv=web_env,
            dbfrom="Nuccore",
            db="Assembly",
            linkname="nuccore_assembly",
        ) as handle:
            linked_records = Entrez.read(handle, validate=False)
    except (TypeError, AttributeError, RuntimeError) as err:
        logger.warning(f"Failed to link nuccore records to assembly records:\n{err}")
        return assembly_ids, failed_batches

    for record in tqdm(linked_records, desc="Getting assembly ids"):
        for index in range(len(record['LinkSetDb'][0]['Link'])):
            assembly_id = record['LinkSetDb'][0]['Link'][index]['Id']
            assembly_ids.add(assembly_id)

    return assembly_ids, failed_batches


def get_assembly_data(assembly_ids, failed_batches, parsed_assembly_ids, args, retry=False):
    """Retrieve the data for assemblies represented by their NCBI Assembly DB ID

    :param assembly_ids: list of assembly ncbi db IDs
    :param failed_batches: dict listing batches for which data was not retrieved
    :param args: cmd-line args parser
    :param parsed_assembly_ids: set of ncbi assembly db IDs that have already been parsed
    :param retry: bool, is this a retry of previously failed retrieval of genome data?

    Return dict of assembly meta data and failed batches dict
    """
    logger = logging.getLogger(__name__)

    genome_dict = {}

    # post assembly IDs
    try:
        query_key, web_env = post_ids(assembly_ids, "Assembly", args)

        if query_key is None:
            failed_batches.append(assembly_ids)
            return genome_dict, failed_batches
    except RuntimeError:
        if retry:
            logger.warning(f"Data for Assembly ID '{assembly_ids[0]}' could not be retrieved from NCBI")
        else:
            logger.warning("Batch contains invalid NCBI Assembly IDs")
            failed_batches['assemblies'].append(assembly_ids)
        return genome_dict, failed_batches

    try:
        with entrez_retry(
            args.retries,
            Entrez.esummary,
            db="Assembly",
            query_key=query_key,
            WebEnv=web_env,
            rettype="docsum",
            retmode="xml",
        ) as record_handle:
            assembly_records = Entrez.read(record_handle, validate=False)
    except (TypeError, AttributeError, RuntimeError) as err:
        logger.warning(f"Failed to retrieve Assembly records:\n{err}")

    for genome in assembly_records['DocumentSummarySet']['DocumentSummary']:
        assembly_name = genome['AssemblyName']

        try:
            genome_dict[assembly_name]
            continue
        except KeyError:
            pass

        gbk_uid = genome['GbUid']
        ref_uid = genome['RsUid']

        try:
            gbk_acc = genome['Synonym']['Genbank']
        except KeyError:
            gbk_acc = None
        
        try:
            refseq_acc = genome['Synonym']['RefSeq']
        except KeyError:
            refseq_acc = None

        try:
            gbk_ftp_path = genome['FtpPath_GenBank']
            file_stem = gbk_ftp_path.split("/")[-1]
            gbk_url = f"{gbk_ftp_path}/{file_stem}_feature_table.txt.gz"
        except KeyError:
            gbk_url = None
        
        try:
            ref_ftp_path = genome['FtpPath_RefSeq']
            file_stem = ref_ftp_path.split("/")[-1]
            ref_url = f"{ref_ftp_path}/{file_stem}_feature_table.txt.gz"
        except KeyError:
            ref_url = None

        genome_dict[assembly_name] = {
            "gbk_acc": gbk_acc,
            "refseq_acc": refseq_acc,
            "gbk_uid": gbk_uid,
            "refseq_uid": ref_uid,
            "gbk_url": gbk_url,
            "refseq_url": ref_url,
        }

        parsed_assembly_ids.add(gbk_uid)
        parsed_assembly_ids.add(ref_uid)

    return genome_dict, failed_batches, parsed_assembly_ids








def get_relationships(genome_dict, protein_accessions, cache_dir, args):
    """Download genomes and track which proteins come from which genomes.
    
    :param genome_dict: 
    :param protein_accessions: list of protein accessions of interest
    :param cache_dir: path to cache directory
    :param args: cmd-line args parser
    
    Return dict {assembly_name: {protein accessions}}
    """
    relationship_dict = {}  # {assembly name: {protein accs}}

    for assembly_name in tqdm(genome_dict, desc="Associating genomes with protein accessions"):
        gbk_feature_table_url = genome_dict[assembly_name]['gbk_url']

        # download file

        # open

        # get protein id column

        # get all ids that are in proteins of interest

        # add key:value to relationships dict

    # if any proteins of interest left
    
    # download refseq feature tables as well

    return relationship_dict


if __name__ == "__main__":
    main()
