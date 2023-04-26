#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
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
from socket import timeout
from typing import List, Optional
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from Bio import Entrez
from tqdm import tqdm
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger

from cazy_webscraper import closing_message, connect_existing_db
from cazy_webscraper.expand import get_chunks_list
from cazy_webscraper.ncbi.genomes import (
    get_nuccore_ids,
    get_assembly_ids,
    get_assembly_data,
)
from cazy_webscraper.sql.sql_interface.add_data.add_genome_data import add_assembly_data
from cazy_webscraper.sql.sql_interface.get_data.get_records import (
    get_user_genbank_sequences,
    get_user_uniprot_sequences
)
from cazy_webscraper.sql.sql_interface.get_data.get_selected_gbks import get_genbank_accessions
from cazy_webscraper.sql.sql_interface.get_data.get_table_dicts import (
    get_gbk_table_dict,
    get_uniprot_table_dict,
)
from cazy_webscraper.sql.sql_interface.get_data.get_assemblies import get_no_assembly_proteins
from cazy_webscraper.utilities import parse_configuration
from cazy_webscraper.utilities.parsers.get_genomes_parser import build_parser


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

    with sql_orm.Session(bind=connection) as session:
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            'NCBI Assembly',
            'Genomic assemly accessions',
            session,
            args,
        )

    logger.info(f"Retrieving Genbank records from the local db:\n{str(args.database)}")

    gbk_dict = get_gbks_of_interest(
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
        connection,
        args,
    )

    genbank_accessions = list(gbk_dict.keys())

    if len(genbank_accessions) == 0:
        logger.warning(
            "No records matching the given criteria found in the local CAZyme database:\n"
            f"{args.database}"
        )
        closing_message("get_genomic_accessions", start_time, args)
        sys.exit(1)

    logger.info(f"Retrieved {len(genbank_accessions)} from the local db")

    # separate GenBank and RefSeq accessions
    gbk_accessions = set()
    refseq_accessions = set()

    for accession in tqdm(genbank_accessions, desc="Separting GenBank and RefSeq accessions"):
        if accession[2] == '_':
            refseq_accessions.add(accession)
        else:
            gbk_accessions.add(accession)

    logger.info(f"Retrieved {len(gbk_accessions)} GenBank accessions from the local db")
    logger.info(f"Retrieved {len(refseq_accessions)} RefSeq accessions from the local db")

    genome_dict = {}
    assembly_dict = {}

    # get data for genbank accessions
    if len(gbk_accessions) != 0:
        new_assembly_dict, new_genome_dict = get_ncbi_assembly_data(
            list(gbk_accessions),
            cache_dir,
            args,
        )
        assembly_dict.update(new_assembly_dict)
        genome_dict.update(new_genome_dict)

    # get data for refseq accessions
    if len(refseq_accessions) != 0:
        new_assembly_dict, new_genome_dict = get_ncbi_assembly_data(
            list(refseq_accessions),
            cache_dir,
            args,
        )
        assembly_dict.update(new_assembly_dict)
        genome_dict.update(new_genome_dict)

    # insert data in to the db
    add_assembly_data(assembly_dict, genome_dict, gbk_dict, connection, args)

    closing_message("Get genomic data", start_time, args)


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

    if args.update:
        gbk_dict = get_no_assembly_proteins(gbk_dict, connection)

    else:
        # only retrieve data for proteins with no assembly data in the local db
        return gbk_dict
        
    return gbk_dict


def get_ncbi_assembly_data(sequence_accessions, cache_dir, args, refseq=False):
    """Param retrieve assembly data for list of proteins
    
    :param sequence_accessions: list of protien accessions
    :param cache_dir: path to cache directory
    :param args: cmd-line args parser
    :param refseq: bool, are the protein accessions from the NCBI RefSeq db?
    
    Return dict of {ass_name: {protein_accessions} and genome_dict, listing assembly meta data
    """
    logger = logging.getLogger(__name__)

    genome_dict, failed_batches, parsed_assembly_ids = get_genomic_assembly_data(sequence_accessions, args)

    if len(failed_batches['proteins']) != 0 or len(failed_batches['nuccores']) != 0 or len(failed_batches['assemblies']) != 0:
        genome_dict.update(retry_failed_batches(
            genome_dict,
            failed_batches,
            parsed_assembly_ids,
            cache_dir,
            args,
        ))

    assembly_dict = {}

    no_urls = cache_dir / "no_assembly_urls.txt"

    failed_feature_tables = {}

    # download assembly feature tables and link assemblies to protein accessions
    for assembly_name in tqdm(genome_dict, desc="Linking proteins to genomic assemblies"):
        feature_table_url = get_feature_table_url(assembly_name, genome_dict[assembly_name], refseq, no_urls)
        if feature_table_url is None:
            continue

        file_name = feature_table_url.split('/')[-1]
        out_file_path = cache_dir / file_name

        successful = download_feature_table(assembly_name, feature_table_url, out_file_path, args)

        if successful is False:
            failed_feature_tables[feature_table_url] = 1
            continue

        feature_table = pd.read_csv(out_file_path, sep="\t")

        product_accessions = set(feature_table['product_accession'])

        for product_accession in tqdm(product_accessions, f"Parsing {assembly_name} feature table"):
            if product_accession in sequence_accessions:
                try:
                    assembly_dict[assembly_name].add(product_accession)
                except KeyError:
                    assembly_dict[assembly_name] = {product_accession}

    if len(list(failed_feature_tables.keys())) != 0:
        logger.info("Retrying failed downloads of feaure tables")
        done_urls = set()

        while len(done_urls) < len(list(failed_feature_tables.keys())):
            for feature_table_url in failed_feature_tables:
                if feature_table_url in done_urls:
                    continue

                # retrieve assembly name from url
                assembly_name = feature_table_url.split("/")[-2].split("_")[-1]

                file_name = feature_table_url.split('/')[-1]
                out_file_path = cache_dir / file_name

                successful = download_feature_table(assembly_name, feature_table_url, out_file_path, args)

                if successful is False:
                    failed_feature_tables[feature_table_url] += 1
                    logger.warning(
                        f"Failed to download feature table from {feature_table_url}\n"
                        f"on the {failed_feature_tables[feature_table_url]} try"
                    )
                    if failed_feature_tables[feature_table_url] > args.retries:
                        done_urls.add(feature_table_url)
                        with open(no_urls, 'a') as fh:
                            fh.write(f"{assembly_name}\n")
                    continue
            
                feature_table = pd.read_csv(out_file_path, sep="\t")

                product_accessions = set(feature_table['product_accession'])

                for product_accession in tqdm(product_accessions, f"Parsing {assembly_name} feature table"):
                    if product_accession in sequence_accessions:
                        try:
                            assembly_dict[assembly_name].add(product_accession)
                        except KeyError:
                            assembly_dict[assembly_name] = {product_accession}

                done_urls.add(feature_table_url)

    return assembly_dict, genome_dict


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

        logger.info(f"Fetching data for {len(nuccore_ids_to_fetch)} nuccore IDs")

        # get assembly IDs linked to the nuccore IDs
        assembly_ids, failed_batches = get_assembly_ids(nuccore_ids_to_fetch, failed_batches, args)

        assembly_ids_to_fetch = [_ for _ in assembly_ids if _ not in parsed_assembly_ids]

        if len(assembly_ids_to_fetch) == 0:
            continue

        logger.info(f"Fetching data for {len(assembly_ids_to_fetch)} assembly IDs")

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

    return genome_dict, failed_batches, parsed_assembly_ids


def retry_failed_batches(genome_dict, failed_batches, parsed_assembly_ids, cache_dir, args):
    """Retry failed batches of protein, nuccore and assembly ids.

    :param genome_dict: dict of assembly meta data (assembly name, ids, version accs, urls)
    :param failed_batches: dict, dict of failed protein, nuccore and assembly ids
    :param parsed_assembly_ids: set of ncbi assembly IDs for which data has already been retrieved
    :param cache_dir: path to cache dir
    :param args: cmd-line args parser
    
    Return genome_dict (dict of assembly meta data from NCBI)
    """
    logger = logging.getLogger(__name__)

    if len(failed_batches['proteins']) != 0:

        failed_proteins_cache = cache_dir / "failed_protein_accessions.txt"

        logger.warning("Retrying failed batches of protein version accessions")

        # Parse protein version accessions individually to find those that are causing the previous issues
        for batch in tqdm(failed_batches['proteins'], desc="Retrying failed protein version acc batches"):
            for prot_ver_acc in batch:
                new_genome_dict, new_failed_batches, new_parsed_assembly_ids = get_genomic_assembly_data([prot_ver_acc], args)

                if len(new_failed_batches['proteins']) != 0:
                    logger.error(
                        "Could not retrieve genomic assembly accessions for the "
                        f"for protein {prot_ver_acc}"
                    )
                    with open(failed_proteins_cache, "a") as fh:
                        fh.write(f"{prot_ver_acc}\n")
                    continue
                
                genome_dict.update(new_genome_dict)

                for assembly_id in new_parsed_assembly_ids:
                    parsed_assembly_ids.add(assembly_id)

    if len(failed_batches['nuccores']) != 0:
        failed_nuccore_cache = cache_dir / "failed_nuccore_ids.txt"
        logger.warning("Retrying failed batches of NCBI Nuccore db IDs")
        for batch in tqdm(failed_batches['nuccores'], desc="Retrying failed nuccore db IDs"):
            # parse nuccore IDs individually to find those that caused issues
            for nuccore_id in batch:
                new_failed_batches = {'nuccores': [], 'assemblies': []}

                # get assembly ids
                assembly_ids, new_failed_batches = get_assembly_ids([nuccore_id], failed_batches, args, retry=True)

                if len(new_failed_batches['nuccores']) != 0:
                    logger.error(
                        "Could not retrieve genomic assembly accessions for the "
                        f"for nuccore ID {nuccore_id}"
                    )
                    with open(failed_nuccore_cache, "a") as fh:
                        fh.write(f"{nuccore_id}\n")
                    continue

                assembly_ids_to_fetch = [_ for _ in assembly_ids if _ not in parsed_assembly_ids]

                if len(assembly_ids_to_fetch) == 0:
                    continue

            new_assembly_data, new_failed_batches, new_parsed_assembly_ids = get_assembly_data(
                assembly_ids_to_fetch,
                failed_batches,
                parsed_assembly_ids,
                args,
            )

            if len(new_failed_batches['nuccores']) != 0 or len(new_failed_batches['assemblies']) != 0:
                logger.error(
                    "Could not retrieve genomic assembly accessions for the "
                    f"for nuccore ID {nuccore_id}"
                )
                continue

            genome_dict.update(new_assembly_data)

            for assembly_id in new_parsed_assembly_ids:
                parsed_assembly_ids.add(assembly_id)

    if len(failed_batches['assemblies']) != 0:
        failed_assembly_cache = cache_dir / "failed_assembly_ids.txt"
        logger.warning("Retrying failed batches of NCBI Assembly db IDs")
        for batch in failed_batches['assemblies']:
            # parse assembly IDs individually to find those that are causing issues
            for assembly_id in batch:
                new_failed_batches = {'nuccores': [], 'assemblies': []}

                new_assembly_data, new_failed_batches, new_parsed_assembly_ids = get_assembly_data(
                    assembly_ids_to_fetch,
                    failed_batches,
                    parsed_assembly_ids,
                    args,
                )

                if len(new_failed_batches['nuccores']) != 0 or len(new_failed_batches['assemblies']) != 0:
                    logger.error(
                        "Could not retrieve genomic assembly accessions for the "
                        f"for assembly ID {nuccore_id}"
                    )
                    with open(failed_assembly_cache, "a") as fh:
                        fh.write(f"{assembly_id}\n")
                    continue

                genome_dict.update(new_assembly_data)

                for assembly_id in new_parsed_assembly_ids:
                    parsed_assembly_ids.add(assembly_id)
    
    return genome_dict


def get_feature_table_url(assembly_name, assembly_data, refseq, cache_path):
    """Retrieve the url for downloading the feature table from the ncbi ftp server
    
    :param assembly_name: str, name of genomic assembly
    :param assembly_data: dict, contains ids, urls and other meta data for the assembly
    :param refseq: bool, is the assembly associated with refseq sequence accessions?
    :param cache_path: path to write out cache
    
    Return url (str), or none if none retrieved
    """
    logger = logging.getLogger(__name__)

    if refseq:
        feature_table_url = assembly_data['refseq_url']
        if len(feature_table_url) == 0:
            logger.warning(
                f"Did not retrieve RefSeq url for assembly {assembly_name}\n"
                "Will use the Gbk feature table instead"
            )
            feature_table_url = assembly_data['gbk_url']
            if len(feature_table_url) == 0:
                logger.warning(
                    f"Did not retrieve RefSeq or Gbk url for assembly {assembly_name}\n"
                    "Cannot retrieve feature table and linke proteins to this assembly"
                )
                with open(cache_path, "a") as fh:
                    fh.write(f"{assembly_name}\n")
    else:
        feature_table_url = assembly_data['gbk_url']
        if len(feature_table_url) == 0:
            logger.warning(
                f"Did not retrieve Gbk url for assembly {assembly_name}\n"
                "Will use the Gbk feature table instead"
            )
            feature_table_url = assembly_data['refseq_url']
            if len(feature_table_url) == 0:
                logger.warning(
                    f"Did not retrieve RefSeq or Gbk url for assembly {assembly_name}\n"
                    "Cannot retrieve feature table and linke proteins to this assembly"
                )
                with open(cache_path, "a") as fh:
                    fh.write(f"{assembly_name}\n")

    if len(feature_table_url) == 0:
        return
    else:
        return feature_table_url


def download_feature_table(assembly_name, feature_table_url, out_file_path, args):
    """Download feature table from the ftp server
    
    :param assembly_name: str, name of assembly in NCBI
    :param feature_table_url: str, url to feature table to be downloaded
    :param out_file_path: path to write out downloaded file to disk
    :param args: cmd-line args parser
    
    Return bool, marking if successful or not
    """
    logger = logging.getLogger(__name__)

    try:
        response = urlopen(feature_table_url, timeout=args.timeout)
    except (HTTPError, URLError, timeout) as e:
        logger.error(
            (
                f"Failed to download the feature table for {assembly_name}:\n"
                f"{feature_table_url}\n"
                "See traceback below\n"
            ), exc_info=1,
        )
        return False

    file_size = int(response.info().get("Content-length"))
    bsize = 1_048_576
    logger.info("Opened URL and parsed metadata")

    try:
        with open(out_file_path, "wb") as out_handle:
            while True:
                buffer = response.read(bsize)
                if not buffer:
                    break
                out_handle.write(buffer)
    except IOError:
        logger.error(f"Download feature table failed for {assembly_name}", exc_info=1)
        return False
    
    # check file size is correct
    if out_file_path.stat().st_size != file_size:
        logger.error(f"Don't not compelte download for {assembly_name}")
        return False
    
    return True


if __name__ == "__main__":
    main()
