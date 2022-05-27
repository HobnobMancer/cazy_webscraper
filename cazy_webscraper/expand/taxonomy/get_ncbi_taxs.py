#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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
"""Produce dataframe listing NCBI species listed in a local CAZyme database, including their full lineage"""


import logging
import time

import pandas as pd

from datetime import datetime
from pathlib import Path
from typing import List, Optional, Type

from Bio import Entrez
from saintBioutils.utilities.logger import config_logger
from saintBioutils.utilities import file_io
from saintBioutils.genbank import entrez_retry
from tqdm import tqdm

from cazy_webscraper import cazy_scraper, closing_message
from cazy_webscraper.expand import get_chunks_list
from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks, get_table_dicts
from cazy_webscraper.sql.sql_interface.get_data.get_records import (
    get_user_genbank_sequences,
    get_user_uniprot_sequences
)
from cazy_webscraper.utilities.parsers.tax_ncbi_parser import build_parser
from cazy_webscraper.utilities.parse_configuration import get_expansion_configuration


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy."""
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # Program preparation
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args, logger_name=__name__)

    connection, logger_name, cache_dir = cazy_scraper.connect_existing_db(args, time_stamp, start_time)
    logger.info(f"Open connection to local cazyme database: {str(args.database)}")

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        logger.info("Building cache dir")
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        logger.info("Building cache dir")
        cache_dir = cache_dir / "ncbi_tax_retrieval"
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)

    # get config data
    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = get_expansion_configuration(args)

    gbk_dict = {}

    if args.genbank_accessions is not None:
        logger.warning(f"Getting GenBank accessions from file: {args.genbank_accessions}")
        with open(args.genbank_accessions, "r") as fh:
            lines = fh.read().splitlines()
        
        accessions = [line.strip() for line in lines]
        accessions = set(accessions)

        gbk_dict = get_selected_gbks.get_ids(accessions, connection)

    if args.uniprot_accessions is not None:
        logger.warning(f"Extracting protein sequences for UniProt accessions listed in {args.uniprot_accessions}")
        gbk_table_dict = get_table_dicts.get_gbk_table_dict(connection)
        uniprot_table_dict = get_table_dicts.get_uniprot_table_dict(connection)
        gbk_dict.update(get_user_uniprot_sequences(gbk_table_dict, uniprot_table_dict, args))

    else:
        gbk_dict.update(get_selected_gbks.get_genbank_accessions(
            class_filters,
            family_filters,
            taxonomy_filter_dict,
            kingdom_filters,
            ec_filters,
            connection,
        ))

    tax_ids, prot_id_dict = get_ncbi_tax_prot_ids(list(gbk_dict.keys()), cache_dir, args)


def get_ncbi_tax_prot_ids(protein_accessions, cache_dir, args):
    """Retrieve NCBI Taxonomy and Protein database IDs

    :param protein_accessions: list of NCBI protein version accessions
    :param cache_dir: Path, path to cache directory
    :param args: cmd-line args parser
    
    Return set of NCBI Tax ids and dict {ncbi prot id: prot acc}
    """
    logger = logging.getLogger(__name__)

    # break up list into batches
    batches = get_chunks_list(protein_accessions, args.batch_size)

    failed_batches = {}  # {batch: # of tries}

    tax_ids = set()

    protein_ncbi_ids = {}

    for batch in tqdm(batches, desc="Getting NCBI Tax IDs"):
        new_tax_ids, new_prot_ids, failed_batches = get_ncbi_ids(batch, args, failed_batches)
        
        tax_ids = tax_ids.union(new_tax_ids)
        protein_ncbi_ids.update(new_prot_ids)
    
    failed_individuals = {}  # prot_acc: tries

    if len(failed_batches.keys()) != 0:
        for batch in failed_batches:
            for prot in batch:
                new_tax_ids, new_prot_ids, failed_individuals = get_ncbi_ids([prot], args, failed_individuals)

                tax_ids = tax_ids.union(new_tax_ids)
                protein_ncbi_ids.update(new_prot_ids)
    
    failed_proteins = set()

    if len(list(failed_individuals.keys())) != 0:
        prots_to_retry = list(failed_individuals.keys())

        while len(list(failed_individuals.keys())) > 0:
            for prot in prots_to_retry:
                try:
                    failed_individuals[prot]
                except KeyError:
                    continue

                new_tax_ids, new_prot_ids, failed_individuals = get_ncbi_ids([prot], args, failed_individuals)

                tax_ids = tax_ids.union(new_tax_ids)
                protein_ncbi_ids.update(new_prot_ids)

                if failed_individuals[prot] >= args.retries:
                    del failed_individuals[prot]
                    failed_proteins.add(prot)
    
    if len(failed_proteins) != 0:
        with open((cache_dir / "failed_protein_accs.txt"), "a") as fh:
            for prot in failed_proteins:
                fh.write(f"{prot}\n")

    return tax_ids, protein_ncbi_ids


def get_ncbi_ids(prots, args, failed_batches):
    """Retrieve NCBI Taxonomy and Protein database IDs, for first time or retried attempts.

    :param protein_accessions: list of NCBI protein version accessions
    :param args: cmd-line args parser
    :param round: int, (1) first time, (2) second, (3) continued until reach retry limit
    :param failed_batches: dict {batch/prots: #ofTries}
    
    Return set of NCBI Tax ids and set of NCBI Prot ids, and failed_batches dict
    """
    logger = logging.getLogger(__name__)

    new_tax_ids, new_protein_ids = set(), {}

    # post protein accessions
    try:
        query_key, web_env = post_ids(prots, 'Protein', args)
    except RuntimeError:
        logger.warning(
            "Batch contained invalid protein version accession"
        )
        try:
            failed_batches[prots] += 1
        except KeyError:
            failed_batches[prots] = 1

        return new_tax_ids, new_protein_ids, failed_batches

    if query_key is None:
        try:
            failed_batches[prots] += 1
        except KeyError:
            failed_batches[prots] = 1

        return new_tax_ids, new_protein_ids, failed_batches
    
    # elink proteins to tax records
    new_tax_ids, new_protein_ids_list, success = link_prot_taxs(query_key, web_env, args)

    if success is False:
        logger.warning(
            f"Could not retrieve tax link for {','.join(prots)}."
        )

        try:
            failed_batches[prots] += 1
        except KeyError:
            failed_batches[prots] = 1

        return new_tax_ids, new_protein_ids, failed_batches

    for i in range(len(new_protein_ids)):
        prot_acc = prots[i]
        prot_id = new_protein_ids_list[i]
        new_protein_ids[prot_id] = prot_acc

    return new_tax_ids, new_protein_ids, failed_batches


def link_prot_taxs(query_key, web_env, args):
    """Get NCBI tax and protein IDs from elinking protein accs to tax records
    
    :param query_key: str, from Entrez epost
    :param web_env: str, from Entrez epost
    :param args: cmd-line args parser
    
    Return set of NCBI tax ids, and set of NCBI protein IDS, and bool TRUE if successful
    """
    logger = logging.getLogger(__name__)

    tax_ids, protein_ids = set(), set()

    try:
        with entrez_retry(
            args.retries,
            query_key=query_key,
            WebEnv=web_env,
            dbfrom="Protein",
            db="Taxonomy",
            linkname="protein_taxonomy",
        ) as handle:
            tax_links = Entrez.read(handle, validate=False)
    except RuntimeError:
        logger.warning("Batch included invalid NCBI protein ids")
        return tax_ids, protein_ids, False

    for result in tax_links:

        for prot_id in result['IdList']:
            protein_ids.add(prot_id)

        for item in result['LinkSetDb']:
            links = item['Link']
            for link in links:
                tax_ids.add(link['Id'])
    
    return tax_ids, list(protein_ids), True


def get_lineage_data(tax_ids, prot_id_dict, gbk_dict, args):
    """Retrieve lineage data from NCBI Taxonomy and associated tax ids with protein ids
    
    :param tax_ids: set of ncbi tax record ids
    :param prot_id_dict: dict {prot_id: prot_acc}
    :param gbk_dict: dict {protein acc: local db id}
    :param args: cmd-line args parser
    
    Return dict {tax_id: {linaege info, 'proteins' {local db protein ids}}
    """
    logger = logging.getLogger(__name__)

    tax_prot_dict = {}

    failed_ids = {}

    for tax_id in tqdm(tax_ids, desc="Retrieving lineages"):
        # get lineage information
        lineage_dict = get_lineage(tax_id, args)

        if lineage_dict is None:
            failed_ids[tax_id] = {'lineage': 1}
        
        # for ncbi db ids for proteins from local db that are linked to the tax record
        tax_prot_dict = get_tax_proteins(tax_id, prot_id_dict, gbk_dict, args)
        # {tax_id: {local db protein ids}}

        if tax_prot_dict is None:
            try:
                failed_ids[tax_id]['proteins'] = 1
            except KeyError:
                failed_ids[tax_id] = {'lineage': lineage_dict, 'proteins': 1}
                continue

        lineage_dict[tax_id]['proteins'] = tax_prot_dict[tax_id]

        tax_prot_dict[tax_id] = lineage_dict[tax_id]

    




def get_lineage(tax_id, args):
    """Retrieve lineage from NCBI taxonomy record, and add to lineage dict
    
    :param tax_id: str, ncbi tax db id
    :param args: cmd-line args parser

    Return dict of lineage data
    """
    logger = logging.getLogger(__name__)

    tax_dict = {}

    try:
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            db="Taxonomy",
            id=tax_id,
        ) as handle:
            tax_records = Entrez.read(handle, validate=False)
    
    except (TypeError, AttributeError) as err:
        logger.warning(f"Failed to fetch tax record from NCBI tax for id '{tax_id}'':\n{err}")
        return
    
    for record in tax_records:
        record_id = record['TaxId']

        # set lineage data to None
        kingdom, phylum, tax_class, order, family, genus, species, strain = None, None, None, None, None, None, None, None

        for i in record['LineageEx']:
            rank = i['Rank']

            if rank == 'superkingdom':
                kingdom = i['ScientificName']

            elif rank == 'phylum':
                phylum = i['ScientificName']

            elif rank == 'class':
                tax_class = i['ScientificName']

            elif rank == 'order':
                order = i['ScientificName']

            elif rank == 'family':
                family = i['ScientificName']

            elif rank == 'genus':
                genus = i['ScientificName']

            elif rank == 'species' or 'species group':
                species = i['ScientificName']

            elif rank == 'serotype' or 'strain':
                strain = i['ScientificName']

        scientific_name = record['ScientificName']

        if genus is not None and species is not None and strain is not None:
            strain = scientific_name.replace(f"{genus} {species}", "").strip()

        elif genus is not None and species is None:
            species = scientific_name.repace(genus, "").strip()

        tax_dict[record_id] = {record_id: {
                'kingdom': kingdom,
                'phylum': phylum,
                'class': tax_class,
                'order': order,
                'family': family,
                'genus': genus,
                'species': species,
                'strain': strain,
            }
        }

    return tax_dict


def get_tax_proteins(tax_id, prot_id_dict, gbk_dict, args):
    """Get the proteins linked to a tax id in NCBI, and link the tax id with the local db protein ids
    
    :param tax_id: str, NCBI tax db id
    :param prot_id_dict: dict {protein ncbi id: prot acc}
    :param gbk_dict: dict, {prot acc: local db id}
    :param args: cmd-line args parser
    
    Return dict {tax_id: {local db protein ids}}, or None if fails
    """
    logger = logging.getLogger(__name__)

    try:
        with entrez_retry(
            args.retries,
            Entrez.elink,
            id=tax_id,
            db="Protein",
            dbfrom="Taxonomy",
            linkname="taxonomy_protein",
        ) as handle:
            tax_links = Entrez.read(handle, validate=False)

    except (AttributeError, TypeError, RuntimeError) as err:
        logger.warning(f"Failed to link NCBI tax id to NCBI Protein db for tax id {tax_id}\n{err}")
        return

    tax_prot_dict = {tax_id: set()}
    
    for result in tax_links:
        for item in result['LinkSetDb']:
            links = item['Link']
            for link in links:
                linked_prot_id = link['Id']

                # check if from the local database
                try:
                    prot_ver_acc = prot_id_dict[linked_prot_id]
                except KeyError:
                    continue

                prot_local_db_id = gbk_dict[prot_ver_acc]

                tax_prot_dict[tax_id].add([prot_local_db_id])

    return tax_prot_dict


if __name__ == "__main__":
    main()
