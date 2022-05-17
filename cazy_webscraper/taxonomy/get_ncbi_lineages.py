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
from cazy_webscraper.sql.sql_interface.get_data.get_taxonomies import get_taxonomies
from cazy_webscraper.utilities.parsers.tax_ncbi import build_parser
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
        config_logger(args)

    connection, logger_name, cache_dir = cazy_scraper.connect_existing_db(args, time_stamp, start_time)
    logger.info("Open connection to local cazyme database:", args.database)

    if str(args.output.parent) != '.':
        file_io.make_output_directory(args.output_dir, args.force, args.nodelete)

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        logger.info("Building cache dir")
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        logger.info("Building cache dir")
        cache_dir = cache_dir / "uniprot_data_retrieval"
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

    organisms = get_taxonomies(
        class_filters,
        family_filters,
        taxonomy_filter_dict,
        kingdom_filters,
        ec_filters,
        connection,
        args,
    )
    logger.info(f"Identified {len(organisms)} matching the provided criteria")

    # convert list of organisms into dict {genus: {species: {strain}}}
    genus_dict = build_genus_dict(organisms)
    genera = list(genus_dict).keys()
    logger.info(f"Retrievied {len(genera)} matching the provided criteria")

    Entrez.email = args.email

    lineage_dict = build_lineage_dict(genera, args)

    lineage_df = build_lineage_df(lineage_dict, genus_dict)
    logger.info(f"Compiled {len(lineage_df)} lineages")

    lineage_df.to_csv(args.output)

    closing_message("NCBI taxonomy summary", start_time, args)


def get_gbk_accessions(
    class_filters,
    family_filters,
    taxonomy_filter_dict,
    kingdom_filters,
    ec_filters,
    connection,
    args,
):
    """Retrieve the genbank accessions of proteins matching the user's specified criteria
    
    :param connection: open connection to local sql database
    :param args: cmd-line args parser
    
    Return list of genbank_accessions
    """
    logger = logging.getLogger(__name__)

    # retrieve genbank accessions of the proteins of interest

    if args.genbank_accessions is not None or args.uniprot_accessions is not None:
        gbk_dict = {}  # {gbk_acc: gbk_id}

        # retrieve GenBank accessions of proteins defined by user
        gbk_table_dict = get_table_dicts.get_gbk_table_dict(connection)
        # {genbank_accession: 'taxa_id': str, 'gbk_id': int}

        if args.genbank_accessions is not None:
            logger.warning(f"Retrieving PDB structures for GenBank accessions listed in {args.genbank_accessions}")
            gbk_dict.update(get_user_genbank_sequences(gbk_table_dict, args))

        if args.uniprot_accessions is not None:
            logger.warning(f"Extracting protein sequences for UniProt accessions listed in {args.uniprot_accessions}")
            uniprot_table_dict = get_table_dicts.get_uniprot_table_dict(connection)
            gbk_dict.update(get_user_uniprot_sequences(gbk_table_dict, uniprot_table_dict, args))

    else:
        gbk_dict = get_selected_gbks.get_genbank_accessions(
            class_filters,
            family_filters,
            taxonomy_filter_dict,
            kingdom_filters,
            ec_filters,
            connection,
        )

    return list(gbk_dict.keys())


def build_genus_dict(organisms):
    """Build a dict of {genus: {species: {strains}}} from list of organisms
    
    :param organisms: list, scientific names of organisms from a local cazyme database
    
    Return dict
    """
    logger = logging.getLogger(__name__)

    logger.info("Building dict of organisms of interest")

    db_tax_dict = {}  # {genus: {species: {strain}}}

    for organism in tqdm(organisms, desc="Compiling tax dict"):
        genus = organism.split(" ")[0]

        species_strain = organism.split(" ")[1]
        
        if species_strain.split(" ")[0] == 'sp.':
            species = species_strain
            strain = ""
        else:
            species = species_strain.split(" ")[0]
            strain = species_strain.replace(species, "").strip()
        
        try:
            db_tax_dict[genus]
            
            try:
                db_tax_dict[genus][species].add(strain)
            
            except KeyError:
                db_tax_dict[genus][species] = {strain}
        
        except KeyError:
            db_tax_dict[genus] = {species: {strain}}
    
    return db_tax_dict


def build_lineage_dict(genera, cache_dir, args):
    """Retrieve full lineage for each organism from the NCBI Taxonomy database
    
    :param genera: set of genera from the local CAZyme database
    :param cache_dir: path to cache dir
    :param args: cmd-line args parser
    
    Return dataframe of lineages, one unique organism per row.
    """
    logger = logging.getLogger(__name__)

    failed_batches = set()

    batches = get_chunks_list(genera, args.batch_size)

    lineage_dict = {}

    for batch in tqdm(batches, desc="Batch retrieving lineages from NCBI"):
        # search NCBI Taxonomy to retrieve tax record ids
        tax_id_dict, failed_genera = get_tax_record_ids(batch, args)
        
        if len(list(tax_id_dict.keys())) == 0:
            failed_batches.add(batch)
            continue
        
        if len(failed_genera) != 0:
            failed_batches.add(list(failed_genera))
        
        # post tax IDS
        query_key, web_env = post_to_entrez(batch, args)

        if query_key is None:
            failed_batches.add(batch)
            continue
            
        # fetch lineage from tax records
        lineage_dict, success = get_lineage(genera, lineage_dict, query_key, web_env, args)

        if success is False:
            failed_batches.add(batch)
    
    # handled failed batches
    if len(failed_batches) != 0:
        lineage_dict = parse_failed_batches(batches, args, lineage_dict, cache_dir)

    return lineage_dict
        

def get_tax_record_ids(genera, args):
    """Retrieve the NCBI Taxonomy db record IDs.
    
    :param genera: list of genera of interest
    :param query_key: str, query key from Entrez.epost
    :param web_env: str, web_env from Entrez.epost
    :param args: cmd-line args parser
    
    Return dict {scientific_name: record-id}
    """
    logger = logging.getLogger(__name__)

    tax_id_dict = {}  # {id: genus}

    failed_genera = set()

    # Entrez does not support posting scientific names to NCBI.Taxonomy
    for genus in genera:
        try:
            with entrez_retry(
                args.retries,
                Entrez.esearch,
                db="Taxonomy",
                term=genus,
            ) as handle:
                tax_record = Entrez.read(handle, validate=False)

        except (TypeError, AttributeError) as err:
            logger.warning(f"Could not retrieve record from NCBI Tax db\n{err}")
            continue

        term = tax_record['TranslationStack'][0]['Term'].split("[")[0]
        tax_id = tax_record['IdList'][0]

        tax_id_dict[tax_id] = term

        time.sleep(0.5)  # to prevent bombarding Entrez

    return tax_id_dict, failed_genera


def post_to_entrez(data, args):
    """Post data to NCBI entrez for batch query.
    
    :param data: list, data to be posted
    :param args: cmd-line args parser
    
    Return query key and web env, or None x2 if fails
    """
    logger = logging.getLogger(__name__)

    try:
        with entrez_retry(
            args.retries,
            Entrez.epost,
            db="Taxonomy",
            id=",".join(data),
        ) as handle:
            posted_record = Entrez.read(handle, validate=False)
    
    # if not recrd is returned
    except (TypeError, AttributeError) as err:
        logger.warning(f"Failed to post data to Entrez:\n{err}")
        return None, None
    
    query_key = posted_record['QueryKey']
    web_env = posted_record['WebEnv']

    return query_key, web_env


def get_lineage(genera, lineage_dict, query_key, web_env, args):
    """Retrieve lineage from NCBI taxonomy record, and add to lineage dict
    
    :param genera: str, list of genera
    :param lineage_dict: {superkingdom: {phylum: {class: {order: {genus}}}}}
    :param query_key: str, query key from Entrez.epost
    :param web_env: str, web_env from Entrez.epost
    :param args: cmd-line args parser

    Return lineage_dict and boolean to reflect successful or failed retrieval of records from NCBI
    """
    logger = logging.getLogger(__name__)

    try:
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            query_key=query_key,
            WebEnv=web_env,
        ) as handle:
            tax_records = Entrez.read(handle, validate=False)
    
    except (TypeError, AttributeError) as err:
        logger.warning("Failed to fetch tax records from NCBI")
        return lineage_dict, False
    
    for record in tqdm(tax_records, desc="Extracing lineages from NCBI Tax records"):
        # set lineage data to None
        superkingdom, phylum, tax_class, order, family = None, None, None, None, None

        # collect the lineage data
        genus = record['ScientificName']

        if genus not in genera:
            continue

        for i in record['LineageEx']:
            rank = i['Rank']

            if rank == 'superkingdom':
                superkingdom = i['ScientificName']

            elif rank == 'phylum':
                phylum = i['ScientificName']

            elif rank == 'class':
                tax_class = i['ScientificName']

            elif rank == 'order':
                order = i['ScientificName']

            elif rank == 'family':
                family = i['ScientificName']

        # add lineage to lineage dict
        try:
            lineage_dict[superkingdom]

            try:
                lineage_dict[superkingdom][phylum]

                try:
                    lineage_dict[superkingdom][phylum][tax_class]

                    try:
                        lineage_dict[superkingdom][phylum][tax_class][order]

                        try:
                            lineage_dict[superkingdom][phylum][tax_class][order][family].add(genus)

                        except KeyError:
                            lineage_dict[superkingdom][phylum][tax_class][order][family] = {genus}

                    except KeyError:
                        lineage_dict[superkingdom][phylum][tax_class][order] = {family: {genus}}

                except KeyError:
                    lineage_dict[superkingdom][phylum][tax_class] = {order: {family: {genus}}}

            except KeyError:
                lineage_dict[superkingdom][phylum] = {family: {order: {family: {genus}}}}

        except KeyError:
            lineage_dict[superkingdom] = {phylum: {family: {order: {family: {genus}}}}}
    
    return lineage_dict, True


def build_lineage_df(lineage_dict, genus_dict):
    """Combine two tax dicts into a single dataframe
    
    :param lineage_dict: {superkingdom: {phylum: {class: {order: {family: {genus}}}}}}
    :param genus_dict: {genus: {species: {strain}}}
    
    Return pandas df
    """
    df_data = []

    for kingdom in lineage_dict:
        kingdom_data = lineage_dict[kingdom]

        for phylum in kingdom_data:
            phylum_data = kingdom_data[phylum]

            for tax_class in phylum_data:
                class_data = phylum_data[tax_class]

                for order in class_data:
                    order_data = class_data[order]

                    for family in order_data:
                        family_data = order_data[family]

                        for genus in family_data:
                            genus_data = genus_dict[genus]

                            for species in genus_data:
                                if len(genus_data[species]) == 0:
                                    data = [kingdom, phylum, tax_class, order, genus, species, None]
                                    df_data.append(data)

                                else:
                                    species_data = genus_data[species]
                                    for strain in species_data:
                                        data = [kingdom, phylum, tax_class, order, genus, species, strain]
                                        df_data.append(data)
    
    df = pd.DataFrame(df_data, ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', ' Genus', 'Species', 'Strain'])

    return df


def parse_failed_batches(batches, args, lineage_dict, cache_dir):
    """Rerun failed batches.
    
    :param bataches: set of batches, each batch is a list of genera
    :param args: cmd-line args parser
    :param lineage_dict: 
    :param cache_dir: path to cache directory
    
    Return lineage dict
    """
    logger = logging.getLogger(__name__)

    failed_log = cache_dir / "failed_organisms.txt"
    failed_organisms = set()

    for batch in tqdm(batches, desc="Rerunning failed batches"):
        for genus in batch:
            tax_id_dict, failed_genera = get_tax_record_ids([genus], args)

            if len(failed_genera) != 0:
                logger.warning(f"Could not retrieve tax data from NCBI for {genus}")
                failed_organisms.add(genus)
                continue

            if len(list(tax_id_dict.keys())) == 0:
                logger.warning(f"Could not retrieve tax data from NCBI for {genus}")
                failed_organisms.add(genus)
                continue

            query_key, web_env = post_to_entrez([genus], args)

            if query_key is None:
                failed_organisms.add(genus)
                continue

            lineage_dict, success = get_lineage([genus], lineage_dict, query_key, web_env, args)

            if success is False:
                failed_organisms.add(genus)
                continue

    if len(failed_organisms) != 0:
        with open(failed_log, "w") as fh:
            for genus in failed_organisms:
                fh.write(f"{genus}\n")
    
    return lineage_dict


if __name__ == "__main__":
    main()
