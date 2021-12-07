#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

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
"""Explore the number of GenBank genomes annotated by CAZy."""


import logging

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio import Entrez
from saintBioutils.genbank import entrez_retry
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from cazy_webscraper import cazy_scraper


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # Program preparation
    if argv is None:
        parser = genbank_parser.build_parser()
        args = parser.parse_args()
    else:
        parser = genbank_parser.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)

    Entrez.email = args.email


    # check if need to build output dir
    ???

    # connect to the local CAZyme database
    connection, logger_name, cache_dir = cazy_scraper.connect_existing_db(args, time_stamp)

    # load Genbank and Kingdom records from the db
    genbank_kingdom_dict = get_gbk_kingdom_dict(connection)

    genomic_assembly_names = get_assebmly_names(genbank_kingdom_dict, args)

    genomic_accession_dict = get_genomic_accessions(genomic_assembly_names, args)

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    if args.verbose:
        logger.info(
            "Finished getting data from UniProt\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_scraper.VERSION_INFO}\n"
            f"Citation: {cazy_scraper.CITATION_INFO}"
        )
    else:
        print(
            "=====================cazy_webscraper=====================\n"
            "Finished getting data from UniProt\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_scraper.VERSION_INFO}\n"
            f"Citation: {cazy_scraper.CITATION_INFO}"
        )


def get_gbk_kingdom_dict(connection):
    """Compile dict of Genbank and Kingdom records
    
    :param connection: open sqlalchemy db connection
    
    Return dict {kingdom: {genbank_accessions}}
    """
    # load the Genbank and Kingdom records from the local CAZyme db
    genbank_kingdom_records = []

    # extract the GenBank accessions and Kingdoms from the records
    genbank_kingdom_dict = {}  # kingdom: {genbank_accessions}
    for record in genbank_kingdom_records:
        kingdom = record[1].kingdom
        gbk_acc = record[0].genbank_accession

        try:
            genbank_kingdom_dict[kingdom].add(gbk_acc)
        except KeyError:
            genbank_kingdom_dict[kingdom] = gbk_acc

    return genbank_kingdom_dict


def get_assebmly_names(genbank_kingdom_dict, args):
    """Retrieve assembly names of the source genomic accessions for the protein accessions in the local db.
    
    :param genbank_kingdom_dict: dict of Genbank and Kingdom records from db
        {kingdom: {genomic_accessions}}
    :param args: cmd-line args parser
    
    Return dict {kingdom: {genomic_acc: {'proteins': [p_acc], 'count':#ofProteins}}}
    """
    logger = logging.getLogger(__name__)

    genomic_assembly_names = {}  # {kingdom: {assembly_name: {protein_accessions}}}
    
    for kingdom in tqdm(genbank_kingdom_dict, desc="Retrieving genomic assembly names per kingdom"):
        gbk_accessions = genbank_kingdom_dict[kingdom]

        # break up the list into a series of smaller lists that can be batched querried
        batch_queries = []

        for batch_query in batch_queries:
            batch_query_ids = ",".join(batch_query)
            with entrez_retry(
                args.retries, Entrez.epost, "Protein", id=batch_query_ids,
            ) as handle:
                batch_post = Entrez.read(handle)


            # eFetch against the PRotein database, retrieve in xml retmode
            with entrez_retry(
                args.retries,
                Entrez.efetch,
                db="Protein",
                query_key=batch_post['QueryKey'],
                WebEnv=batch_post['WebEnv'],
                retmode="xml",
            ) as handle:
                batch_fetch = Entrez.read(handle)
            
            for record in tqdm(batch_fetch, desc="Retrieving assembly name from protein records"):
                assembly_name = None
                protein_accession = record['GBSeq_accession-version']

                for item in record['GBSeq_comment'].split(";"):
                    if item.strip().startswith("Assembly Name ::"):
                        assembly_name = item.strip()[len("Assembly Name :: "):]

                if assembly_name is None:
                    logger.warning(
                        f"Could not retrieve genome assembly name for {protein_accession}"
                    )
                    continue
                
                try:
                    genomic_assembly_names[kingdom]

                    try:
                        genomic_assembly_names[kingdom][assembly_name].add(protein_accession)
                    except KeyError:
                        genomic_assembly_names[kingdom][assembly_name] = {protein_accession}

                except KeyError:
                    genomic_assembly_names[kingdom] = {assembly_name: {protein_accession}}

    return genomic_assembly_names


def get_genomic_accessions(genomic_assembly_names, args):
    """Retrieve genomic accessions for the genomic assemblies

    :param genomic_assembly_names: dict
        {kingdom: {genomic_acc: {'proteins': [p_acc], 'count':#ofProteins}}}
    :param args: cmd-line args parser
    
    Return dict, {kingdom: {genomic_acc: {'proteins': [p_acc], 'count':#ofProteins}}}"""

    genomic_accession_dict = {}  #  {kingdom: {genomic_acc: {'proteins': [p_acc], 'count':#ofProteins}}}

    for kingdom in tqdm(genomic_assembly_names, desc='Retrieving genomic accessions per kingdom'):
        assembly_names = list(genomic_assembly_names[kingdom].keys())

        # break up the list into a series of smaller lists that can be batched querried
        batch_queries = []

        for batch_query in batch_queries:
            batch_query_ids = ",".join(batch_query)

        # retrieve the records IDs for the assembly names
        with entrez_retry(
            args.retries,
            Entrez.esearch,
            "Assembly",
            id=batch_query_ids,
        ) as handle:
            batch_post = Entrez.read(handle)
        
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            db="Assembly",
            query_key=batch_post['QueryKey'],
            WebEnv=batch_post['WebEnv'],
            retmode="xml",
        ) as handle:
            batch_fetch = Entrez.read(handle)

        genomic_accessions = {}

        for genome_record in tqdm(batch_fetch, desc="Retrieving assembly IDs"):
            index = 0
            accessions = set()

            for index in range(len(genome_record['IdList'])):
                with entrez_retry(
                    10, Entrez.efetch, db="Assembly", id=genome_record['IdList'][index], retmode="xml", rettype="docsum",
                ) as handle:
                    result = Entrez.read(handle)
                genomic_accession = result['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
                assembly_name = result['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']

                genomic_accessions[genomic_accession] = assembly_name

            accessions = list(genomic_accessions.keys)
            accessions.sort(reverse=True)
            latest_accession = accessions[0]
            latest_assembly_name = genomic_accessions[latest_accession]

            # replace assemlby name for genomic accession
            try:
                protein_accessions = genomic_assembly_names[kingdom][latest_assembly_name]

                try:
                    genomic_accession_dict[kingdom]

                except KeyError:
                    genomic_accession_dict[kingdom] = {genomic_accession = {
                        "proteins": protein_accessions,
                        "count": 
                    }}

            except KeyError:
                logger.warning(f"Retrieved assembly name {latest_assembly_name}, but not retrieved previously")

    return genomic_accession_dict


if __name__ == "__main__":
    main()
