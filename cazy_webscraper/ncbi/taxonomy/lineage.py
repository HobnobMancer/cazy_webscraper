#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Parse taxonomy data from NCBI and CAZy to replace when multiple taxa are retrieved from CAZy"""


import logging

from copy import copy

from saintBioutils.genbank import entrez_retry
from tqdm import tqdm

from Bio import Entrez


def fetch_lineages(tax_ids, query_key, web_env, args):
    """Fetch lineage data from NCBI Taxonomy database.

    :param tax_ids: list of NCBI tax ids
    :param query_key: str, query key from ncbi.entrez.epost
    :param web_env: str, web environment from ncbi.entrez.epost
    :param args: CLI args parser

    Return dict of lineage data {ncbi_tax_id: {rank: str/lineage}}
        or None if connection fails
    """
    lineage_dict = {}  # {ncbi_tax_id: {rank: str/lineage}}

    try:
        with entrez_retry(
            10,
            Entrez.efetch,
            db="Taxonomy",
            query_key=query_key,
            WebEnv=we,
        ) as handle:
            batched_tax_records = Entrez.read(handle, validate=False)

    except (TypeError, AttributeError) as err:
        logger.warning(f"Failed to fetch tax record from NCBI tax for id '{tax_id}'':\n{err}")
        return

    for record in batched_tax_records:
        record_id = record['TaxId']
        if record_id not in tax_ids:
            continue

        # set all lineage data to None
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

        if genus is not None:
            # drop genus from species name
            if species is not None:
                species = species.replace(genus, "").strip()

            # extract strain from scientific name if not retrieved as rank
            if species is not None and strain is None:
                strain = scientific_name.replace(f"{genus} {species}", "").strip()

            # extract species from the scientific name if not retrieved as rank
            elif species is None:
                species = scientific_name.replace(genus, "").strip()

        lineage_dict[record_id] = {
            'kingdom': kingdom,
            'phylum': phylum,
            'class': tax_class,
            'order': order,
            'family': family,
            'genus': genus,
            'species': species,
            'strain': strain,
        }

    return lineage_dict


def get_taxid_lineages(batches, tax_ids, args):
    """Batch query NCBI to retrieve lineage data for each tax_id

    :param batches: list of nested lists
    :param tax_ids: list of str (ncbi tax id)
    :param args: CLI args parser

    return dict {ncbi tax id: {rank: str}}
    """
    failed_batches = {}  # {str(batch): {'batch': list(batch), 'tries': int(num of failed attempts)}}
    unlisted_id_batches = []

    tax_lineage_dict = {}  # {ncbi_tax_id: {rank: str/lineage}}

    for batch in tqdm(batches, desc="Batch retrieving lineage data from NCBI"):
        # post the ids
        try:
            query_key, web_env = post_ids(batch, 'Taxonomy', args)
        except (RuntimeError) as err:
            logger.warning(
                "RunTime error occurred when querying NCBI.\n"
                "Likely caused by batch containing invalid ID.\n"
                "Will try individual IDs later\n"
                f"Error message:\n{err}"
            )
            unlisted_id_batches += batch
            continue
        
        if query_key is None:
            failed_batches[str(batch)] = {'batch': batch, 'tries': 1}
            continue
        
        # fetch the lineage data
        try:
            lineage_dict = fetch_lineages(tax_ids, query_key, web_env, args)  # {ncbi_tax_id: {rank: str/lineage}}
        except RuntimeError as errr:
            logger.warning(
                "RunTime error occurred when querying NCBI.\n"
                "Likely caused by batch containing invalid ID.\n"
                "Will try individual IDs later\n"
                f"Error message:\n{err}"
            )
            unlisted_id_batches += batch
            continue
            
        if lineage_dict is None:  # failed connection
            failed_batches[str(batch)] = {'batch': batch, 'tries': 1}
            continue

        tax_lineage_dict.update(lineage_dict)

    return tax_lineage_dict, failed_batches, unlisted_id_batches


def retry_get_tax_lineages(failed_batches, tax_ids, args):
    """Retry batches were the connection files previously

    :param failed_batches: dict, {str(batch): {'batch': list(batch), 'tries': int(num of failed attempts)}}
    :param tax_ids: list of tax ids
    :param args: CLI args parser

    Return 
    * lineage_dict {ncbi tax id: {rank: str}}
    * list of unlisted IDs
    * list of tax ids for who data could not be retrieved from NCBI
    """
    tax_lineage_dict = {}
    batches = [failed_batches[batch_name]['batch'] for batch_name in failed_batches] # [[b1], [b2]]
    unlisted_id_batches = []
    failed_ids = []

    while len(list(failed_batches.keys())) != 0:
        for batch in tqdm(batches, "Retrying failed batches"):
            # check if batch in failed_batches
            try:
                failed_batches[str(batch)]
            except KeyError:
                continue  # batch no longer in failed_batches

            # batch still in failed_batches so need to try to gather data again

            # post the ids
            try:
                query_key, web_env = post_ids(batch, 'Taxonomy', args)
            except (RuntimeError) as err:
                logger.warning(
                    "RunTime error occurred when querying NCBI.\n"
                    "Likely caused by batch containing invalid ID.\n"
                    "Will try individual IDs later\n"
                    f"Error message:\n{err}"
                )
                unlisted_id_batches += batch
                continue
            
            if query_key is None:
                failed_batches[str(batch)]['tries'] += 1

                if failed_batches[str(batch)]['tries'] >= args.retries:
                    logger.warning(
                        f"Reached maximum number of connection reattempts ({args.retries} tries) for batch\n"
                        f"of tax ids to retrieve the lineage data from NCBI:\n{batch}\n"
                        "But could not connect to NCBI.\n"
                        "Therefore, not retrieving lineage data for these taxonomy ids"
                    )
                    del failed_batches[str(batch)]
                    failed_ids += batch
                    continue

            # fetch the lineage data
            try:
                lineage_dict = fetch_lineages(tax_ids, query_key, web_env, args)  # {ncbi_tax_id: {rank: str/lineage}}
            except RuntimeError as errr:
                logger.warning(
                    "RunTime error occurred when querying NCBI.\n"
                    "Likely caused by batch containing invalid ID.\n"
                    "Will try individual IDs later\n"
                    f"Error message:\n{err}"
                )
                unlisted_id_batches += batch
                continue
                
            if lineage_dict is None:  # failed connection
                failed_batches[str(batch)]['tries'] += 1

                if failed_batches[str(batch)]['tries'] >= args.retries:
                    logger.warning(
                        f"Reached maximum number of connection reattempts ({args.retries} tries) for batch\n"
                        f"of tax ids to retrieve the lineage data from NCBI:\n{batch}\n"
                        "But could not connect to NCBI.\n"
                        "Therefore, not retrieving lineage data for these taxonomy ids"
                    )
                    del failed_batches[str(batch)]
                    failed_ids += batch
                    continue

            tax_lineage_dict.update(lineage_dict)

    return tax_lineage_dict, unlisted_id_batches, failed_ids


def parse_unlised_taxid_lineages(batches, tax_ids, args):
    """Batch query NCBI to retrieve lineage data for each tax_id

    :param batches: list of nested lists, one item per nested list
    :param tax_ids: list of str (ncbi tax id)
    :param args: CLI args parser

    return
    * dict {ncbi tax id: {rank: str}}
    * list of unlisted ids
    """
    failed_ids = {}  # {id: int(num of failed attempts)}}
    unlisted_ids = []
    tax_lineage_dict = {}  # {ncbi_tax_id: {rank: str/lineage}}

    all_ids_to_parse = [_[0] for _ in batches]

    for batch in batches:
        failed_ids[batch[0]] = 0  # {id: num of tries}

    while len(all_ids_to_parse) != 0:
        for batch in tqdm(failed_ids, desc="Batch retrieving lineage data from NCBI"):
            # post the ids
            try:
                query_key, web_env = post_ids([batch], 'Taxonomy', args)
            except (RuntimeError) as err:
                logger.warning(
                    "RunTime error occurred when queryig NCBI.\n"
                    f"Likely caused by ID {batch} containing invalid ID.\n"
                    "Will not retrieve linage data for this ID\n"
                    f"Error message:\n{err}"
                )
                unlisted_ids += batch
                all_ids_to_parse.remove(batch)
                del failed_ids[batch]
                continue
            
            if query_key is None:
                failed_ids[batch] += 1

                if failed_ids[batch] >= args.retries:
                    logger.warning(
                        f"Reached maximum number of connection reattempts ({args.retries} tries) when\n"
                        f"retrieving lineage data for tax id {batch}\n"
                        "Could not connect to NCBI. Therefore, not retrieving lineage data for these taxonomy ids"
                    )
                    unlisted_ids += batch
                    all_ids_to_parse.remove(batch)
                    del failed_ids[batch]

                continue
            
            # fetch the lineage data
            try:
                lineage_dict = fetch_lineages(tax_ids, query_key, web_env, args)  # {ncbi_tax_id: {rank: str/lineage}}
            except (RuntimeError) as err:
                    logger.warning(
                        "RunTime error occurred when queryig NCBI.\n"
                        f"Likely caused by ID {batch} containing invalid ID.\n"
                        "Will not retrieve linage data for this ID\n"
                        f"Error message:\n{err}"
                    )
                    unlisted_ids += batch
                    all_ids_to_parse.remove(batch)
                    del failed_ids[batch]
                    continue
                    
            if lineage_dict is None:  # failed connection
                failed_ids[batch] += 1

                if failed_ids[batch] >= args.retries:
                    logger.warning(
                        f"Reached maximum number of connection reattempts ({args.retries} tries) when\n"
                        f"retrieving lineage data for tax id {batch}\n"
                        "Could not connect to NCBI. Therefore, not retrieving lineage data for these taxonomy ids"
                    )
                    unlisted_ids += batch
                    all_ids_to_parse.remove(batch)
                    del failed_ids[batch]

                continue

            tax_lineage_dict.update(lineage_dict)

    return tax_lineage_dict, unlisted_ids
