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

from Bio import Entrez
from saintBioutils.genbank import entrez_retry
from tqdm import tqdm
from http.client import IncompleteRead
from Bio.Entrez.Parser import NotXMLError

from cazy_webscraper.expand import get_chunks_list
from cazy_webscraper.ncbi import post_ids


def link_single_prot_tax(prot_id, args):
    """Get NCBI tax IDs from linking a single NCBI protein ID

    :param prot_id: str, ncbi protein record id
    :param args: cmd-line args parser

    Return set of NCBI tax ids or None if fails
    """
    logger = logging.getLogger(__name__)

    tax_ids = set()

    try:
        with entrez_retry(
            args.retries,
            Entrez.elink,
            id=prot_id,
            dbfrom="Protein",
            db="Taxonomy",
            linkname="protein_taxonomy",
        ) as handle:
            tax_links = Entrez.read(handle, validate=False)
    except RuntimeError as err:
        logger.warning(
            f"When retrieving taxonomy IDs, the NCBI protein id {prot_id} is invalid\n"
            f"Error:\n{err}"
        )
        return "invalid"
    except (TypeError, AttributeError, NotXMLError, IncompleteRead) as err:
        logger.warning(
            "Failed to connect to NCBI.Entrez and download the Entrez.elink result\n"
            "Will retry later"
            f"ID:\n{prot_id}\n"
            f"Error:\n{err}"
        )
        return "connection"


    for result in tax_links:
        for item in result['LinkSetDb']:
            links = item['Link']
            for link in links:
                tax_ids.add(link['Id'])

    return tax_ids
    

def get_lineage_data(tax_ids, args):
    """Coordinate retrieving all lineage data for all tax_ids from NCBI Taxonomy DB

    :param tax_ids: list of tax ids
    :param args: CLI args parser

    Return 
    * {ncbi tax id: {rank: str}}
    * set of tax ids for which data could not be retrieved from NCBI
    """
    logger = logging.getLogger(__name__)
    all_tax_lineage_dict = {}  # {ncbi tax id: {rank: str}}
    og_batches = get_chunks_list(list(tax_ids), args.batch_size)
    all_failed_ids = set()

    tax_lineage_dict, failed_batches, unlisted_id_batches = get_taxid_lineages(
        og_batches,
        tax_ids,
        args,
    )
    all_tax_lineage_dict.update(tax_lineage_dict)

    if len(list(failed_batches.keys())) != 0:
        logger.warning("Retring failed connections")
        tax_lineage_dict, new_unlisted_id_batches, failed_ids = retry_get_tax_lineages(
            failed_batches,
            tax_ids,
            args,
        )
        all_tax_lineage_dict.update(tax_lineage_dict)
        unlisted_id_batches += new_unlisted_id_batches
        for failed_id in failed_ids:
            all_failed_ids.add(failed_id)
    
    if len(unlisted_id_batches) != 0:
        logger.warning("Retrying batches with unlised IDs")
        # makke list of nested lists [[id1], [id2]]
        individual_ids = [[_] for _ in unlisted_id_batches] 

        tax_lineage_dict, failed_ids = parse_unlised_taxid_lineages(
            individual_ids,
            tax_ids,
            args,
        )
        
        all_tax_lineage_dict.update(tax_lineage_dict)
        for tax_id in failed_ids:
            all_failed_ids.add(tax_id)

    return all_tax_lineage_dict, all_failed_ids


def fetch_lineages(tax_ids, query_key, web_env, args):
    """Fetch lineage data from NCBI Taxonomy database.

    :param tax_ids: list of NCBI tax ids
    :param query_key: str, query key from ncbi.entrez.epost
    :param web_env: str, web environment from ncbi.entrez.epost
    :param args: CLI args parser

    Return dict of lineage data {ncbi_tax_id: {rank: str/lineage}}
        or None if connection fails
    """
    logger = logging.getLogger(__name__)
    lineage_dict = {}  # {ncbi_tax_id: {rank: str/lineage}}

    try:
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            db="Taxonomy",
            query_key=query_key,
            WebEnv=web_env,
        ) as handle:
            batched_tax_records = Entrez.read(handle, validate=False)

    except (TypeError, AttributeError) as err:
        logger.warning(f"Failed to fetch tax records from NCBI tax for ids '{tax_ids}'':\n{err}")
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
    logger = logging.getLogger(__name__)
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
    logger = logging.getLogger(__name__)
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
    logger = logging.getLogger(__name__)
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


def retry_tax_retrieval(failed_accessions, failed_ids, prot_id_dict, tax_prot_dict, args, cache_dir):
    """Attempt to link protein IDs to tax records, where information was not retrieved previously
    
    :param failed_accessions: ncbi prot accessions for whom tax data was not retrieved
    :param failed_ids: protein ncbi ids for failed accessions
    :param prot_id_dict: {ncbi prot id: ncbi prot accession}
    :param tax_prot_dict: {tax_id: {linaege info, 'proteins' {local db protein ids}}
    :param args: CLI args parser
    :param cache_dir: path to cache directory
    
    return updated tax_prot_dict
    """
    logger = logging.getLogger(__name__)

    failed_connections = []
    invalid_ids = []
    invalid_tax_ids = []

    for protein_accession in tqdm(failed_accessions, desc="Retrying failed accessions"):
        if protein_accession not in list(prot_id_dict.values()):
            logger.warning(
                f"Protein {protein_accession} not listed in NCBI\n"
                "Not retrieving tax data for this protein"
            )
            continue
    
        prot_id = None
        for temp_id in prot_id_dict:
            if prot_id_dict[temp_id] == protein_accession:
                prot_id = temp_id
                break
        
        if prot_id is None:
            logger.warning(
                f"Could not retrieved local db ID for protein {protein_accession}\n"
                "Not retrieving tax data for this protein"
            )
            continue

        if prot_id in failed_ids:
            logger.warning(
                f"Protein accession {protein_accession} ID {prot_id} listed as invalid\n"
                "Not retrieving tax data for this protein"
            )
            continue
            
        tax_ids = link_single_prot_tax(prot_id, args)
        if tax_ids == "connection":
            failed_connections.append((protein_accession, prot_id))
            logger.warning(
                f"Connection to NCBI for protein accession {protein_accession} failed\n"
                "Not retrieving tax data for this protein"
            )
            continue

        if tax_ids == "invalid":
            invalid_ids.append((protein_accession, prot_id))
            logger.warning(
                f"NCBI for protein accession {protein_accession} (ID:{prot_id}) could not be linked to a tax\n"
                "record in NCBI. Potentially an invalid ID or accession.\n"
                "Not retrieving tax data for this protein"
            )
            continue
        
        taxs_to_retrieve = []
        for tax_id in tax_ids:
            if tax_id not in list(tax_prot_dict.keys()):
                taxs_to_retrieve.append(tax_id)
                continue
            
            tax_prot_dict[tax_id]['proteins'].add(prot_id)
            
        for tax_id in taxs_to_retrieve:
            # retrieve tax data from ncbi
            tax_lineage_dict, failed_id = get_lineage_data([tax_id], args)
            # {ncbi tax id: {rank: str}}

            if len(failed_id) > 0:
                invalid_tax_ids.append((tax_id))
                logger.warning(
                    f"Could not retrieved record from NCBI for NCBI Tax ID {tax_id}\n"
                    f"Not retrieving lineage for linked protein {protein_accession}"
                )
                continue

            # add to tax_prot_dict
            # tax_prot_dict = {tax_id: {linaege info, 'proteins' {local db protein ids}} 
            tax_prot_dict[tax_id] = {'proteins': {prot_id}}

            for rank in tax_lineage_dict[tax_id]:
                tax_prot_dict[tax_id][rank] = tax_lineage_dict[tax_id][rank]

    if len(failed_connections) != 0:
        cache_file = cache_dir / "failed_connections_on_reattempt"
        logger.warning(
            "Writing to cache the protein accessions and ids whose connection to NCBI failed\n"
            f"Cache: {cache_file}"
        )
        with open(cache_file, 'a') as fh:
            for prot_acc_id_pair in failed_connections:
                fh.write(f"{prot_acc_id_pair[0]}\t{prot_acc_id_pair[1]}\n")

    if len(invalid_ids) != 0:
        cache_file = cache_dir / "invalid_prot_ids_on_reattempt"
        logger.warning(
            "Writing to cache the protein ids that could not be linked to a taxonomy record in NCBI\n"
            f"Cache: {cache_file}"
        )
        with open(cache_file, 'a') as fh:
            for prot_id in invalid_ids:
                fh.write(f"{prot_id}\n")

    if len(invalid_tax_ids) != 0:
        cache_file = cache_dir / "invalid_tax_ids_on_reattempt"
        logger.warning(
            "Writing to cache the tax ids whose records could not be retrieved from NCBI\n"
            f"Cache: {cache_file}"
        )
        with open(cache_file, 'a') as fh:
            for tax_id in invalid_tax_ids:
                fh.write(f"{tax_id}\n")

    return tax_prot_dict
