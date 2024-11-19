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
"""Get lineage database from NCBI"""


import json
import logging
import sys

import pandas as pd

from datetime import datetime
from http.client import IncompleteRead
from typing import List, Optional

from Bio import Entrez
from Bio.Entrez.Parser import NotXMLError
from saintBioutils.utilities.logger import config_logger
from saintBioutils.utilities import file_io
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.genbank import entrez_retry
from tqdm import tqdm

from cazy_webscraper import closing_message
from cazy_webscraper.cazy_scraper import connect_existing_db
from cazy_webscraper.expand import get_chunks_list
from cazy_webscraper.ncbi import post_ids
from cazy_webscraper.ncbi.taxonomy.lineage import (
    retry_tax_retrieval,
    get_lineage_data,
)
from cazy_webscraper.sql import sql_orm, sql_interface
from cazy_webscraper.sql.sql_interface.add_data.add_ncbi_tax_data import (
    add_ncbi_taxonomies,
    update_genbank_ncbi_tax,
)
from cazy_webscraper.sql.sql_interface.get_data.get_selected_gbks import get_ids, get_genbank_accessions
from cazy_webscraper.sql.sql_interface.get_data.get_records import get_user_uniprot_sequences
from cazy_webscraper.sql.sql_interface.get_data.get_table_dicts import (
    get_no_tax_gbk_table_dict,
    get_gbk_table_dict,
    get_uniprot_table_dict,
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

    Entrez.email = args.email

    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)
    logger.info(f"Open connection to local cazyme database: {str(args.database)}")

    logger.info("Adding log of scrape to the local CAZyme database")
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
            ec_filters,  # ec_filters not applied when scraping CAZy
            'NCBI Taxonomy',
            'NCBI taxonomc lineages',
            session,
            args,
        )

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        logger.info("Building cache dir")
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        logger.info("Building cache dir")
        cache_dir = cache_dir / "ncbi_tax_retrieval"
        make_output_directory(cache_dir, args.force, args.nodelete_cache)

    # Use Cache
    if args.use_lineage_cache is not None:
        gbk_dict, tax_ids, prot_id_dict = None, None, None
        logger.info("Adding cached lineages to local CAZyme db")
        try:
            with open(args.use_lineage_cache, "r") as fh:
                tax_prot_dict = json.load(fh)
        except FileNotFoundError:
            logger.error(
                "Could not find lineage cache at \n"
                f"{str(args.use_lineage_cache)}\n"
                "Check path is correct\n"
                "Terminating program"
            )
            sys.exit(1)

    # Don't use cache
    else:
        gbk_dict = get_db_proteins(
            class_filters,
            family_filters,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            connection,
            args,
        )

        tax_ids, prot_id_dict = get_ncbi_ids(gbk_dict, cache_dir, args)

        tax_prot_dict = get_lineage_protein_data(tax_ids, prot_id_dict, gbk_dict, cache_dir, args)
        # {tax_id: {linaege info, 'proteins' {local db protein ids}}

        cache_taxonomy(tax_prot_dict, cache_dir)

    add_ncbi_taxonomies(tax_prot_dict, connection, args)
    logger.info("Added lineage data to db")

    update_genbank_ncbi_tax(tax_prot_dict, connection, args)
    logger.info("Added lineage data to protein records in db")

    closing_message("Get NCBI Taxonomy data", start_time, args)


def cache_taxonomy(tax_prot_dict, cache_dir):
    """Cache retrieved taxonomy data

    :param tax_prot_dict: dict, {tax_id: {linaege info, 'proteins' {local db protein ids}}}
    :param cache_dir: Path, path to cache dir

    Return nothing
    """
    logger = logging.getLogger(__name__)

    cache_tax_dict = {}

    for tax_id in tax_prot_dict:
        cache_tax_dict[tax_id] = {}
        for key in tax_prot_dict[tax_id]:
            if key == 'proteins':
                cache_tax_dict[tax_id]['proteins'] = list(tax_prot_dict[tax_id]['proteins'])
            else:
                cache_tax_dict[tax_id][key] = tax_prot_dict[tax_id][key]

    with open((cache_dir/"lineage_data.json"), "w") as fh:
        json.dump(cache_tax_dict, fh)
    logger.info(f"Cached lineage data")


def get_db_proteins(
    class_filters,
    family_filters,
    kingdom_filters,
    taxonomy_filter_dict,
    ec_filters,
    connection,
    args,
):
    """Retrieve proteins from local CAZyme database matching user criteria.

    :param connection: open connection to SQLite db engine
    :param args: cmd-line args parser

    Return dict {}
    """
    logger = logging.getLogger(__name__)

    gbk_dict = {}

    if args.genbank_accessions is not None:
        logger.warning(f"Getting GenBank accessions from file: {args.genbank_accessions}")
        with open(args.genbank_accessions, "r") as fh:
            lines = fh.read().splitlines()

        accessions = [line.strip() for line in lines]
        accessions = set(accessions)

        gbk_dict = get_ids(accessions, connection)

    if args.uniprot_accessions is not None:
        logger.warning(
            "Extracting protein sequences for UniProt "
            f"accessions listed in {args.uniprot_accessions}"
        )
        gbk_table_dict = get_gbk_table_dict(connection)
        uniprot_table_dict = get_uniprot_table_dict(connection)

        gbk_dict.update(get_user_uniprot_sequences(gbk_table_dict, uniprot_table_dict, args))

    if args.genbank_accessions is None and args.uniprot_accessions is None:
        # get user config data
        gbk_dict.update(get_genbank_accessions(
            class_filters,
            family_filters,
            taxonomy_filter_dict,
            kingdom_filters,
            ec_filters,
            connection,
        ))

    if args.update_gbk is False:
        # filter gbk accessions to those proteins with no ncbi tax data
        gbk_db_ids = get_no_tax_gbk_table_dict(connection)

        filtered_gbk_dict = {}

        for gbk_acc in gbk_dict:
            if gbk_dict[gbk_acc] in gbk_db_ids:
                filtered_gbk_dict[gbk_acc] = gbk_dict[gbk_acc]

        gbk_dict = filtered_gbk_dict

    return gbk_dict


def get_ncbi_ids(gbk_dict, cache_dir, args):
    """Get NCBI Protein and Tax IDs from cache and/or NCBI

    :param gbk_dict: dict of gbk accession and local db IDs
    :param cahce_dir: Path to cache dir
    :param args: cmd-line args parser

    Return set of NCBI Tax IDs, and dict of NCBI Prot ID valued by Gbk sequence accession
    """
    tax_ids, prot_id_dict = set(), {}
    logger = logging.getLogger(__name__)

    if (args.use_tax_ids is not None) and (args.use_protein_ids is not None):
        # only use ids in files
        tax_ids = load_tax_ids(args.use_tax_ids)
        prot_id_dict = load_protein_ids(args.use_protein_ids)

        return tax_ids, prot_id_dict

    tax_ids, prot_id_dict = get_ncbi_tax_prot_ids(list(gbk_dict.keys()), cache_dir, args)

    if args.use_tax_ids is not None:
        tax_ids = load_tax_ids(args.use_tax_ids)

    if args.use_protein_ids is not None:
        prot_id_dict = load_protein_ids(args.use_protein_ids)

    logger.info("caching retrieved NCBI Taxonomy and Protein IDs")
    with open((cache_dir/"tax_ids.out"), "a") as fh:
        for tax_id in tax_ids:
            fh.write(f"{tax_id}\n")

    if prot_id_dict is not None:
        with open((cache_dir/"protein_ncbi_ids.out"), "a") as fh:
            for ncbi_prot_id in prot_id_dict:
                fh.write(f"{ncbi_prot_id}\t{prot_id_dict[ncbi_prot_id]}\n")

    return tax_ids, prot_id_dict


def load_tax_ids(tax_id_path):
    """Load NCBI taxonomy ids from plain text file

    :param tax_id_path: Path, path to file of tax ids

    Return list of tax ids
    """
    logger = logging.getLogger(__name__)
    logger.warning(f"Use Tax IDs in {tax_id_path}")
    try:
        with open(tax_id_path, "r") as fh:
            tax_ids = [_.replace("\n", "") for _ in fh.read().splitlines() if _ != "\n"]
    except FileNotFoundError:
        logger.error(
            f"Could not find Tax ID cache at:\n{str(tax_id_path)}\n"
            "Check path is correct\n"
            "Terminating program"
        )
        sys.exit(1)

    return tax_ids


def load_protein_ids(prot_id_path):
    """Load protein IDs and accessions from tab delimited file

    :param prot_id_path: Path, path to tab delimited file

    Return dict {ID: accession}
    """
    logger = logging.getLogger(__name__)
    logger.warning(f"Use Protein IDs and accessions in {prot_id_path}")
    try:
        with open(prot_id_path, "r") as fh:
            cached_prot_data = fh.read().splitlines()
    except FileNotFoundError:
        logger.error(
            f"Could not find NCBI Protein ID cache at:\n{str(prot_id_path)}\n"
            "Check path is correct\n"
            "Terminating program"
        )
        sys.exit(1)

    prot_id_dict = {}
    for line in cached_prot_data:
        if line == "\n":
            continue
        prot_id_dict[line.split("\t")[0].strip()] = line.split("\t")[1].strip().replace("\n", "")

    prot_id_dict


def get_ncbi_tax_prot_ids(protein_accessions, cache_dir, args):
    """Retrieve NCBI Taxonomy and Protein database IDs

    Post the protein accession.
    Then entrez.efetch the protein records to associate each accession with the correct id
    Then entrez.elink to retrieve all linked taxonomy records - do no associated with the correct 
    protein accesion at this point.

    :param protein_accessions: list of NCBI protein version accessions
    :param cache_dir: Path, path to cache directory
    :param args: cmd-line args parser

    Return set of NCBI Tax ids and dict {ncbi prot id: prot acc}
    """
    logger = logging.getLogger(__name__)

    # break up list into batches
    batches = get_chunks_list(protein_accessions, args.batch_size)

    failed_batches = {}  # {protein: tries}

    tax_ids = set()  # set of all retrieve tax ids

    protein_ncbi_ids = {}  # {ncdi prot id : ncbi prot acc}

    for batch in tqdm(batches, desc="Getting NCBI Protein and Taxonomy IDs"):
        # post the IDs
        try:
            query_key, web_env = post_ids(batch, 'Protein', args)
        except RuntimeError as err:
            logger.warning(
                "The batch of protein accessions contains an accession not in NCBI\n"
                "Will parse the accessions individually later to identify the offending accession\n"
                f"Batch:\n{batch}"
            )
            for protein in batch:
                failed_batches[protein] = 0
            continue
        
        if query_key is None:
            for protein in batch:
                failed_batches[protein] = 0
            continue
        
        # fetch protein ids - associate prot id with prot acc
        new_protein_ids = get_prot_ids(query_key, web_env, args)

        if new_protein_ids is None:
            for protein in batch:
                failed_batches[protein] = 0
            continue

        protein_ncbi_ids.update(new_protein_ids)

        new_tax_ids = link_prot_taxs(query_key, web_env, args, batch)

        if new_tax_ids is None:
            for protein in batch:
                failed_batches[protein] = 0
            continue
        
        tax_ids = tax_ids.union(new_tax_ids)

    # retry failed batches
    if len(list(failed_batches.keys())) != 0:
        logger.warning(
            f"Retrying to retrieve Protein IDs and Taxonomy IDs for {len(list(failed_batches.keys()))} proteins"
        )
        while len(list(failed_batches.keys())):
            for protein in failed_batches:
                # post ids
                try:
                    query_key, web_env = post_ids([protein], 'Protein', args)
                except RuntimeError as err:
                    logger.warning(
                        f"The protein accessions {protein} is not listed in NCBI\n"
                        "Therefore not retrieve taxonomy data for this protein.\n"
                        f"Error:\n{err}"
                    )
                    del failed_batches[protein]
                    continue
                
                if query_key is None:
                    failed_batches[protein] += 1
                    if failed_batches[protein] >= args.retries:
                        logger.warning(
                            f"Ran out of connection reattempt for protein accession {protein} while posting IDs\n"
                            "Therefore not retrieve taxonomy data for this protein."
                        )
                        del failed_batches[protein]
                    continue

                # fetch prot ids - linked to acc
                new_protein_ids = get_prot_ids(query_key, web_env, args)

                if new_protein_ids is None:
                    failed_batches[protein] += 1
                    if failed_batches[protein] >= args.retries:
                        logger.warning(
                            f"Ran out of connection reattempt for protein accession {protein} while fetching Protein IDs\n"
                            "Therefore not retrieve taxonomy data for this protein."
                        )
                        del failed_batches[protein]
                    continue

                protein_ncbi_ids.update(new_protein_ids)

                # fetch linked tax ids - not linked to acc
                new_tax_ids = link_prot_taxs(query_key, web_env, args, batch)

                if new_tax_ids is None:
                    failed_batches[protein] += 1
                    if failed_batches[protein] >= args.retries:
                        logger.warning(
                            f"Ran out of connection reattempt for protein accession {protein} while fetching Tax IDs\n"
                            "Therefore not retrieve taxonomy data for this protein."
                        )
                        del failed_batches[protein]
                    continue
                
                tax_ids = tax_ids.union(new_tax_ids)

    failed_prot_accs = [_ for _ in protein_accessions if _ not in set(protein_ncbi_ids.values())]

    if len(failed_prot_accs) != 0:
        cache_file = cache_dir / "failed_to_retrieve_ids"
        logger.warning(
            f"Attempted to retrieve Protein IDs for {len(protein_accessions)} proteins accessions\n"
            f"Retrieved {len(list(protein_ncbi_ids.keys()))} protein IDs\n"
            f"Therefore, protein IDs were not retrieved for {len(failed_prot_accs)} protein accessions\n"
            f"Writing accessions to cache file:{cache_file}"
        )

        with open(cache_file,'w') as fh:
            for acc in failed_prot_accs:
                fh.write(f"{acc}\n")
    
    else:
        logger.warning(
            f"Success: Attempted to retrieve Protein IDs for {len(protein_accessions)} proteins accessions\n"
            f"and retrieved {len(list(protein_ncbi_ids.keys()))} protein IDs\n"
            f"Therefore, retrieved a protein ID for all protein accessions"
        )

    return tax_ids, protein_ncbi_ids


def get_prot_ids(query_key, web_env, args):
    """Retrieve protein records ids from the NCBI Protein database

    :param query_key: str, query key from entrez.epost of protein accessions
    :param web_env: str, web environment from entrez.epost of protein accessions
    :param args: CLI args parser

    Return dict {protein ID: protein accession} or None if fails
    """
    logger = logging.getLogger(__name__)
    new_prot_ids = {}  # {id: acc}
    try:
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            query_key=query_key,
            WebEnv=web_env,
            db="Protein",
            rettype='docsum',
        ) as handle:
            prot_records = Entrez.read(handle, validate=False)

    except (TypeError, AttributeError, RuntimeError,  NotXMLError, IncompleteRead) as err:
        logger.warning(
            "Could not rerieve records from the NCBI Protein database.\n"
            "Will retry individual records later\n"
            f"Error message retrieved:\n{err}"
        )
        return None

    for record in prot_records:
        prot_id = record['Id']
        prot_acc = record['AccessionVersion']
        new_prot_ids[prot_id] = prot_acc

    return new_prot_ids


def link_prot_taxs(query_key, web_env, args, batch):
    """Get NCBI tax and protein IDs from elinking protein accs to tax records

    :param query_key: str, from Entrez epost
    :param web_env: str, from Entrez epost
    :param args: cmd-line args parser
    :param batch: list of ncbi protein accs - used in error messages only

    Return set of NCBI tax ids or None if fails
    """
    logger = logging.getLogger(__name__)

    tax_ids = set()

    try:
        with entrez_retry(
            args.retries,
            Entrez.elink,
            query_key=query_key,
            WebEnv=web_env,
            dbfrom="Protein",
            db="Taxonomy",
            linkname="protein_taxonomy",
        ) as handle:
            tax_links = Entrez.read(handle, validate=False)
    except RuntimeError as err:
        logger.warning(
            "When retrieving taxonomy IDs, the batch included invalid NCBI protein ids\n"
            "Will try batch components individually later, to identify the offending id\n"
            f"Batch:\n{batch}\n"
            f"Error:\n{err}"
        )
        return None
    except (TypeError, AttributeError, NotXMLError, IncompleteRead) as err:
        logger.warning(
            "Failed to connect to NCBI.Entrez and download the Entrez.elink result\n"
            "Will retry later"
            f"Batch:\n{batch}\n"
            f"Error:\n{err}"
        )
        return None


    for result in tax_links:
        for item in result['LinkSetDb']:
            links = item['Link']
            for link in links:
                tax_ids.add(link['Id'])

    return tax_ids


def get_lineage_protein_data(tax_ids, prot_id_dict, gbk_dict, cache_dir, args):
    """Retrieve lineage data from NCBI Taxonomy and associated tax ids with protein ids

    :param tax_ids: set of ncbi tax record ids
    :param prot_id_dict: dict {prot_id: prot_acc}
    :param gbk_dict: dict {protein acc: local db id}
    :param cache_dir: path to cache directory
    :param args: cmd-line args parser

    Return dict {tax_id: {linaege info, 'proteins' {local db protein ids}}
    """
    logger = logging.getLogger(__name__)
    tax_prot_dict = {}  # {ncbi tax id: {phylo_rank: str, proteins: [NCBI protein IDs]}}
    all_failed_ids = set()  # tax ids for whom data could not be retrieved

    # retrieve lineage data from the NCBI Taxonomy database
    lineage_dict, failed_ids = get_lineage_data(tax_ids, args)
    logger.warning(
        f"Queried NCBI with {len(tax_ids)} tax ids\n"
        f"Retrieved lineages for {len(list(lineage_dict.keys()))} tax ids"
    )

    if len(list(lineage_dict.keys())) == 0:
        return lineage_dict

    all_failed_ids = all_failed_ids.union(failed_ids)

    # cache the retrieved lineages
    with open((cache_dir/"ncbi_lineages.json"), "w") as fh:
        json.dump(lineage_dict, fh)

    # retrieve proteins linked to taxon record in NCBI - for only tax ids where lineage data was retrieved
    tax_prot_dict = {}  # {tax_id: {local db protein ids}}
    failed_linked_ids = {}

    # must query each tax id individually
    # Batch querying Entrez.elink combines all the protein ids together for all tax ids
    # therefore, could not separate out the protein ids for each tax id
    # batch query returns ... DbFrom': 'taxonomy', 'IdList': ['55207', '204038', '1421545']}]
    for tax_id in tqdm(lineage_dict, desc="Link proteins to NCBI Tax record"):    
        new_tax_prot_dict, success = get_tax_proteins(
            tax_id,
            tax_prot_dict,
            prot_id_dict,
            gbk_dict,
            cache_dir,
            args,
        )
        # {tax_id: {local db protein ids}}
        tax_prot_dict.update(new_tax_prot_dict)        

        if success is False:
            failed_linked_ids[tax_id] = 1  # first attempt to conenct failed

    if len(failed_linked_ids) != 0:
        logger.info(f"Retrying retrieving linked proteins for failed NCBI Tax IDs")
        while len(list(failed_linked_ids.keys())) > 0:
            for tax_id in tqdm(failed_linked_ids, desc="Retrying retrieving linked proteins"):
                lineage_dict, success = get_tax_proteins(
                    tax_id,
                    tax_prot_dict,
                    prot_id_dict,
                    gbk_dict,
                    cache_dir,
                    args,
                )

                if success is False:
                    failed_linked_ids[tax_id] += 1  # number of attempts that have failed

                if failed_linked_ids[tax_id] > args.retries:
                    logger.warning(
                        f"Ran out of reattempts to get linked proteins data for ncbi tax {tax_id}"
                    )
                    del failed_linked_ids[tax_id]
                    all_failed_ids.add(f"{tax_id}\tCould not retrieve linked proteins")

    # combine lineage data and proteins into a single dict
    for tax_id in tax_ids:
        try:
            lineage_dict[tax_id]['proteins'] = tax_prot_dict[tax_id]
        except KeyError:
            logger.error(
                f"Did not retrieve lineage and/or linked proteins for ncbi tax id {tax_id}"
            )
            if tax_id not in all_failed_ids:
                all_failed_ids.add(
                    f"{tax_id}\tCould not lineage data and/or retrieve linked proteins"
                )

    if len(all_failed_ids) != 0:
        with open((cache_dir / "failed_tax_ids.txt"), "a") as fh:
            for tax_id in all_failed_ids:
                fh.write(f"{tax_id}\n")

    return lineage_dict


def get_tax_proteins(tax_id, tax_prot_dict, prot_id_dict, gbk_dict, cache_dir, args):
    """Get the proteins linked to a tax id in NCBI, and link the tax id with the local db protein ids

    :param tax_id: str, NCBI tax db id
    :param tax_prot_dict: {ncbi tax id: {local db protein ids}}
    :param prot_id_dict: dict {protein ncbi id: prot acc}
    :param gbk_dict: dict, {prot acc: local db id}
    :param cache_dir: Path, path to cache dir
    :param args: cmd-line args parser

    Return dict {tax_id: {local db protein ids}} and bool (True=success, False=failed)
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
        return tax_prot_dict, False

    try:
        tax_prot_dict[tax_id]
    except KeyError:
        tax_prot_dict[tax_id] = set()

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

                try:
                    prot_local_db_id = gbk_dict[prot_ver_acc]
                except KeyError:
                    logger.error(
                        "Did not previously retrieved data from the local "
                        f"db for {prot_local_db_id}\n"
                        "Caching and skipping protein"
                    )
                    with open((cache_dir/"failed_local_db_retrieval.out"), "a") as fh:
                        fh.write(f"{prot_local_db_id}\n")
                    continue

                tax_prot_dict[tax_id].add(prot_local_db_id)

    return tax_prot_dict, True


def check_for_failed_accessions(tax_prot_dict, gbk_dict, cache_dir):
    """Check to see for how many proteins, taxonomy inforamtion could not be retrieved from NCBI
    
    :param tax_prot_dict: {tax_id: {linaege info, 'proteins' {local db protein ids}}
    :param gbk_dict: dict of ncbi accessions paired with local db IDs
    :param cache_dir: Path to cache directory
    
    Return list of accessions for whom tax data was not retrieved from ncbi
    """
    logger = logging.getLogger(__name__)
    protein_accessions = list(gbk_dict.keys())

    # tax_prot_dict = {tax_id: {linaege info, 'proteins' {local db protein ids}}
    retrieved_protein_acc = []
    for tax_id in tax_prot_dict:
        for protein in tax_prot_dict[tax_id]['proteins']:
            db_id = protein
            for acc in gbk_dict:
                if gbk_dict[acc] == db_id:
                    protein_acc = acc
                    break
            
            retrieved_protein_acc.append(protein_acc)

    failed_accessions = [_ for _ in protein_accessions if _ not in retrieved_protein_acc]
    if len(failed_accessions) != 0:
        cache_file = cache_dir / "failed_to_retrieve_tax_data"
        logger.warning(
            f"Attempted to retrieve tax data for {len(protein_accessions)} proteins accessions\n"
            f"Retrieved tax data for {len(retrieved_protein_acc)} proteins\n"
            f"Therefore, tax data was not retrieved for {len(failed_accessions)} protein accessions\n"
            f"Writing accessions to cache file:{cache_file}"
        )

        with open(cache_file,'w') as fh:
            for acc in failed_accessions:
                fh.write(f"{acc}\n")
    
    else:
        logger.warning(
            f"Success: Attempted to retrieve tax data for {len(protein_accessions)} proteins accessions\n"
            f"and retrieved tax data for {len(retrieved_protein_acc)} proteins"
        )

    return failed_accessions


def retry_failed_linkages(
    class_filters,
    family_filters,
    kingdom_filters,
    taxonomy_filter_dict,
    ec_filters,
    connection,
    prot_id_dict,
    cache_dir,
    tax_prot_dict,
    args,
):
    """Retry connecting proteins to taxonomy as single linkages
    Used mostly for debugging retrieval of taxonomy info for all protein accessions
    """
    logger = logging.getLogger(__name__)
    if gbk_dict is None:
        logger.warning("Loading Genbanks table into dict to check for failed retrievals of tax data")
        gbk_dict = get_db_proteins(
            class_filters,
            family_filters,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            connection,
            args,
        )
        # gbk_dict = {gbk acc: db id}

    # check to see for how many proteins a ncbi taxonomy record could not be retrieved
    failed_accessions = check_for_failed_accessions(tax_prot_dict, gbk_dict, cache_dir)

    if len(failed_accessions) != 0:
        # check if failure was due to not getting a protein id from ncbi
        
        if prot_id_dict is not None:
            failed_ids = [_ for _ in failed_accessions if _ not in list(prot_id_dict.values())]
        else:
            failed_ids = []
        
        if len(failed_ids) != 0:
            logger.warning(f"Did not retrieved protein IDs for {len(failed_ids)} proteins:\n{failed_ids}")

        # tax_prot_dict = {tax_id: {linaege info, 'proteins' {local db protein ids}}
        # prot_id_dict = {ncbi prot id: prot acc}
        tax_prot_dict = retry_tax_retrieval(failed_accessions, failed_ids, prot_id_dict, tax_prot_dict, args, cache_dir)

    return tax_prot_dict


if __name__ == "__main__":
    main()
