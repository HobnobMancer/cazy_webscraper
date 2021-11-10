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
"""Retrieve proteins sequences from GenBank and populate the local database and write to FASTA"""


import logging
import math
import sys

from tqdm import tqdm

from scraper.expand import get_genbank_sequences
from scraper.expand.get_genbank_sequences.from_sql_db import query_sql_db
from scraper.expand.get_genbank_sequences.ncbi import query_entrez
from scraper.sql.sql_orm import get_db_session
from scraper.utilities import file_io, parse_configuration


def sequences_for_proteins_from_db(date_today, args):
    """Coordinate retrievel of protein sequences for proteins in a SQL database.
    
    :param date_today: str, date script was invoked, used for naming files
    :param args: cmd-line args parser
    
    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # get database session
    try:
        session = get_db_session(args)
    except Exception as err:
        logger.error(
            "Could not connect to local CAZyme database.\nThe following error was raised:\n"
            f"{err}\nTerminating program\n"
        )
        sys.exit(1)

    file_io.make_output_directory(args.fasta, args.force, args.nodelete)

    if args.blastdb is not None:  # build directory to store FASTA file for BLAST db
        file_io.make_output_directory(args.blastdb, args.force, args.nodelete)

    # retrieve configuration data, as a dict of CAZy classes and families to retrieve seqs for
    parse_configuration_path = parse_configuration.__file__
    (
        config_dict, taxonomy_filters, kingdoms, ec_filters,
    ) = parse_configuration.parse_configuration_for_cazy_database(
        parse_configuration_path,
        args,
    )

    genbank_accessions = get_genbank_accessions(
        args,
        session,
        config_dict,
        taxonomy_filters,
        kingdoms,
        ec_filters,
    )

    # break up protein_list into multiple, smaller lists for batch querying Entrez
    # batches of greater than 200 can be rejected by Entrez during busy periods
    # args.epost=size of chunks

    accessions_lists_for_individual_queries = []

    for accession_list in tqdm(
        get_genbank_sequences.get_accession_chunks(genbank_accessions, args.epost),
        desc="Batch retrieving sequences from NCBI",
        total=(math.ceil(len(genbank_accessions) / args.epost)),
    ):
        try:
            accession_list.remove("NA")
        except ValueError:
            pass

        try:
            query_entrez.get_sequences_add_to_db(accession_list, date_today, args)
        except RuntimeError as err:  # typically Some IDs have invalid value and were omitted.
            logger.warning(
                "RuntimeError raised for accession list. Will query accessions individualy after.\n"
                f"The following error was raised:\n{err}"
            )
            with open("legihton_error.txt", "a") as fh:
                fh.write(f"{err}\n{str(accession_list)}")
            accessions_lists_for_individual_queries.append(accession_list)

    if len(accessions_lists_for_individual_queries) != 0:
        for accession_list in tqdm(
            accessions_lists_for_individual_queries,
            desc="Performing individual queries to parse GenBank accessions without records",
        ):
            for accession in tqdm(accession_list, desc="Retrieving individual sequences"):
                try:
                    query_entrez.get_sequences_for_dict([accession], date_today, args)
                except RuntimeError as err:
                    logger.warning(
                        f"Querying NCBI for {accession} raised the following RuntimeError:\n"
                        f"{err}"
                    )
    return

    return


def get_genbank_accessions(args, session, config_dict, taxonomy_filters, kingdoms, ec_filters):
    """Retrieve the GenBank accessions for all proteins for which a sequence will be retrieved.

    :param args: cmd-line argument parser
    :param session: open SQLite db session
    :param config_dict: dict, defines CAZy classes and families to get sequences for
    :param taxonomy_filters: set of genera, species and strains to retrieve sequences for
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param ec_filters: set of EC numbers annotations CAZymes must have at least one to retrieve
        a sequence

    Return a list of GenBank accessions, containing no duplicate GenBank accessions
    """
    logger = logging.getLogger(__name__)

    if config_dict:  # there are specific CAZy classes/families to retrieve sequences for
        genbank_query_class, genbank_query_family = get_cazy_class_fam_genbank_records(
            session,
            config_dict,
        )

        class_genbank_accessions = parse_genbank_query(
            genbank_query_class,
            taxonomy_filters,
            kingdoms,
            ec_filters,
            session,
        )

        family_genbank_accessions = parse_genbank_query(
            genbank_query_family,
            taxonomy_filters,
            kingdoms,
            ec_filters,
            session,
        )

        genbank_accessions = class_genbank_accessions + family_genbank_accessions

    else:
        if args.update:  # retrieve all GenBank accessions

            if args.primary:
                logger.warning(
                    "Retrieving sequences for all PRIMARY GenBank accessions that:\n"
                    "Do not have a sequence in the db OR the sequence has been updated in NCBI"
                )
                genbank_query = query_sql_db.get_all_prim_genbank_acc(session)

            else:
                logger.warning(
                    "Retrieving sequences for all PRIMARY GenBank accessions that\n"
                    "do not have a sequence in the db"
                )
                genbank_query = query_sql_db.get_all_genbank_acc(session)

        else:  # retrieve GenBank accesions of records that don't have a sequence
            if args.primary:
                logger.warning(
                    "Retrieving sequences for all PRIMARY GenBank accessions that\n"
                    "do not have a sequence in the db"
                )
                genbank_query = query_sql_db.get_prim_genbank_accessions_with_no_seq(session)

            else:
                logger.warning(
                    "Retrieving sequences for ALL GenBank accessions that\n"
                    "do not have a sequence in the db"
                )
                genbank_query = query_sql_db.get_genbank_accessions_with_no_seq(session)

        genbank_accessions = parse_genbank_query(
            genbank_query,
            taxonomy_filters,
            kingdoms,
            ec_filters,
            session,
        )

    return list(set(genbank_accessions))  # prevent quering the same accession multiple times


def get_cazy_class_fam_genbank_records(args, session, config_dict):
    """GenBank acc query results from the local CAZyme database for CAZyme from specific classes/fams

    :param args: cmd-line argument parser
    :param session: open SQLite db session
    :param config_dict: dict, defines CAZy classes and families to get sequences for

    Return CAZy class and CAZy family GenBank accession query results
    """
    logger = logging.getLogger(__name__)
    if args.update:  # retrieve all GenBank accessions
        if args.primary:
            logger.warning(
                "Retrieving sequences for PRIMARY GenBank accessions that:\n"
                "belong to specific CAZy classes/families AND\n"
                "do not have a sequence in the db OR the sequence has been updated in NCBI"
            )
            (
                genbank_query_class,
                genbank_query_family,
            ) = query_sql_db.get_prim_gnbk_acc_from_clss_fams(
                session,
                config_dict,
            )

        else:
            logger.warning(
                "Retrieving sequences for PRIMARY GenBank accessions that:\n"
                "belong to specific CAZy classes/families AND\n"
                "do not have a sequence in the db OR the sequence has been updated in NCBI"
            )
            (
                genbank_query_class,
                genbank_query_family,
            ) = query_sql_db.get_all_gnbk_acc_from_clss_fams(
                session,
                config_dict,
            )

    else:  # retrieve GenBank accesions of records that don't have a sequence
        if args.primary:
            logger.warning(
                "Retrieving sequences for PRIMARY GenBank accessions that:\n"
                "belong to specific CAZy classes/families AND do not have a sequence in the db"
            )
            (
                genbank_query_class,
                genbank_query_family,
            ) = query_sql_db.get_prim_gnbk_acc_from_clss_fams_no_seq(
                session,
                config_dict,
            )

        else:
            logger.warning(
                "Retrieving sequences for PRIMARY GenBank accessions that:\n"
                "belong to specific CAZy classes/families AND do not have a sequence in the db"
            )
            (
                genbank_query_class,
                genbank_query_family,
            ) = query_sql_db.get_all_gnbk_acc_from_clss_fams_no_seq(
                session,
                config_dict,
            )

    return genbank_query_class, genbank_query_family


def parse_genbank_query(genbank_query, taxonomy_filters, kingdoms, ec_filters, session):
    """Parse SQL query result and retrieve GenBank accessions of CAZymes that meet the user cirteria

    :param:

    Return list of GenBank accessions.
    """
    if genbank_query is None:
        return []

    if (taxonomy_filters is None) and (kingdoms is None) and (ec_filters is None):
        accessions = [item[0] for item in genbank_query]
        return [x for x in accessions if "NA" != x]
    
    genbank_accessions = []

    for item in genbank_query:
        if item[0] != "NA":  # if GenBank accession not stored as 'NA'

            # check if CAZyme records meets the taxonomy criteria
            source_organism = item[-2].genus + item[-2].species
            if any(filter in source_organism for filter in taxonomy_filters):
                genbank_accessions.append(item[0])
                continue

            # check if CAZyme records meets the kingdom requirement
            if item[-1].kingdom in kingdoms:
                genbank_accessions.append(item[0])
                continue

    if ec_filters is None:
        return genbank_accessions
    
    # check if the parent CAZymes of the GenBank accessions meet the EC number filter
    filtered_genbank_accessions = []
    for i in tqdm(range(len(genbank_accessions), desc="Applying EC number filter")):
        ec_annotations = query_sql_db.query_sql_db.query_ec_number(session, genbank_accessions[i])
        if (set(ec_annotations) and set(ec_filters)):
            filtered_genbank_accessions.append(genbank_accessions[i])

    return filtered_genbank_accessions
