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
"""Retrieve proteins sequences from GenBank and populate the local database and write to FASTA"""


import logging
import os
import re
import sys
import time

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio import Entrez, SeqIO
from tqdm import tqdm

from scraper.sql.sql_orm import (
    Cazyme,
    CazyFamily,
    Cazymes_Genbanks,
    Genbank,
    Kingdom,
    Taxonomy,
    get_db_session,
)
from scraper.utilities import config_logger, file_io, parse_configuration
from scraper.utilities.parsers import build_genbank_sequences_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)
    date_today = datetime.now().strftime("%Y/%m/%d")  # used as seq_update_date in the db

    # parse cmd-line arguments
    if argv is None:
        parser = build_genbank_sequences_parser()
        args = parser.parse_args()
    else:
        args = build_genbank_sequences_parser(argv).parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)

    # check database was passed
    if os.path.isfile(args.database) is False:
        logger.error(
            "Could not find local CAZy database. Check path is correct. Terminating programme."
        )
        sys.exit(1)

    Entrez.email = args.email

    # create session to local database
    session = get_db_session(args)

    # retrieve configuration data
    file_io_path = file_io.__file__
    config_dict, taxonomy_filters, kingdoms = parse_configuration.get_configuration(
        file_io_path,
        args,
    )

    if config_dict is None:
        if args.update:
            # get sequence for everything without a sequence and those with newer remote sequence
            add_and_update_all_sequences(date_today, taxonomy_filters, kingdoms, session, args)

        else:
            # get sequences for everything without a sequence
            get_missing_sequences_for_everything(
                date_today,
                taxonomy_filters,
                kingdoms,
                session,
                args,
            )

    else:
        # get sequences for only specified classes/families
        if args.update:
            update_sequences_for_specific_records(
                date_today,
                config_dict,
                taxonomy_filters,
                kingdoms,
                session,
                args,
            )
        else:
            get_missing_sequences_for_specific_records(
                date_today,
                config_dict,
                taxonomy_filters,
                kingdoms,
                session,
                args,
            )

    if args.blastdb is not None:
        file_io.build_blast_db(args)

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    end_time = pd.to_datetime(start_time)
    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    logger.info(
        "Finished populating local CAZy database with GenBank protein sequences. "
        "Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}"
    )

    print(
        "=====================cazy_webscraper-expand-genank_sequences=====================\n"
        "Finished populating local CAZy database with GenBank protein sequences. "
        "Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}\n"
    )


# The folowing functions are for querying the local database to get GenBank accessions


def get_missing_sequences_for_everything(date_today, taxonomy_filters, kingdoms, session, args):
    """Retrieve protein sequences for all CAZymes in the local CAZy database that don't have seq.

    :param date_today: str, today's date, used for logging the date the seq is retrieved in the db
    :param taxonomy_filters: set of genera, species and strains to restrict sequence retrieval
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param session: open SQLite db session
    :param args: cmd-line argument parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # retrieve only sequences for primary GenBank accessions, and those without sequences
    if args.primary is True:
        logger.warning(
            "Retrieving sequences for all primary GenBank accessions that do not have sequences"
        )
        genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
            join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
            join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
            join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            filter(Cazymes_Genbanks.primary == True).\
            filter(Genbank.sequence == None).\
            all()

    # retrieve sequences for all GenBank accessions without sequences
    else:
        logger.warning(
            "Retrieving sequences for all GenBank accessions that do not have sequences"
        )
        genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
            join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
            join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
            join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            filter(Genbank.sequence == None).\
            all()

    # retrieve the genbank_accessions
    accessions = extract_accessions(genbank_query, taxonomy_filters)

    if len(accessions) == 0:
        logger.warning(
            "Did not retrieve any GenBank accessions from the local database\n"
            "that have sequences missing. Not adding sequences to the local database."
        )
        return

    # separate accesions in to separate lists of length args.epost, epost doesn't like more than 200
    accessions = get_accession_chunks(accessions, args.epost)  # args.epost = number per chunk
    for lst in accessions:
        get_sequences_add_to_db(lst, date_today, session, args)
    return


def add_and_update_all_sequences(date_today, taxonomy_filters, kingdoms, session, args):
    """Retrieve sequences for all proteins in the database.

    For records with no sequences, add the retrieved sequence.
    For records with a sequence, check if the remove sequence is more recent than the existing
    sequence. It it is, update the local sequence.

    :param date_today: str, today's date, used for logging the date the seq is retrieved in the db
    :param taxonomy_filters: set of genera, species and strains to retrieve sequences for
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param session: open SQLite db session
    :param args: cmd-line argument parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # retrieve only sequences for primary GenBank accessions, and those without sequences
    if args.primary is True:
        logger.warning(
            "Retrieving sequences for all primary GenBank accessions that do not have sequences\n"
            "and those whose sequences have been updated in NCBI "
            "since they were retrieved previously"
        )
        genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
            join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
            join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
            join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            filter(Cazymes_Genbanks.primary == True).\
            all()

    # retrieve sequences for all GenBank accessions
    else:
        logger.warning(
            "Retrieving sequences for all GenBank accessions that do not have sequences\n"
            "and those whose sequences have been updated in NCBI "
            "since they were retrieved previously"
        )
        genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
            join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
            join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
            join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            all()

    # create dictionary of {genbank_accession: 'sequence update date' (str)}
    accessions = extract_accessions_and_dates(genbank_query, taxonomy_filters)

    if len(accessions.keys()) == 0:
        logger.warning(
            "Did not retrieve any GenBank accessions from the local database.\n"
            "Not adding sequences to the local database."
        )
        return

    accessions = get_accessions_for_new_sequences(accessions)  # list of genkbank_accession

    if len(accessions) == 0:
        logger.warning(
            "Did not retrieve any GenBank accessions whose sequences need updating.\n"
            "Not adding sequences to the local database."
        )
        return

    # separate accesions in to separate lists of length args.epost, epost doesn't like more than 200
    accessions = get_accession_chunks(accessions, args.epost)  # args.epost = number per chunk
    for lst in accessions:
        get_sequences_add_to_db(lst, date_today, session, args)
    return


def get_missing_sequences_for_specific_records(
    date_today,
    config_dict,
    taxonomy_filters,
    kingdoms,
    session,
    args,
):
    """Coordinate getting the sequences for specific CAZymes, not with seqs in the db.

    :param date_today: str, today's date, used for logging the date the seq is retrieved in the db
    :param config_dict: dict, defines CAZy classes and families to get sequences for
    :param taxonomy_filters: set of genera, species and strains to restrict sequence retrieval
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param session: open SQL database session
    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    logger.warning(
        "Retrieving sequences for GenBank accessions that do not have a sequence in the database"
    )

    # start with the classes
    if len(config_dict["classes"]) != 0:
        # retrieve list of CAZy classes to get sequences for
        cazy_classes = config_dict["classes"]

        for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):
            # retrieve class name abbreviation
            cazy_class = cazy_class[((cazy_class.find("(")) + 1):((cazy_class.find(")")) - 1)]

            # get the CAZymes within the CAZy class
            class_subquery = session.query(Cazyme.cazyme_id).\
                join(CazyFamily, Cazyme.families).\
                filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
                subquery()

            # retrieve the GenBank accessions of the CAZymes in the CAZy class without seqs
            if args.primary:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(class_subquery)).\
                    filter(Cazymes_Genbanks.primary == True).\
                    filter(Genbank.sequence == None).\
                    all()
            else:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(class_subquery)).\
                    filter(Genbank.sequence == None).\
                    all()

            # retrieve the genbank_accessions from the sql collection object returned by the query
            accessions = extract_accessions(genbank_query, taxonomy_filters)

            if len(accessions) == 0:
                logger.warning(
                    f"Did not retrieve any GenBank accessions for the CAZy class {cazy_class}\n"
                    "that have missing sequences. Not adding sequences to the local database."
                )
                continue

            # separate accesions in to separate lists of length args.epost
            # epost doesn't like posting more than 200 at once
            accessions = get_accession_chunks(accessions, args.epost)  # args.epost = number/chunk
            for lst in accessions:
                get_sequences_add_to_db(lst, date_today, session, args)
            continue

    # Retrieve protein sequences for specified families
    for key in config_dict:
        if key == "classes":
            continue
        if config_dict[key] is None:
            continue  # no families to parse

        for family in tqdm(config_dict[key], desc=f"Parsing families in {key}"):
            if family.find("_") != -1:  # subfamily
                # Retrieve GenBank accessions catalogued under the subfamily
                family_subquery = session.query(Cazyme.cazyme_id).\
                    join(CazyFamily, Cazyme.families).\
                    filter(CazyFamily.subfamily == family).\
                    subquery()

            else:  # family
                # Retrieve GenBank accessions catalogued under the family
                family_subquery = session.query(Cazyme.cazyme_id).\
                    join(CazyFamily, Cazyme.families).\
                    filter(CazyFamily.family == family).\
                    subquery()

            # get the GenBank accessions of thes CAZymes, without sequences
            if args.primary:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(family_subquery)).\
                    filter(Cazymes_Genbanks.primary == True).\
                    filter(Genbank.sequence == None).\
                    all()
            else:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(family_subquery)).\
                    filter(Genbank.sequence == None).\
                    all()

            # retrieve a list of GenBank accessions from the sql collection returned from the query
            accessions = extract_accessions(genbank_query, taxonomy_filters)

            if len(accessions) == 0:
                logger.warning(
                    f"Did not retrieve any GenBank accessions for the CAZy class {family}\n"
                    "that have missing sequences. Not adding sequences to the local database."
                )
                continue
            # separate accesions in to separate lists of length args.epost
            # epost doesn't like posting more than 200 at once
            accessions = get_accession_chunks(accessions, args.epost)  # args.epost = acc/chunk
            for lst in accessions:
                get_sequences_add_to_db(lst, date_today, session, args)

    return


def update_sequences_for_specific_records(
    date_today,
    config_dict,
    taxonomy_filters,
    kingdoms,
    session,
    args,
):
    """Coordinate getting the sequences for specific CAZymes, not with seqs in the db nad those
    whose seq in NCBI has been updated since the last retrieval.

    For records with no sequences, add the retrieved sequence.
    For records with a sequence, check if the remove sequence is more recent than the existing
    sequence. It it is, update the local sequence.

    :param date_today: str, today's date, used for logging the date the seq is retrieved in the db
    :param config_dict: dict, defines CAZy classes and families to get sequences for
    :param taxonomy_filters: set of genera, species and strains to restrict sequence retrieval
    :param kingdoms: set of taxonomy Kingdoms to retrieve sequences for
    :param session: open SQL database session
    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    logger.warning(
        "Retrieving sequences for GenBank accessions that do not have a sequence in the database,\n"
        "and those whose sequence in NCBI has been updated since they were previously retrieved."
    )

    # start with the classes
    if len(config_dict["classes"]) != 0:
        # retrieve list of CAZy classes to get sequences for
        cazy_classes = config_dict["classes"]

        for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):
            # retrieve class name abbreviation
            cazy_class = cazy_class[((cazy_class.find("(")) + 1):((cazy_class.find(")")) - 1)]

            # get the CAZymes within the CAZy class
            class_subquery = session.query(Cazyme.cazyme_id).\
                join(CazyFamily, Cazyme.families).\
                filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
                subquery()

            # retrieve the GenBank accessions of the CAZymes in the CAZy class without seqs
            if args.primary:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(class_subquery)).\
                    filter(Cazymes_Genbanks.primary == True).\
                    all()
            else:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(class_subquery)).\
                    all()

            # create dictionary of genbank_accession: 'sequence update date' (str)
            accessions = extract_accessions_and_dates(genbank_query, taxonomy_filters)

            if len(accessions.keys()) == 0:
                logger.warning(
                    f"Did not retrieve any GenBank accessions for the CAZy class {cazy_class}.\n"
                    "Not adding sequences to the local database."
                )
                continue

            accessions = get_accessions_for_new_sequences(accessions)  # list of genkbank_accession

            if len(accessions) == 0:
                logger.warning(
                    "Did not retrieve any GenBank accessions whose sequences need updating for "
                    f"the CAZy class {cazy_class}.\n"
                    "Not adding sequences to the local database."
                )
                continue
            # separate accesions in to separate lists of length args.epost
            # epost doesn't like posting more than 200 at once
            accessions = get_accession_chunks(accessions, args.epost)  # args.epost = acc/chunk
            for lst in accessions:
                get_sequences_add_to_db(lst, date_today, session, args)

    # Retrieve protein sequences for specified families
    for key in config_dict:
        if key == "classes":
            continue
        if config_dict[key] is None:
            continue  # no families to parse

        for family in tqdm(config_dict[key], desc=f"Parsing families in {key}"):
            if family.find("_") != -1:  # subfamily
                # Retrieve GenBank accessions catalogued under the subfamily
                family_subquery = session.query(Cazyme.cazyme_id).\
                    join(CazyFamily, Cazyme.families).\
                    filter(CazyFamily.subfamily == family).\
                    subquery()

            else:  # family
                # Retrieve GenBank accessions catalogued under the family
                family_subquery = session.query(Cazyme.cazyme_id).\
                    join(CazyFamily, Cazyme.families).\
                    filter(CazyFamily.family == family).\
                    subquery()

            # get the GenBank accessions of thes CAZymes, without sequences
            if args.primary:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(family_subquery)).\
                    filter(Cazymes_Genbanks.primary == True).\
                    filter(Genbank.sequence == None).\
                    all()
            else:
                genbank_query = session.query(Genbank, Cazymes_Genbanks, Cazyme, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(Cazymes_Genbanks, (Cazymes_Genbanks.cazyme_id == Cazyme.cazyme_id)).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    filter(Cazyme.cazyme_id.in_(family_subquery)).\
                    filter(Genbank.sequence == None).\
                    all()

            # create dictionary of {genbank_accession: 'sequence update date' (str)}
            accessions = extract_accessions_and_dates(genbank_query, taxonomy_filters)

            if len(accessions.keys()) == 0:
                logger.warning(
                    f"Did not retrieve any GenBank accessions for the CAZy class {family}.\n"
                    "Not adding sequences to the local database."
                )
                continue

            accessions = get_accessions_for_new_sequences(accessions)  # list of genkbank_accession

            if len(accessions) == 0:
                logger.warning(
                    "Did not retrieve any GenBank accessions whose sequences need updating for "
                    f"the CAZy class {family}.\n"
                    "Not adding sequences to the local database."
                )
                continue
            # separate accesions in to separate lists of length args.epost
            # epost doesn't like posting more than 200 at once
            accessions = get_accession_chunks(accessions, args.epost)  # args.epost = acc/chunk
            for lst in accessions:
                get_sequences_add_to_db(lst, date_today, session, args)

    return


# The following functions are retrieving the list of Genbank accessions to retrieve sequences for #


def extract_accessions(genbank_query, taxonomy_filters):
    """The query contains GenBank accessions and Cazymes_Genbanks records, retrieve the accessions.

    :param genbank_query: sql collection
    :param taxonomy_filters: set of genera, species and strains to restrict retrieval of sequences

    Return a list of GenBank accessions. Each element is a string of a unique accession.
    """
    if taxonomy_filters is None:
        accessions = [item[0] for item in genbank_query]
        return [x for x in accessions if "NA" != x]

    else:
        accessions = []
        for item in genbank_query:
            if item[0] != "NA":  # if GenBank accession not stored as 'NA'
                source_organism = item[-1].genus + item[-1].species
                if any(filter in source_organism for filter in taxonomy_filters):
                    accessions.append(item[0])
    return accessions


def extract_accessions_and_dates(genbank_query, taxonomy_filters):
    """Retrieve the GenBank accessions and retrieval dates of existing sequences from the db query.

    :param genbank_query: sql collection
    :param taxonomy_filters: set of genera, species and strains to restrict retrieval of sequences

    Return a dict {GenBank_accession: retrieval_date}
    """
    accessions = {}

    if taxonomy_filters is None:
        for item in genbank_query:
            if item[0].genbank_accession == "NA":  # no GenBank accession stored in CAZy
                continue
            accessions[item[0].genbank_accession] = item[0].seq_update_date

    else:
        for item in genbank_query:
            if item[0].genbank_accession == "NA":  # no GenBank accession stored in CAZy
                continue
            source_organism = item[-1].genus + item[-1].species
            if any(filter in source_organism for filter in taxonomy_filters):
                accessions[item[0].genbank_accession] = item[0].seq_update_date

    return accessions


def get_accessions_for_new_sequences(accessions):
    """Get the GenBank accessions of sequences to be added to the local database.

    For records currently with no protein sequence, the retrieved protein sequence will be added
    to the record. For records with a sequence, the 'UpdateDate' for the sequence from NCBI will
    be compared against the  'seq_update_date' in the local database. The 'seq_update_date' is the
    'UpdateDate' previosuly retrieved from NCBI. If the NCBI sequence is newer,
    the local database will be updated with the new sequence.

    :param accessions: dict, {GenBank accessions (str):sequence retrieval data (str)}
    :param session: open SQL database session

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    accessions_list = list(accessions.keys())
    accessions_string = ",".join(accessions_list)
    # perform batch query of Entrez
    epost_result = Entrez.read(
        entrez_retry(
            Entrez.epost, "Protein", id=accessions_string, retmode="text",
        )
    )
    # retrieve the web environment and query key from the Entrez post
    epost_webenv = epost_result["WebEnv"]
    epost_query_key = epost_result["QueryKey"]

    # retrieve summary docs to check the sequence 'UpdateDates' in NCBI
    with entrez_retry(
        Entrez.efetch,
        db="Protein",
        query_key=epost_query_key,
        WebEnv=epost_webenv,
        rettype="docsum",
        retmode="xml",
    ) as handle:
        summary_docs = Entrez.read(handle)

    for doc in summary_docs:
        try:
            temp_accession = doc["AccessionVersion"]  # accession of the current working protein
        except KeyError:
            logger.warning(
                f"Retrieved protein with accession {temp_accession} but this accession is not in "
                "the local database.\n"
                "Not retrieving a sequence for this accession."
            )
            continue
        previous_data = accessions[temp_accession]
        if previous_data is not None:
            # sequence retrieved previosuly, thus check if the NCBI seq has been updated since
            previous_data = previous_data.split("/")  # Y=[0], M=[1], D=[]
            update_date = doc["UpdateDate"]
            update_date = update_date.split("/")  # Y=[0], M=[1], D=[]
            if datetime.date(
                previous_data[0], previous_data[1], previous_data[2],
            ) < datetime.data(
                update_date[0], update_date[1], update_date[2],
            ) is False:
                # the sequence at NCBI has not been updated since the seq was retrieved
                # thus no need to retrieve it again
                accessions_list.remove(temp_accession)

    return accessions_list


def get_accession_chunks(lst, chunk_length):
    """Separate the long list into separate chunks.

    :param lst: list to be separated into smaller lists (or chunks)
    :param chunk_length: int, the length of the lists the longer list is to be split up into

    Return a generator object containing lists.
    """
    for i in range(0, len(lst), chunk_length):
        yield lst[i:i + chunk_length]


# The following functions are for retrieving sequences, adding to the db and writing fasta files


def get_sequences_add_to_db(accessions, date_today, session, args):
    """Retrieve protein sequences from Entrez and add to the local database.

    :param accessions: list, GenBank accessions
    :param date_today: str, YYYY/MM/DD
    :param session: open SQL database session
    :param args: cmb-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    # perform batch query of Entrez
    accessions_string = ",".join(accessions)
    epost_result = Entrez.read(
        entrez_retry(
            Entrez.epost, "Protein", id=accessions_string,
        )
    )
    # retrieve the web environment and query key from the Entrez post
    epost_webenv = epost_result["WebEnv"]
    epost_query_key = epost_result["QueryKey"]

    # retrieve the protein sequences
    with entrez_retry(
        Entrez.efetch,
        db="Protein",
        query_key=epost_query_key,
        WebEnv=epost_webenv,
        rettype="fasta",
        retmode="text",
    ) as seq_handle:
        for record in SeqIO.parse(seq_handle, "fasta"):
            # retrieve the accession of the record
            temp_accession = record.id  # accession of the current working protein record

            if temp_accession.find("|") != -1:  # sometimes multiple items are listed
                success = False   # will be true if finds protein accession
                temp_accession = temp_accession.split("|")
                for item in temp_accession:
                    # check if a accession number
                    try:
                        re.match(
                            (
                                r"(\D{3}\d{5,7}\.\d+)|(\D\d(\D|\d){3}\d)|"
                                r"(\D\d(\D|\d){3}\d\D(\D|\d){2}\d)"
                            ),
                            item,
                        ).group()
                        temp_accession = item
                        success = True
                        break
                    except AttributeError:  # raised if not an accession
                        continue

            else:
                success = True  # have protein accession number

            if success is False:
                logger.error(
                    f"Could not retrieve accession from {record.id}, therefore, "
                    "protein sequence not added to the database,\n"
                    "because cannot retrieve the necessary CAZyme record"
                )
                continue

            # check the retrieve protein accession is in the list of retrieved accession
            if temp_accession not in accessions:
                logger.warning(
                    f"Retrieved the accession {temp_accession} from the record id={record.id}, "
                    "but this accession is not in the database.\n"
                    "Therefore, not adding this protein seqence to the local database"
                )
                continue

            # retrieve the GenBank record from the local data base to add the seq to
            genbank_record = session.query(Genbank).\
                filter(Genbank.genbank_accession == temp_accession).first()

            retrieved_sequence = str(record.seq)  # convert to a string becuase SQL expects a string
            genbank_record.sequence = retrieved_sequence
            genbank_record.seq_update_date = date_today
            session.commit()

            if args.fasta is not None:
                file_io.write_out_fasta(record, temp_accession, args)

            if args.blastdb is not None:
                file_io.write_fasta_for_db(record, temp_accession, args)

            # remove the accession from the list
            accessions.remove(temp_accession)

    if len(accessions) != 0:
        logger.warning(
            "Protein sequences were not retrieved for the following CAZyme in the local database"
        )
        for acc in accessions:
            logger.warning(f"GenBank accession: {acc}")

    return


def entrez_retry(entrez_func, *func_args, **func_kwargs):
    """Call to NCBI using Entrez.

    Maximum number of retries is 10, retry initated when network error encountered.

    :param logger: logger object
    :param retries: parser argument, maximum number of retries excepted if network error encountered
    :param entrez_func: function, call method to NCBI

    :param *func_args: tuple, arguments passed to Entrez function
    :param ** func_kwargs: dictionary, keyword arguments passed to Entrez function

    Returns record.
    """
    logger = logging.getLogger(__name__)
    record, retries, tries = None, 10, 0

    while record is None and tries < retries:
        try:
            record = entrez_func(*func_args, **func_kwargs)

        except IOError:
            # log retry attempt
            if tries < retries:
                logger.warning(
                    f"Network error encountered during try no.{tries}.\nRetrying in 10s",
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered too many times. Exiting attempt to call to NCBI"
        )
        return

    return record


if __name__ == "__main__":
    main()
