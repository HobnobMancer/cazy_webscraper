#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
"""Retrieve proteins sequences from GenBank and populate the local database and write to FASTA"""


import logging
import os
import time

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio import Entrez, SeqIO
from tqdm import tqdm

from scraper import file_io
from scraper.sql.sql_orm import (
    Cazyme,
    CazyFamily,
    Cazymes_Genbanks,
    Genbank,
    get_db_session,
)
from scraper.utilities import build_logger, build_genbank_sequences_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # parse cmd-line arguments
    if argv is None:
        parser = build_genbank_sequences_parser()
        args = parser.parse_args()
    else:
        args = build_genbank_sequences_parser(argv).parse_args()

    # build logger
    logger = build_logger("expand.genbank_sequences", args)

    # check database was passed
    if os.path.isfile(args.database) is False:
        logger.error(
            "Could not find local CAZy database. Check pack is correct. Terminating programme."
        )

    Entrez.email = args.email

    # create session to local database
    session = get_db_session(args, logger)

    # check if any classes or families were specified to retrieve the sequences only for them
    if (args.classes is None) and (args.families is None):
        config_dict = None
    
    else:
        # create dictionary of CAZy classes/families to retrieve sequences for
        file_io_path = file_io.__file__
        cazy_dict, std_class_names = file_io.get_cazy_dict_std_names(file_io_path, logger)
        config_dict = file_io.get_cmd_defined_fams_classes(cazy_dict, std_class_names, args, logger)

    if config_dict is None:
        # get sequences for everything
        get_everything_sequences(session, args, logger)

    else:
        # get sequences for only specified classes/families
        if args.primary is True:
            get_specific_proteins_sequencse_primary_only(config_dict, session, args, logger)
        else:
            get_specific_proteins_sequencse(config_dict, session, args, logger)

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


def get_everything_sequences(session, args, logger):
    """Retrieve protein sequences for all CAZymes in the local CAZy database.

    :param session: open SQLite db session
    :param args: cmd-line argument parser
    :param logger: logger object

    Return nothing.
    """
    # retrieve only sequences for primary GenBank accessions
    if args.primary is True:
        # retrieve all primary GenBank accessions
        genbank_query = session.query(Genbank, Cazymes_Genbanks).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            filter(Cazymes_Genbanks.primary == True).\
            all()

    # retrieve sequences for all GenBank accessions
    else:
        # retrieve all GenBank accessions
        genbank_query = session.query(Genbank, Cazymes_Genbanks).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            all()

    for genbank_object in tqdm(genbank_query, desc="Retrieving sequences"):
        genbank_accession = genbank_object[0].genbank_accession

        # check if a protein sequence is already stored in the local database
        if genbank_object[0].sequence is None:
            get_new_protein_sequence(genbank_accession, genbank_object[0], session, args, logger)

        elif (genbank_object[0].sequence is not None) and (args.update is True):
            update_protein_sequence(genbank_accession, genbank_object[0], session, args, logger)

    return


def get_specific_proteins_sequencse_primary_only(config_dict, session, args, logger):
    """Retrieve protein sequences for primary GenBank accessions for items in the config_dict.

    :param config_dict: dict, defines CAZy classes and families to retrieve accessions from
    :param session: open SQLite db session
    :param args: cmd-line argument parser
    :param logger: logger object

    Return nothing.
    """
    # start with the classes
    if len(config_dict["classes"]) != 0:
        # create a dictionary to convert full class name to abbreviation
        cazy_classes = config_dict["classes"]
        for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):
            # retrieve class name abbreviation
            cazy_class = cazy_class[cazy_class.find("("):cazy_class.find(")")]

            # retrieve all GenBank accessions catalogued in the CAZy class
            class_query = session.query(CazyFamily, Cazyme, Cazymes_Genbanks, Genbank).\
                join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                join(Cazyme, (Cazyme.cazyme_id == Cazymes_Genbanks.cazyme_id)).\
                filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
                filter(Cazymes_Genbanks.primary==True).\
                all()

            for query_result in tqdm(class_query, desc=f"Parsing families in {cazy_class}"):
                genbank_accession = query_result[-1].genbank_accession

                # check if a protein sequence is already stored in the local database
                if query_result[-1].sequence is None:
                    get_new_protein_sequence(
                        genbank_accession,
                        query_result[-1],
                        session,
                        args,
                        logger,
                    )

                elif (query_result[-1].sequence is not None) and (args.update is True):
                    update_protein_sequence(
                        genbank_accession,
                        query_result[-1],
                        session,
                        args,
                        logger,
                    )

    # Retrieve protein sequences for specified families
    for key in config_dict:
        if key == "classes":
            continue

        for family in tqdm(config_dict[key], desc=f"Parsing families in {key}"):
            if family.find("_") != -1:  # subfamily
                # Retrieve GenBank accessions catalogued under the subfamily
                family_query = session.query(CazyFamily, Cazyme, Cazymes_Genbanks, Genbank).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    join(Cazyme, (Cazyme.cazyme_id == Cazymes_Genbanks.cazyme_id)).\
                    filter(CazyFamily.subfamily == family).\
                    filter(Cazymes_Genbanks.primary==True).\
                    all()

            else:  # family
                # Retrieve GenBank accessions catalogued under the family
                family_query = session.query(CazyFamily, Cazyme, Cazymes_Genbanks, Genbank).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    join(Cazyme, (Cazyme.cazyme_id == Cazymes_Genbanks.cazyme_id)).\
                    filter(CazyFamily.family == family).\
                    filter(Cazymes_Genbanks.primary == True).\
                    all()

            for query_result in tqdm(family_query, desc=f"Parsing GenBank accessions in {family}"):
                genbank_accession = query_result[-1].genbank_accession

                if query_result[-1].sequence is None:
                    get_new_protein_sequence(
                        genbank_accession,
                        query_result[-1],
                        session,
                        args,
                        logger,
                    )

                elif (query_result[-1].sequence is not None) and (args.update is True):
                    update_protein_sequence(
                        genbank_accession,
                        query_result[-1],
                        session,
                        args,
                        logger,
                    )

    return


def get_specific_proteins_sequencse(config_dict, session, args, logger):
    """Retrieve protein sequences for only CAZymes specified in config_dict.+

    :param config_dict: dict, defines CAZy classes and families to retrieve accessions from
    :param session: open SQLite db session
    :param args: cmd-line argument parser
    :param logger: logger object

    Return nothing.
    """
    # start with the classes
    if len(config_dict["classes"]) != 0:
        # create a dictionary to convert full class name to abbreviation
        cazy_classes = config_dict["classes"]
        for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):
            # retrieve class name abbreviation
            cazy_class = cazy_class[cazy_class.find("("):cazy_class.find(")")]

            # retrieve all GenBank accessions catalogued in the CAZy class
            class_query = session.query(CazyFamily, Cazyme, Cazymes_Genbanks, Genbank).\
                join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                join(Cazyme, (Cazyme.cazyme_id == Cazymes_Genbanks.cazyme_id)).\
                filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
                all()

            for query_result in tqdm(class_query, desc=f"Parsing families in {cazy_class}"):
                genbank_accession = query_result[-1].genbank_accession

                # check if a protein sequence is already stored in the local database
                if query_result[-1].sequence is None:
                    get_new_protein_sequence(
                        genbank_accession,
                        query_result[-1],
                        session,
                        args,
                        logger,
                    )

                elif (query_result[-1].sequence is not None) and (args.update is True):
                    update_protein_sequence(
                        genbank_accession,
                        query_result[-1],
                        session,
                        args,
                        logger,
                    )

    # Retrieve protein sequences for specified families
    for key in config_dict:
        if key == "classes":
            continue

        for family in tqdm(config_dict[key], desc=f"Parsing families in {key}"):
            if family.find("_") != -1:  # subfamily
                # Retrieve GenBank accessions catalogued under the subfamily
                family_query = session.query(CazyFamily, Cazyme, Cazymes_Genbanks, Genbank).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    join(Cazyme, (Cazyme.cazyme_id == Cazymes_Genbanks.cazyme_id)).\
                    filter(CazyFamily.subfamily == family).\
                    all()

            else:  # family
                # Retrieve GenBank accessions catalogued under the family
                family_query = session.query(CazyFamily, Cazyme, Cazymes_Genbanks, Genbank).\
                    join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
                    join(Cazyme, (Cazyme.cazyme_id == Cazymes_Genbanks.cazyme_id)).\
                    filter(CazyFamily.family == family).\
                    all()

            for query_result in tqdm(family_query, desc=f"Parsing GenBank accessions in {family}"):
                if query_result[-1].sequence is None:
                    get_new_protein_sequence(
                        genbank_accession,
                        query_result[-1],
                        session,
                        args,
                        logger,
                    )

                elif (query_result[-1].sequence is not None) and (args.update is True):
                    update_protein_sequence(
                        genbank_accession,
                        query_result[-1],
                        session,
                        args,
                        logger,
                    )

    return


def get_new_protein_sequence(genbank_accession, genbank_record, session, args, logger):
    """Retrieve protein sequence from NCBI.

    :param genbank_accession: str, GenBank protein accession
    :param genbank_record: Genbank instance, existing Genbank record in the local database
    :param args: cmd-line argument parser
    :param logger: logger object

    Return nothing.
    """
    # retrieve protein record from NCBI Protein database
    handle = entrez_retry(
        logger,
        Entrez.efetch,
        db="Protein",
        id=genbank_accession,
        rettype="fasta",
        retmode="text",
    )

    if handle is None:
        logger.warning(f"Failed to retrieve protein sequence for {genbank_accession}")
        return

    # add sequence to the local database
    for record in SeqIO.parse(handle, "fasta"):
        retrieved_sequence = str(record.seq)  # convert to a string becuase SQL expects a string
        genbank_record.sequence = retrieved_sequence
        session.commit()

    # write out sequence to FASTA file if enabled
    if args.write is not None:
        write_out_fasta(handle, genbank_accession, args)

    return


def update_protein_sequence(genbank_accession, genbank_record, session, args, logger):
    """Retrieve protein sequence from NCBI and update local database IF the sequence is different.

    :param genbank_accession: str, GenBank protein accession
    :param genbank_record: Genbank instance, existing Genbank record in the local database
    :param args: cmd-line argument parser
    :param logger: logger object

    Return nothing.
    """
    # retrieve protein record from NCBI Protein database
    handle = entrez_retry(
        logger,
        Entrez.efetch,
        db="Protein",
        id=genbank_accession,
        rettype="fasta",
        retmode="text",
    )

    if handle is None:
        logger.warning(f"Failed to retrieve protein sequence for {genbank_accession}")
        return

    # check if the local and retrieved sequences are the same, update sequence if different
    existing_sequence = genbank_record.sequence

    for record in SeqIO.parse(handle, "fasta"):
        retrieved_sequence = record.seq
        if retrieved_sequence != existing_sequence:
            # update the local record
            genbank_record.sequence = retrieved_sequence
            session.commit()

    # write out sequence to FASTA file if enabled
    if args.write is not None:
        write_out_fasta(handle, genbank_accession, args)

    return


def entrez_retry(logger, entrez_func, *func_args, **func_kwargs):
    """Call to NCBI using Entrez.

    Maximum number of retries is 10, retry initated when network error encountered.

    :param logger: logger object
    :param retries: parser argument, maximum number of retries excepted if network error encountered
    :param entrez_func: function, call method to NCBI

    :param *func_args: tuple, arguments passed to Entrez function
    :param ** func_kwargs: dictionary, keyword arguments passed to Entrez function

    Returns record.
    """
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


def write_out_fasta(handle, genbank_accession, args):
    """Write out GenBank protein record to a FASTA file.

    :param handle: Entrez record instance
    :param genbank_accession: str, GenBank protein accession number
    :param args: cmd-line arguments parser

    Return nothing.
    """
    for record in SeqIO.parse(handle, "fasta"):
        fasta_name = f"{genbank_accession}.fasta"
        fasta_name = args.write / fasta_name
        with open(fasta_name, "w") as fh:
            SeqIO.write(record, fh, "fasta")
    return

if __name__ == "__main__":
    main()
