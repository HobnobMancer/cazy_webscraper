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
"""Submodule for handling inputs and outputs"""


import logging
import shutil
import sys

from Bio import SeqIO


def make_output_directory(output, force, nodelete):
    """Create output directory for genomic files.

    :param output: path, path of dir to be created
    :param force: bool, enable/disable creating dir if already exists
    :param nodelete: bool, enable/disable deleting content in existing dir

    Raises FileExistsError if an attempt is made to create a directory that already
    exist and force (force overwrite) is False.

    Return Nothing
    """
    logger = logging.getLogger(__name__)

    if output.exists():
        if force is True:

            if nodelete is True:
                logger.warning(
                    f"Output directory {output} exists, nodelete is {nodelete}. "
                    "Adding output to output directory."
                )
                return

            else:
                logger.warning(
                    f"Output directory {output} exists, nodelete is {nodelete}. "
                    "Deleting content currently in output directory."
                )
                shutil.rmtree(output)
                output.mkdir(exist_ok=force)
                return

        else:
            logger.warning(
                f"Output directory {output} exists. 'force' is False, cannot write to existing "
                "output directory.\nTerminating program."
            )
            sys.exit(1)

    else:
        output.mkdir(exist_ok=force)
        logger.warning(f"Built output directory: {output}")

    return


def write_out_failed_scrapes(failed_urls, time_stamp, args):
    """Write out the URLs for which a connection to CAZy failed.
    :param failed_urls: list, contains the URL and reason for the failed scrape
    :param args: cmd args parser
    Return nothing.
    """
    logger = logging.getLogger(__name__)

    if args.output is not sys.stdout:
        output_path = args.output / f"failed_cazy_scrapes_{time_stamp}.txt"

        with open(output_path, "a") as fh:
            for url in failed_urls:
                fh.write(f"{url}\n")

    else:
        logger.error("The following items were not scraped:")
        for url in failed_urls:
            logger.error(url)

    return


def write_out_failed_proteins(sql_failures, time_stamp, args):
    """Write out the names of proteins which raised errors when being added to the local db.
    :param sql_failures: list, the names of proteins that were unsuccessfully added to the db
    :param args: cmd args parser
    Return nothing.
    """
    logger = logging.getLogger(__name__)

    if args.output is not sys.stdout:
        output_path = args.output / f"failed_db_protein_additions_{time_stamp}.txt"

        with open(output_path, "a") as fh:
            for fail in sql_failures:
                fh.write(f"{fail}\n")

    else:
        logger.error("The following proteins were not entered into database:")
        for fail in sql_failures:
            logger.error(fail)

    return


def write_out_fasta(record, genbank_accession, args):
    """Write out GenBank protein record to a FASTA file.

    :param record: SeqIO parsed record
    :param genbank_accession: str, accession number of the protein sequence in NCBI.GenBank
    :param args: cmd-line arguments parser

    Return nothing.
    """
    fasta_name = f"{genbank_accession}.fasta"
    fasta_name = args.write / fasta_name

    with open(fasta_name, "w") as fh:
        SeqIO.write(record, fh, "fasta")

    return
