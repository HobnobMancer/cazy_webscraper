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
