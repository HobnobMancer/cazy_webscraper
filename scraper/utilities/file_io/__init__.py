#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
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
"""Submodule for handling inputs and outputs"""


import logging
import shutil
import sys

from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline



def make_target_directory(output_path, force):
    """Create target output directory for SQL database file and/or JSON file.

    :param output_path: Path to target output file
    :param force: bool, whether to overwrite an existing output file

    Raises FileExistsError if target database file output path exists and 
    force (force overwrite) is False.

    Return Nothing
    """
    logger = logging.getLogger(__name__)

    if output_path.exists():
        if force is True:

            logger.warning(
                "Target output path for database exists.\n"
                "Forced overwritting enabled. Overwritting existing CAZome database."
            )
            return
        
        else:
            logger.warning(
                "Target output path for database exists.\n"
                "Forced overwritting NOT enabled\n."
                "To overwrite the existing CAZome database\n,"
                "rerun cazy_webscraper with the -f flag.\n"
                "Terminating program."
            )
            sys.exit(1)

    # make output directory
    output_path.parent.mkdir(parents=True, exist_ok=force)
    logger.warning(f"Built output directory: {output_path.parent}")

    return


def make_output_directory(output_path, force, nodelete):
    """Create output directory to which multiple output files will be written

    :param output_path: Path to target output dir
    :param force: bool, whether to overwrite an existing output file
    :param nodelete: bool, whether to delete contents in an existing target output dir

    Raises FileExistsError if target database file output path exists and 
    force (force overwrite) is False.

    Return Nothing
    """
    logger = logging.getLogger(__name__)

    if output_path.exists():

        if nodelete is True:
            logger.warning(
                "Target output path for database exists.\n"
                "Forced writing in existing output directory.\n"
                "NOT deleting existing content in the output directory."
            )
            return

        else:
            logger.warning(
                "Target output path for database exists.\n"
                "Forced writing in existing output directory.\n"
                "DELETING existing content in the output directory."
            )

    # make output directory
    try:
        output_path.mkdir(parents=True, exist_ok=force)
    except FileExistsError:
        pass
    logger.warning(f"Built directory to store family files: {output_path}")

    return


def write_out_fasta(record, genbank_accession, args):
    """Write out GenBank protein record to a FASTA file.

    :param record: SeqIO parsed record
    :param genbank_accession: str, accession number of the protein sequence in NCBI.GenBank
    :param args: cmd-line arguments parser

    Return nothing.
    """
    if args.fasta == 'separate':
        fasta_name = f"{genbank_accession}.fasta"
        fasta_name = args.write / fasta_name

        with open(fasta_name, "w") as fh:
            SeqIO.write(record, fh, "fasta")

    else:  # add sequences to FASTA file
        with open(args.fasta, "a") as fh:
            SeqIO.write(record, fh, "fasta")

    return


def write_fasta_for_db(record, genbank_accession, args):
    """Write out protein sequences to FASTA file for building a BLAST db of all retrieved sequences.

    :param record: SeqIO parsed record
    :param genbank_accession: str, accession number of the protein sequence in NCBI.GenBank
    :param args: cmd-line arguments parser

    Return nothing.
    """
    fasta_name = args.blastdb
    fasta_name = fasta_name / "blast_db.fasta"

    with open(args.fasta, "a") as fh:
        SeqIO.write(record, fh, "fasta")

    return


def build_blast_db(args):
    """Build BLAST database of sequences retrieved from GenBank.

    :param args: cmd-line arguments parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    fasta_name = args.blastdb
    fasta_name = fasta_name / "blast_db.fasta"

    # build the command
    cmd_makedb = NcbimakeblastdbCommandline(cmd='makeblastdb', dbtype='prot', input_file=fasta_name)
    # invoke the command
    stdout, stderr = cmd_makedb()

    # check the command was successfully exectured
    if len(stderr) != 0:
        logger.warning()
        print(f"Could not build non-CAZyme db.\nstdout={stdout}\nstderr={stderr}")

    return
