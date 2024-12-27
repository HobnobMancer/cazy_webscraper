#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
# Author:
# Emma E. M. Hobbs

# Contact
# ehobbs@ebi.ac.uk

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
"""Read and write cache data retrieved from the remote NCBI database"""


import argparse
import logging
import json

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from src import closing_message
from src.ncbi import get_protein_accession
from src.ncbi.sequences import get_protein_accession


logger = logging.getLogger(__name__)


def get_cache_seqs(
    start_time: str,
    args: argparse.ArgumentParser
) -> dict[str, Seq]:
    """Extract protein sequences from FASTA and/or JSON file, which will be added to the
    local CAZyme database
    """
    seq_dict = {}

    if args.seq_dict:
        seq_dict = load_json_seqs(args.seq_dict, seq_dict)
        if not seq_dict:
            closing_message("Get GenBank seqs", start_time, args, early_term=True)

    if args.seq_file:
        seq_dict = load_fasta_seq(args.seq_file, seq_dict)
        if not seq_dict:
            closing_message("Get GenBank seqs", start_time, args, early_term=True)

    logger.warning("Retrieved %s from cache", len(seq_dict.keys()))

    return seq_dict


def load_json_seqs(seq_dict_path: Path, seq_dict: dict[str, Seq]) -> dict[str, Seq]:
    logger.warning("Getting sequences from JSON cache:\n%s", seq_dict_path)
    try:
        with open(seq_dict_path, "r") as fh:
            cache_dict = json.load(fh)
    except FileNotFoundError:
        logger.error(
            "Could not find JSON file of protein sequences at:%s\n"
            "Check the path is correct. Terminating program",
            seq_dict_path
        )
        return

    for key in cache_dict:  # convert strs to SeqRecords
        seq_dict[key] = Seq(cache_dict[key])
    return seq_dict


def load_fasta_seq(fasta_path: Path, seq_dict: dict[str, Seq]) -> dict[str, Seq]:
    logger.warning("Getting sequences from FASTA cache:\n%s", fasta_path)
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            retrieved_accession = get_protein_accession(record)

            if retrieved_accession is None:
                logger.error(
                    "Could not retrieve a NCBI protein version accession from cache\n"
                    "from the record id '%s'\n"
                    "The sequence from this record will not be added to the db",
                    record.id
                )
                continue

            try:
                if seq_dict[retrieved_accession] != record.seq:
                    logger.warning(
                        "Retrieved seq for %s from JSON file which does NOT match "
                        "the seq in the FASTA file.\n"
                        "Adding seq from the FASTA file to the local CAZyme database\n"
                        "JSON seq: %s\n"
                        "FASTA seq: %s",
                        retrieved_accession,
                        seq_dict[retrieved_accession],
                        record.seq
                    )
                    seq_dict[retrieved_accession] = record.seq
            except KeyError:
                seq_dict[retrieved_accession] = record.seq

    except FileNotFoundError:
        logger.error(
            "Could not find FASTA file of protein sequences at:%s\n"
            "Check the path is correct. Terminating program",
            fasta_path
        )
        return

    return seq_dict


def cache_connection_errors(batches: list[list[str]], cache_dir: Path):
    cache_path = cache_dir / "failed_connections.out"
    with open(cache_path, "w") as fh:
        for batch in batches:
            text = '\n'.join(batch)
            fh.write(f"{text}\n")
