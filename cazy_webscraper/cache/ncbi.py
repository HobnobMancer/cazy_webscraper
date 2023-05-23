#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Cache data retrieved from the remove NCBI database"""


import logging
import json

from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord


def get_cache_seqs(start_time, args):
    """Extract protein sequences from FASTA and/or JSON file, which will be added to the
    local CAZyme database

    :param seq_dict: dict, {genbank_acc: Bio.Seq}

    Return update seq_dict and list of SeqRecords
    """
    logger = logging.getLogger(__name__)

    seq_dict = {}
    seq_records = []

    if args.seq_dict:
        logger.warning(f"Getting sequences from JSON cache:\n{args.seq_dict}")

        try:
            with open(args.seq_dict, "r") as fh:
                cache_dict = json.load(fh)

        except FileNotFoundError:
            logger.error(
                f"Could not find JSON file of protein sequences at:\n"
                f"{args.seq_dict}\n"
                "Check the path is correct"
                "Terminating program"
            )
            closing_message("Get GenBank seqs", start_time, args, early_term=True)

        # convert strs to SeqRecords
        for key in cache_dict:
            seq_dict[key] = Seq(cache_dict[key])

    if args.seq_file:
        logger.warning(f"Getting sequences from FASTA cache:\n{args.seq_file}")

        try:
            for record in SeqIO.parse(args.seq_file, "fasta"):
                retrieved_accession = get_protein_accession(record)

                if retrieved_accession is None:
                    logger.error(
                        "Could not retrieve a NCBI protein version accession from cache\n"
                        f"from the record id '{record.id}'\n"
                        "The sequence from this record will not be added to the db"
                    )
                    continue

                try:
                    seq_dict[retrieved_accession]
                    if seq_dict[retrieved_accession] != record.seq:
                        logger.warning(
                            f"Retrieved seq for {retrieved_accession} from JSON file which does NOT match "
                            "the seq in the FASTA file.\n"
                            "Adding seq from the FASTA file to the local CAZyme database\n"
                            f"JSON seq: {seq_dict[retrieved_accession]}\n"
                            f"FASTA seq: {record.seq}"
                        )
                        seq_dict[retrieved_accession] = record.seq
                except KeyError:
                    seq_dict[retrieved_accession] = record.seq

        except FileNotFoundError:
            logger.error(
                f"Could not find FASTA file of protein sequences at:\n"
                f"{args.seq_file}\n"
                "Check the path is correct"
                "Terminating program"
            )
            closing_message("Get GenBank seqs", start_time, args, early_term=True)

    for key in seq_dict:
        seq_records.append(SeqRecord(id=key, seq=Seq(seq_dict[key])))

    logger.warning(f"Retrieved {len(seq_records)} from cache")

    return seq_dict, seq_records
