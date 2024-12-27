#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
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
#
# Bio.PDB reference:
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
# implemented in Python. Bioinformatics 19: 2308–2310
"""Retrieve records from db for accessions provided in a list, i.e. a file"""


import logging
import sys

from tqdm import tqdm


def get_user_genbank_sequences(gbk_table_dict, args):
    """Extract protein sequences for GenBank accessions listed in a file

    :param gbk_table_dict: dict {genbank_accession: 'taxa_id': int, 'gbk_id': int}
    :param args: cmd-line args parser

    Return dict {gbk_acc: db_id}
    """
    logger = logging.getLogger(__name__)

    try:
        with open(args.genbank_accessions, "r") as fh:
            lines = fh.read().splitlines()
    except FileNotFoundError:
        logger.warning(
            f"Could not find file of GenBank accessions at {args.genbank_accessions}\n"
            "Check the path is correct\n"
            "Terminating program"
        )
        sys.exit(1)

    gbk_accessions = [line.strip() for line in lines]

    gbk_dict = {}  # {accession: id}

    for gbk_accession in tqdm(gbk_accessions, desc="Getting database IDs for provided GenBank IDs"):
        try:
            gbk_dict[gbk_accession] = gbk_table_dict[gbk_accession]
        except KeyError:
            logging.warning(
                f"Genbank accession {gbk_accession} retrieved from list in file\n"
                "But accession not in the local CAZyme database\n"
                f"Not extracted protein sequences for {gbk_accession}"
            )

    return gbk_dict


def get_user_uniprot_sequences(gbk_table_dict, uniprot_table_dict, args):
    """Extract protein sequences for UniProt accessions listed in a file, and get the corresponing
    GenBank accession and local db GenBank id

    :param gbk_table_dict: dict {genbank_accession: 'taxa_id': int, 'gbk_id': int}
    :param uniprot_table_dict: dict {}
    :param args: cmd-line args parser

    Return dict {gbk_acc: db_id}
    """
    logger = logging.getLogger(__name__)

    try:
        with open(args.uniprot_accessions, "r") as fh:
            lines = fh.read().splitlines()
    except FileNotFoundError:
        logger.warning(
            f"Could not find file of UniProt accessions at {args.uniprot_accessions}\n"
            "Check the path is correct\n"
            "Terminating program"
        )
        sys.exit(1)

    uniprot_accessions = [line.strip() for line in lines]

    gbk_dict = {}

    gbk_table_ids = list(gbk_table_dict.values())
    gbk_table_accs = list(gbk_table_dict.keys())

    for uniprot_accession in tqdm(uniprot_accessions, desc="Getting database Ids for provided UniProt IDs"):
        try:
            gbk_id = uniprot_table_dict[uniprot_accession]['genbank_id']
        except KeyError:
            logging.warning(
                f"UniProt accession {uniprot_accession} retrieved from list in file\n"
                "But accession not in the local CAZyme database\n"
                f"Not extracted protein sequences for {uniprot_accession}"
            )
            continue

        position = gbk_table_ids.index(gbk_id)
        gbk_accession = gbk_table_accs[position]

        gbk_dict[gbk_accession] = gbk_id

    return gbk_dict
