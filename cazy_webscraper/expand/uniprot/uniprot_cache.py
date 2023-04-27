#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
#
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
"""Retrieve data from cache files"""


import json
import logging


def get_uniprot_cache(gbk_dict, args):
    """Get cached UniProt data, or return empty dict and set.
    
    :param gbk_dict: {gbk acc: local db id}
    :param args: CLI args parser

    Return
    Dict:
    uniprot_data[ncbi_acc] = {
        'uniprot_acc': uniprot_acc,
        'uniprot_entry_id': uniprot_entry_id,
        'protein_name': protein_name,
        'ec_numbers': ec_numbers,
        'sequence': sequence,
        'pdbs': all_pdbs,
    }
    and list of GenBank accessions to download UniProt data for
    """
    logger = logging.getLogger(__name__)

    # if using cachce skip accession retrieval
    # uniprot_data[ncbi_acc] = {
    #     'uniprot_acc': uniprot_acc,
    #     'uniprot_entry_id': uniprot_entry_id,
    #     'protein_name': protein_name,
    #     'ec_numbers': ec_numbers,
    #     'sequence': sequence,
    #     'pdbs': all_pdbs,
    # }
    uniprot_dict = {}
    gbk_data_to_download = []

    if args.use_uniprot_cache is not None:
        logger.warning(f"Getting UniProt data from cache: {args.use_uniprot_cache}")

        with open(args.use_uniprot_cache, "r") as fh:
            uniprot_dict = json.load(fh)
    
    if args.skip_download:  # only use cached data
        return uniprot_dict, gbk_data_to_download

    # else: check for which GenBank accessions data still needs be retrieved from UniProt
    # if some of the data is used from a cache, if no data is provided from a cache
    # retrieve data for all GenBank accesisons matching the provided criteria
    if len(list(uniprot_dict.keys())) != 0: 
        for uniprot_acc in tqdm(uniprot_dict):
            gbk_data_to_download.append(uniprot_dict[uniprot_acc]['genbank_accession'])
    else:  # get data for all GenBank accessions from the local db matching the user criteria
        gbk_data_to_download = list(gbk_dict.keys())
    
    return uniprot_dict, gbk_data_to_download


def cache_uniprot_data(uniprot_dict, cache_dir, time_stamp):
    """Cache data retrieved from UniProt.

    :param

    Return nothing
    """
    # uniprot_dict[ncbi_acc] = {
    #     'uniprot_acc': uniprot_acc,
    #     'uniprot_entry_id': uniprot_entry_id,
    #     'protein_name': protein_name,
    #     'ec_numbers': ec_numbers,
    #     'sequence': sequence,
    #     'sequence_date': sequence_date,
    #     'pdbs': all_pdbs,
    # }

    # cache updated UniProt data
    for uniprot_accession in uniprot_dict:
        try:
            uniprot_dict[uniprot_accession]['ec_numbers'] = list(uniprot_dict[uniprot_accession]['ec_numbers'])
        except KeyError:
            pass
        try:
            uniprot_dict[uniprot_accession]['pdbs'] = list(uniprot_dict[uniprot_accession]['pdbs'])
        except KeyError:
            pass

    uniprot_acc_cache = cache_dir / f"uniprot_data_{time_stamp}.json"
    with open(uniprot_acc_cache, "w") as fh:
        json.dump(uniprot_dict, fh) 
