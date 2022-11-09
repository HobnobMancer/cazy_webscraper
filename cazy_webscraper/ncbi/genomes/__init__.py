#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
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
"""Module for interacting with the NCBI databases to retrieve genomic data."""


import logging

from Bio import Entrez
from saintBioutils.genbank import entrez_retry
from tqdm import tqdm

from cazy_webscraper.ncbi import post_ids


def get_nuccore_ids(batch, failed_batches, args, retry=False):
    """Retrieve the IDs of nuccore records linkded to a list of protein version accessions

    :param batch: list of GenBank protein accessions
    :param failed_batches: dict listing batches for which data was not retrieved
    :param args: cmd-line args parser
    :param retry: bool, is this a retry of previously failed retrieval of genome data?

    Return set() of nuccore db IDs and failed_batches
    """
    logger = logging.getLogger(__name__)

    nuccore_ids = set()

    logger.info("Posting protein accessions to NCBI")

    try:
        query_key, web_env = post_ids(batch, "Protein", args)

        if query_key is None:
            failed_batches['proteins'].append(batch)
            return nuccore_ids, failed_batches

    except RuntimeError:
        if retry:
            logger.warning(f"{str(batch[0])} is not listed in NCBI", exc_info=1)
        else:
            logger.warning("Batch contains invalid NCBI Protein db accessions", exc_info=1)
            failed_batches['proteins'].append(batch)

        return nuccore_ids, failed_batches

    try:
        with entrez_retry(
            args.retries,
            Entrez.elink,
            query_key=query_key,
            WebEnv=web_env,
            dbfrom="Protein",
            db="Nuccore",
            linkname="protein_nuccore",
        ) as handle:
            nuccore_records = Entrez.read(handle, validate=False)

    except (TypeError, AttributeError, RuntimeError):
        logger.warning(
            f"Entrez failed to link Protein records to nuccore records numbers\n",
            exc_info=1,
        )
        return nuccore_ids, failed_batches

    for record in tqdm(nuccore_records, desc="Get Nuccore record IDs"):
        if len(record['LinkSetDb']) != 0:
            for nuc_id in record['LinkSetDb'][0]['Link']:
                nuccore_ids.add(nuc_id['Id'])
    
    return nuccore_ids, failed_batches


def get_assembly_ids(nuccore_ids, failed_batches, args, retry=False):
    """Retrieve the IDs of assembly records linkded to a list of nuccore record Ids

    :param nuccore_ids: list of uccore ncbi db IDs
    :param failed_batches: dict listing batches for which data was not retrieved
    :param args: cmd-line args parser
    :param retry: bool, is this a retry of previously failed retrieval of genome data?

    Return set() of assebly db IDs and failed_batches
    """
    logger = logging.getLogger(__name__)

    assembly_ids = set()

    logger.info("Posting nuccore IDs")
    # post nuccore IDs
    try:
        query_key, web_env = post_ids(nuccore_ids, "Nuccore", args)

        if query_key is None:
            failed_batches['nuccores'].append(nuccore_ids)
            return assembly_ids, failed_batches

    except RuntimeError:
        if retry:
            logger.warning(f"Nuccore ID '{nuccore_ids[0]}' is not listed in NCBI")
        else:
            logger.warning("Batch contains invalid NCBI Nuccore db IDs")
            failed_batches['nuccores'].append(nuccore_ids)

        return assembly_ids, failed_batches

    logger.info("Getting linked assembly IDs")
    try:
        logger.info("Try")
        with entrez_retry(
            args.retries,
            Entrez.elink,
            query_key=query_key,
            WebEnv=web_env,
            dbfrom="Nuccore",
            db="Assembly",
            linkname="nuccore_assembly",
        ) as handle:
            linked_records = Entrez.read(handle, validate=False)
    except (TypeError, AttributeError, RuntimeError) as err:
        logger.warning(f"Failed to link nuccore records to assembly records:\n{err}")
        return assembly_ids, failed_batches

    for record in linked_records:
        for index in range(len(record['LinkSetDb'][0]['Link'])):
            assembly_id = record['LinkSetDb'][0]['Link'][index]['Id']
            assembly_ids.add(assembly_id)

    return assembly_ids, failed_batches


def get_assembly_data(assembly_ids, failed_batches, parsed_assembly_ids, args, retry=False):
    """Retrieve the data for assemblies represented by their NCBI Assembly DB ID

    :param assembly_ids: list of assembly ncbi db IDs
    :param failed_batches: dict listing batches for which data was not retrieved
    :param args: cmd-line args parser
    :param parsed_assembly_ids: set of ncbi assembly db IDs that have already been parsed
    :param retry: bool, is this a retry of previously failed retrieval of genome data?

    Return dict of assembly meta data and failed batches dict
    """
    logger = logging.getLogger(__name__)

    genome_dict = {}

    # post assembly IDs
    try:
        query_key, web_env = post_ids(assembly_ids, "Assembly", args)

        if query_key is None:
            failed_batches.append(assembly_ids)
            return genome_dict, failed_batches
    except RuntimeError:
        if retry:
            logger.warning(f"Data for Assembly ID '{assembly_ids[0]}' could not be retrieved from NCBI")
        else:
            logger.warning("Batch contains invalid NCBI Assembly IDs")
            failed_batches['assemblies'].append(assembly_ids)
        return genome_dict, failed_batches

    try:
        with entrez_retry(
            args.retries,
            Entrez.esummary,
            db="Assembly",
            query_key=query_key,
            WebEnv=web_env,
            rettype="docsum",
            retmode="xml",
        ) as record_handle:
            assembly_records = Entrez.read(record_handle, validate=False)
    except (TypeError, AttributeError, RuntimeError) as err:
        logger.warning(f"Failed to retrieve Assembly records:\n{err}")
        return genome_dict, failed_batches

    for genome in assembly_records['DocumentSummarySet']['DocumentSummary']:
        assembly_name = genome['AssemblyName']

        try:
            genome_dict[assembly_name]
            continue
        except KeyError:
            pass

        gbk_uid = genome['GbUid']
        ref_uid = genome['RsUid']

        try:
            gbk_acc = genome['Synonym']['Genbank']
        except KeyError:
            gbk_acc = None
        
        try:
            refseq_acc = genome['Synonym']['RefSeq']
        except KeyError:
            refseq_acc = None

        try:
            gbk_ftp_path = genome['FtpPath_GenBank']
            file_stem = gbk_ftp_path.split("/")[-1]
            gbk_url = f"{gbk_ftp_path}/{file_stem}_feature_table.txt.gz"
        except KeyError:
            gbk_url = None
        
        try:
            ref_ftp_path = genome['FtpPath_RefSeq']
            file_stem = ref_ftp_path.split("/")[-1]
            ref_url = f"{ref_ftp_path}/{file_stem}_feature_table.txt.gz"
        except KeyError:
            ref_url = None

        genome_dict[assembly_name] = {
            "gbk_acc": gbk_acc,
            "refseq_acc": refseq_acc,
            "gbk_uid": gbk_uid,
            "refseq_uid": ref_uid,
            "gbk_url": gbk_url,
            "refseq_url": ref_url,
        }

        parsed_assembly_ids.add(gbk_uid)
        parsed_assembly_ids.add(ref_uid)

    return genome_dict, failed_batches, parsed_assembly_ids
