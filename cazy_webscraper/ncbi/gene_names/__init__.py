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
"""Get the protein records for a provided list of gene names"""


import logging

from Bio import Entrez
from saintBioutils.genbank import entrez_retry
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm


def get_linked_ncbi_accessions(uniprot_dict, args):
    """Get the NCBI protein version accessions for the gene names retrieved from UniProt.

    :param uniprot_dict: {uniprot_acc: {gene_name: str, protein_name: str, pdb: set, ec: set, sequence:str, seq_data:str}}
    :param args: CLI parser

    Return uniprot_dict with the ncbi protein accessions added.
    """
    Entrez.email = args.email
    logger = logging.getLogger(__name__)

    gene_names = {}  # to process
    for uniprot_acc in uniprot_dict:
        try:
            uniprot_dict[uniprot_acc]['genbank_accession']  # already have the accession
        except KeyError:
            if uniprot_dict[uniprot_acc]['gene_name'] != 'nan':
                gene_names[uniprot_dict[uniprot_acc]['gene_name']] = uniprot_acc

    if len(list(gene_names.keys())) == 0:  # do not need to query NCBI to get the genbank accessions
        return uniprot_dict

    ncbi_queries = get_chunks_list(
        list(gene_names.keys()),
        args.ncbi_batch_size,
    )

    invalid_gene_names = {'nan'}   # gene names no longer listed in NCBI
    failed_batches = []
    failed_names = {}

    for batch in tqdm(ncbi_queries, desc="Batch quering NCBI to get protein accesions for gene names"):
        batch = list(set(batch).difference(set(invalid_gene_names)))
        uniprot_dict, invalid_gene_names, failed_batches = process_batch(
            batch,
            gene_names,
            uniprot_dict,
            invalid_gene_names,
            failed_batches,
        )

    if len(failed_batches) != 0:
        # remove invalid IDs
        for batch in tqdm(failed_batches, desc="Retrying failed batches"):
            batch = list(set(batch).difference(set(invalid_gene_names)))

            uniprot_dict, invalid_gene_names, processed_batch = process_batch(
                batch,
                gene_names,
                uniprot_dict,
                invalid_gene_names,
                [],  # pass an empty list
            )

            if len(processed_batch) != 0:  # batch still causing issues, retry names individually
                for name in batch:
                    failed_names[name] = 0

    if len(list(failed_names.keys())) != 0:
        while len(list(failed_names.keys())) != 0:
            # remove invalid ids and don't retry
            for name in invalid_gene_names:
                try:
                    failed_names[name]
                    logger.warning(
                        f"Gene name {name} retrieved from UniProt record {gene_names[name]}\n"
                        "is no longer in NCBI, therefore, not adding protein data to the\n"
                        f"local CAZyme database for UniProt record {gene_names[name]}"
                    )
                    del failed_names[name]
                    del uniprot_dict[gene_names[name]]
                except KeyError:
                    continue

            names_to_process = list(failed_names.keys())
            for name in tqdm(names_to_process, desc="Retrying failed gene names"):
                uniprot_dict, invalid_gene_names, processed_batch = process_batch(
                    [name],  # gene name must be in a list
                    gene_names,
                    uniprot_dict,
                    invalid_gene_names,
                    [],
                    uniprot_acc=gene_names[name],
                )

                if len(processed_batch) != 0:
                    failed_names[name] += 1

                if failed_names[name] >= args.retries:
                    del failed_names[name]
                    logger.error(
                        f"Could not retrieved NCBI protein version accession for gene name {name}\n"
                        "from NCBI\n"
                        f"Not adding protein data from UniProt accession {gene_names[name]} because\n"
                        f"could not link its gene name {name} to a NCBI proten record"
                    )

    return uniprot_dict


def process_batch(batch, gene_names, uniprot_dict, invalid_gene_names, failed_batches, uniprot_acc=None):
    """Coordinate processing batch query results.

    :param batch: list of gene_names
    :param uniprot_dict: {uniprot_acc: {gene_name: str, protein_name: str, pdb: set, ec: set, sequence:str, seq_data:str}}
    :param invalid_gene_names: set of gene names not in NCBI
    :param failed_batches: list of failed batches
    :param UniProt_acc: str uniprot record entry accession

    Return uniprot_dict, invalid_gene_names and failed_batches
    """
    record, success = query_ncbi(batch, gene_names=",".join(batch))

    if success == 'invalid ID':
        invalid_gene_names.add(batch[0])
        return uniprot_dict, invalid_gene_names, failed_batches
    
    elif success == 'retry':
        failed_batches.append(batch)
        return uniprot_dict, invalid_gene_names, failed_batches

    epost_webenv = record["WebEnv"]
    epost_query_key = record["QueryKey"]

    record, success = query_ncbi(batch, query_key=epost_query_key, webEnv=epost_webenv)

    if success == 'invalid ID':
        invalid_gene_names.add(batch[0])
        return uniprot_dict, invalid_gene_names, failed_batches
    
    elif success == 'retry':
        failed_batches.append(batch)
        return uniprot_dict, invalid_gene_names, failed_batches

    for prot_record in tqdm(record, desc="Parsing query output"):
        gene_name = None
        uniprot_acc = None
        genbank_accession = None

        for i in prot_record['GBSeq_feature-table']:
            if i['GBFeature_key'] == 'CDS':
                for j in i['GBFeature_intervals']:
                    gene_name = j['GBInterval_accession'].split(".")[0]
                    try:
                        uniprot_acc = gene_names[gene_name]

                    except KeyError:
                        return uniprot_dict, invalid_gene_names, failed_batches
                    
                for k in i['GBFeature_intervals']:
                    genbank_accession = k['GBInterval_accession'].split(".")[0]
                    
        if gene_name is not None:
            uniprot_dict[uniprot_acc]['genbank_accession'] = genbank_accession

    return uniprot_dict, invalid_gene_names, failed_batches


def query_ncbi(batch, gene_names=None, query_key=None, webEnv=None, uniprot_acc=None):
    """Post IDs or retrieved results from positing IDs to NCBI.

    :param batch: list of gene names
    :param gene_names: str of items separted with single comma
    :param query_key: str, from posting IDs
    :param WebEnv: str, from posting IDs
    :param UniProt_acc: str uniprot record entry accession

    Return the Entrez record, str indicating if the gene id is valid or the query should 
    be retried
    """
    logger = logging.getLogger(__name__)
    success = None
    record = None

    try:
        if gene_names is not None:  # post IDs
            process = "ePost"
            record = Entrez.read(
                entrez_retry(
                    10,
                    Entrez.epost,
                    db="Protein",
                    id=gene_names,
                ),
                validate=False,
            )

        else:
            process = "eFetch"
            with entrez_retry(
                10,
                Entrez.efetch,
                db="Protein",
                query_key=query_key,
                WebEnv=webEnv,
                retmode='xml',
            ) as record_handle:
                record = Entrez.read(record_handle, validate=False)

    except RuntimeError as err:
        if repr(err).startswith("RuntimeError('Some IDs have invalid value and were omitted.") or repr(err).startswith("RuntimeError('Empty ID list; Nothing to store')"):
            if len(batch) == 1:
                logger.warning(
                    f"Gene name {batch[0]} retrieved for UniProt entry '{uniprot_acc}' "
                    "is no longer listed in CAZy\n"
                    f"Not adding protein data for UniProt accessions {uniprot_acc}"
                )
                success = "invalid ID"
            
            else:
                logger.warning(
                    "Batch contains a gene name not in NCBI\n"
                    "Will identify invalid gene names later"
                )
                success = "retry"
        else:
            logger.warning(
                f"Runtime error occurenced when batch quering NCBI ({process})\n"
                "Will retry batch later\n"
                "Error retrieved:\n"
                f"{repr(err)}\n"
                f"Batch:\n{batch}"
            )
            success = "retry"

    except (TypeError, AttributeError) as err:  # if no record is returned from call to Entrez
        logger.warning(
            f"Error occurenced when batch quering NCBI ({process})\n"
            "Will retry batch later\n"
            "Error retrieved:\n"
            f"{repr(err)}\n"
            f"Batch:\n{batch}"
        )
        success = "retry"

    return record, success
