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

from saintBioutils.genbank import entrez_retry
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm


def get_linked_ncbi_accessions(uniprot_dict):
    """Get the NCBI protein version accessions for the gene names retrieved from UniProt.

    :param uniprot_dict: {uniprot_acc: {gene_name: str, protein_name: str, pdb: set, ec: set, sequence:str, seq_data:str}}

    Return uniprot_dict with the ncbi protein accessions added.
    """
    logger = logging.getLogger(__name__)

    gene_names = {}
    for uniprot_acc in uniprot_dict:
        gene_names[uniprot_dict[uniprot_acc]['gene_name']] = uniprot_acc

    ncbi_queries = get_chunks_list(
        list(gene_names.keys()),
        args.ncbi_batch_size,
    )

    invalid_gene_names = set()   # gene names no longer listed in NCBI
    failed_batches = []

    for batch in tqdm(ncbi_queries, desc="Batch quering NCBI to get protein accesions for gene names"):
        uniprot_dict, invalid_gene_names, failed_batches = process_batch(
            batch,
            uniprot_dict,
            invalid_gene_names,
            failed_batches,
        )

    if len(failed_batches) != 0:
        failed_names = {}

        # remove invalid IDs
        for batch in failed_batches:
            batch = list(set(batch).difference(set(invalid_ids)))

            uniprot_dict, invalid_gene_names, processed_batch = process_batch(
                batch,
                uniprot_dict,
                invalid_gene_names,
                [],  # pass an empty list
            )

            if len(processed_batch) != 0:  # batch still causing issues, retry names individually
                for name in batch:
                    failed_names[name] = 0

    if len(list(failed_names.keys())) != 0:
        while len(list(failed_names.keys())) != 0:
            names_to_process = list(failed_names.keys())
            for name in names_to_process:
                uniprot_dict, invalid_gene_names, processed_batch = process_batch(
                    [name],  # gene name must be in a list
                    uniprot_dict,
                    invalid_gene_names,
                    [],
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


def process_batch(batch, uniprot_dict, invalid_gene_names, failed_batches):
    """Coordinate processing batch query results.

    :param batch: list of gene_names
    :param uniprot_dict: {uniprot_acc: {gene_name: str, protein_name: str, pdb: set, ec: set, sequence:str, seq_data:str}}
    :param invalid_gene_names: set of gene names not in NCBI
    :param failed_batches: list of failed batches

    Return uniprot_dict, invalid_gene_names and failed_batches
    """
    record, success = query_ncbi(id=",".join(batch))

    if success == 'invalid ID':
        invalid_gene_names.add(batch[0])
        return uniprot_dict, invalid_gene_names, failed_batches
    
    elif success == 'retry':
        failed_batches.append(batch)
        return uniprot_dict, invalid_gene_names, failed_batches

    epost_webenv = epost_result["WebEnv"]
    epost_query_key = epost_result["QueryKey"]

    record, success = query_ncbi(query_key=epost_query_key, WebEnv=epost_webenv)

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
                    genebank_accession = k['GBInterval_accession'].split(".")[0]
                    
        if gene_name is not None:
            uniprot_dict[uniprot_acc]['genbank_accession'] = genebank_accession

    return uniprot_dict, invalid_gene_names, failed_batches


def query_ncbi(id=None, query_key=None, WebEnv=None):
    """Post IDs or retrieved results from positing IDs to NCBI.

    :param id: str of items separted with single comma
    :param query_key: str, from posting IDs
    :param WebEnv: str, from posting IDs

    Return the Entrez record, str indicating if the gene id is valid or the query should 
    be retried
    """
    logger = logging.getLogger(__name__)
    success = None
    record = None

    try:
        if id is not None:  # post IDs
            epost_result = Entrez.read(
                entrez_retry(
                    10,
                    Entrez.epost,
                    db="Protein",
                    id=",".join(batch),
                ),
                validate=False,
            )

        else:
            with entrez_retry(
                10,
                Entrez.efetch,
                db="Protein",
                query_key=epost_query_key,
                WebEnv=epost_webenv,
                retmode='xml',
            ) as record_handle:
                record = Entrez.read(record_handle, validate=False)

    except RuntimeError as err:
        if repr(err).startswith("RuntimeError('Some IDs have invalid value and were omitted."):
            if len(batch) == 1:
                logger.warning(
                    f"Gene name {batch[0]} retrieved for UniProt entry {gene_names[batch[0]]} "
                    "is no longer listed in CAZy\n"
                    f"Not adding protein data for UniProt accessions {gene_names[batch[0]]}"
                )
                success = "invalid ID"
            
            else:
                failed_batches.append(batch)
                logger.warning(
                    "Batch contains a gene name not in NCBI\n"
                    "Will identify invalid gene names laters"
                )
                success = "retry"
        else:
            logger.warning(
                "Error occurenced when batch quering NCBI\n"
                "Will retry batch later\n"
                "Error retrieved:\n"
                f"{repr(err)}"
            )
            success = "retry"

    except (TypeError, AttributeError) as err:  # if no record is returned from call to Entrez
        logger.warning(
            "Error occurenced when batch quering NCBI\n"
            "Will retry batch later\n"
            "Error retrieved:\n"
            f"{repr(err)}"
        )
        success = "retry"

    return record, success
