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

from argparse import Namespace
from http.client import IncompleteRead
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Entrez.Parser import NotXMLError, CorruptedXMLError
from Bio.Seq import Seq
from saintBioutils.genbank import entrez_retry
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm

from src.cache.ncbi import cache_connection_errors
from src.databases.ncbi import call_entrez, get_protein_accession


logger = logging.getLogger(__name__) 


def get_seqs_from_ncbi(
    protein_accs: list[str],
    cache_dir: Path,
    args: Namespace
) -> dict[str, Seq]:
    batches = get_chunks_list(protein_accs, args.batch_size)
    seq_dict = {}
    connection_err_cache = []
    for batch in tqdm(batches, desc="Querying NCBI"):

        validated_batch, connection_err_cache = validate_batch(batch, connection_err_cache, args)

        if validated_batch:
            epost_webenv, epost_query_key, invalid_err, connection_err = post_acc_to_entrez(validated_batch, args)
            if connection_err:
                connection_err_cache.append(validated_batch)
                continue

            new_seqs, invalid_err, connection_err = fetch_seqs_from_entrez(epost_webenv, epost_query_key, protein_accs, args)
            if connection_err:
                connection_err_cache.append(validated_batch)
                continue
            seq_dict.update(new_seqs)

    if connection_err_cache:
        cache_connection_errors(connection_err_cache, cache_dir)

    return seq_dict


def validate_batch(
    batch: list[str],
    connection_err_cache: list[list[str]],
    args: Namespace
) -> list[str]:
    """Valididate a batch by attempting to post and splitting recursively if it fails.
    Batches may contain IDs that are not in NCBI (termed in valid accessions)."""
    try:
        epost_webenv, epost_query_key, invalid_err, connection_err = post_acc_to_entrez(batch, args)
        if connection_err:
            connection_err_cache.append(batch)
            return [], connection_err_cache
        if invalid_err:
            if len(batch) == 1:
                logger.warning("Invalid ID. Accessions '%s' not listed in NCBI", batch[0])
                return []
            raise ValueError("Batch post failed")

    except ValueError:
        mid = len(batch) // 2
        valid_left = validate_batch(batch[:mid], connection_err_cache, args)
        valid_right = validate_batch(batch[mid:], connection_err_cache, args)
        return valid_left + valid_right, connection_err_cache

    return batch, connection_err_cache


def post_acc_to_entrez(batch: list[str], args: Namespace) -> tuple[any, any, bool, bool]:
    """Post NCBI protein accessions to NCBI via Entrez, and capture error message if one returned.
    Return Entrez ePost web environment and query key.
    """
    invalid_err, connection_err = False, False
    epost_webenv, epost_query_key = None, None

    epost_result, invalid_err, connection_err = call_entrez(
        Entrez.epost,
        retries=args.retries,
        db="Protein",
        id=",".join(batch)
    )

    if not invalid_err or not connection_err:
        try:
            epost_webenv = epost_result["WebEnv"]
            epost_query_key = epost_result["QueryKey"]
        except (TypeError, AttributeError) as err:
            logger.warning("Error occurred when query NCBI. Error:\n%sCaching this batches accessions", err)
            connection_err = True

    return epost_webenv, epost_query_key, invalid_err, connection_err


def fetch_seqs_from_entrez(
    epost_webenv,
    epost_query_key,
    all_acc: list[str],
    args: Namespace
):
    """Retrieve protein sequences from NCBI using ePost results."""
    invalid_err = False
    connection_err = False
    new_seqs = {}

    try:
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            db="Protein",
            query_key=epost_query_key,
            WebEnv=epost_webenv,
            rettype="fasta",
            retmode="text"
        ) as seq_handle:
            for record in SeqIO.parse(seq_handle, "fasta"):
                retrieved_acc = get_protein_accession(record)
                if not retrieved_acc:
                    logger.warning(
                        "Could not extract protein accession from '%s'."
                        "Sequence will not be added to the database", record.id
                    )
                    continue
                if retrieved_acc not in all_acc:
                    logger.warning(
                        "Protein accession from NCBI record ID '%s' (id: '%s) does not match any IDs in the batches.\n"
                        "Sequence will not be added to the database", record.id, retrieved_acc
                    )
                    continue
                new_seqs[retrieved_acc] = Seq(record.seq)

    except RuntimeError as err:
        logger.warning("Runtime error occurred during Entrez call. Error returned:\n%s\nRecord:%s", err, record.id)
        invalid_err = True

    except (IncompleteRead, CorruptedXMLError) as err:
        logger.warning("IncompleteRead or CorruptedXMLError during Entrez call:\n%s\nRecord:%s", err, record.id)
        connection_err = True

    except NotXMLError as err:
        logger.warning("NotXMLError during Entrez call:\n%s\nRecord:%s", err, record.id)
        connection_err = True

    except (TypeError, AttributeError) as err:
        logger.warning("TypeError or AttributeError during Entrez call:\n%s\nRecord:%s", err, record.id)
        connection_err = True

    except Exception as err:
        logger.warning("Unhandled exception during Entrez call:\n%s\nRecord:%s", err, record.id)
        connection_err = True

    return new_seqs, invalid_err, connection_err
