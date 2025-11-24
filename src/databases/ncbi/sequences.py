#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) James Hutton Institute 2024
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
# Author:
# Emma E. M. Hobbs
#
# Contact
# ehobbs@ebi.ac.uk
#
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
"""Get sequences from NCBI GenBank"""

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
from typing import Tuple

from src.cache.ncbi import cache_connection_errors, cache_invalid_ids
from src.databases.ncbi import post_acc_to_entrez, get_protein_accession, validate_batch
from src.sql.interface.add_data.add_ncbi_seqs import update_ncbi_seqs

logger = logging.getLogger(__name__) 


def get_seqs_from_ncbi(
    protein_accs: list[str],
    cache_dir: Path,
    time_stamp: str,
    args: Namespace
) -> dict[str, int]:
    """Get sequences from NCBI and persist them batch by batch.

    Args:
        protein_accs: List of protein accessions to retrieve
        cache_dir: Directory for caching failed requests
        time_stamp: Timestamp for the operation
        args: Command line arguments
 
    Returns:
        Dictionary with statistics about the operation
    """
    batches = get_chunks_list(protein_accs, args.batch_size)
    stats = {
        'sequences_added': 0,
        'sequences_updated': 0,
        'batches_processed': 0,
        'failed_batches': 0
    }
    connection_err_cache = []
    invalid_id_cache = []

    def process_batch(batch):
        nonlocal stats, connection_err_cache, invalid_id_cache

        validated_batch, connection_err_cache = validate_batch(batch, connection_err_cache, args)

        if not validated_batch:
            stats['failed_batches'] += 1
            return

        epost_webenv, epost_query_key, invalid_err, connection_err = post_acc_to_entrez(validated_batch, args)
        if connection_err:
            connection_err_cache.append(validated_batch)
            stats['failed_batches'] += 1
            return
        if invalid_err:
            # If batch size is 1, cache as invalid ID
            if len(validated_batch) == 1:
                invalid_id_cache.append(validated_batch[0])
            else:
                # Split batch in half and retry
                mid = len(validated_batch) // 2
                process_batch(validated_batch[:mid])
                process_batch(validated_batch[mid:])
            stats['failed_batches'] += 1
            return

        batch_seqs, invalid_err, connection_err = fetch_seqs_from_entrez(
            epost_webenv, epost_query_key, protein_accs, args
        )
        if connection_err:
            connection_err_cache.append(validated_batch)
            stats['failed_batches'] += 1
            return
        if invalid_err:
            if len(validated_batch) == 1:
                invalid_id_cache.append(validated_batch[0])
            else:
                mid = len(validated_batch) // 2
                process_batch(validated_batch[:mid])
                process_batch(validated_batch[mid:])
            stats['failed_batches'] += 1
            return

        # Persist batch to database immediately
        if batch_seqs:
            batch_stats = update_ncbi_seqs(batch_seqs, args.database, time_stamp, args.update)
            stats['sequences_added'] += batch_stats['added']
            stats['sequences_updated'] += batch_stats['updated']

        stats['batches_processed'] += 1
        logger.info("Processed batch %d/%d: %d sequences retrieved", 
                   stats['batches_processed'], len(batches), len(batch_seqs))

    for i, batch in enumerate(tqdm(batches, desc="Querying NCBI and persisting sequences")):
        process_batch(batch)

    # Cache any errors
    if connection_err_cache:
        cache_connection_errors(connection_err_cache, cache_dir)

    if invalid_id_cache:
        cache_invalid_ids(invalid_id_cache, cache_dir)

    logger.info("Sequence retrieval completed:")
    logger.info("  - Sequences added: %d", stats['sequences_added'])
    logger.info("  - Sequences updated: %d", stats['sequences_updated'])
    logger.info("  - Batches processed: %d", stats['batches_processed'])
    logger.info("  - Failed batches: %d", stats['failed_batches'])

    return stats


def fetch_seqs_from_entrez(
    epost_webenv,
    epost_query_key,
    all_acc: list[str],
    args: Namespace
) -> Tuple[dict[str, Seq], bool, bool]:
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
        logger.warning("Runtime error occurred during Entrez call. Error: %s", err)
        invalid_err = True

    except (IncompleteRead, CorruptedXMLError) as err:
        logger.warning("IncompleteRead or CorruptedXMLError during Entrez call: %s", err)
        connection_err = True

    except NotXMLError as err:
        logger.warning("NotXMLError during Entrez call: %s", err)
        invalid_err = True

    except (TypeError, AttributeError) as err:
        logger.warning("TypeError or AttributeError during Entrez call: %s", err)
        invalid_err = True

    except Exception as err:
        logger.warning("Unhandled exception during Entrez call: %s", err)
        invalid_err = True

    return new_seqs, invalid_err, connection_err
