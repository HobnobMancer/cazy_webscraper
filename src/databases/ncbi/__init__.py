#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) James Hutton Institute 2024
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
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
"""Module for interacting with the NCBI databases."""


import logging
import re

from argparse import Namespace
from http.client import IncompleteRead
from typing import Optional

from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Entrez.Parser import NotXMLError, CorruptedXMLError
from saintBioutils.genbank import entrez_retry


logger = logging.getLogger(__name__)


def handle_entrez_errors(err) -> None:
    if isinstance(err, RuntimeError):
        logger.warning("Runtime error occurred during Entrez call. Error returned:\n%s", err)
    elif isinstance(err, (IncompleteRead, CorruptedXMLError)):
        logger.warning("IncompleteRead or CorruptedXMLError during Entrez call:\n%s", err)
    elif isinstance(err, NotXMLError):
        logger.warning("NotXMLError during Entrez call:\n%s", err)
    elif isinstance(err, (TypeError, AttributeError)):
        logger.warning("TypeError or AttributeError during Entrez call:\n%s", err)
    else:
        logger.warning("Unhandled exception during Entrez call:\n%s", err)


def call_entrez(func, retries: int, **kwargs):
    """
    Wrapper for Entrez calls with exception handling.

    Parameters:
    - func: The Entrez function to be called (e.g., Entrez.epost, Entrez.efetch).
    - retries: Number of retries for the Entrez function.
    - kwargs: Parameters to pass to the Entrez function.

    Returns:
    - result: The result of the Entrez function call, or None if errors occur.
    - invalid_err: Flag indicating invalid IDs.
    - connection_err: Flag indicating connection-related errors.
    """
    invalid_err = False
    connection_err = False
    result = None

    try:
        result = Entrez.read(entrez_retry(retries, func, **kwargs), validate=False)

    except RuntimeError as err:
        if repr(err).startswith("RuntimeError('Some IDs have invalid value and were omitted.") \
                or repr(err).startswith("RuntimeError('Empty ID list; Nothing to store')"):
            invalid_err = True
        else:
            logger.warning(
                "Runtime error occurred during Entrez call. Error returned:\n%s", err
            )
            invalid_err = True

    except (IncompleteRead, CorruptedXMLError) as err:
        logger.warning(
            "IncompleteRead or CorruptedXMLError during Entrez call:\n%s", err
        )
        connection_err = True

    except NotXMLError as err:
        logger.warning(
            "NotXMLError during Entrez call:\n%s", err
        )
        connection_err = True

    except (TypeError, AttributeError) as err:
        logger.warning(
            "TypeError or AttributeError during Entrez call:\n%s", err
        )
        connection_err = True

    except Exception as err:
        logger.warning(
            "Unhandled exception during Entrez call:\n%s", err
        )
        connection_err = True

    return result, invalid_err, connection_err


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
                return [], connection_err_cache
            raise ValueError("Batch post failed")

    except ValueError:
        mid = len(batch) // 2
        valid_left, connection_err_cache = validate_batch(batch[:mid], connection_err_cache, args)
        valid_right, connection_err_cache = validate_batch(batch[mid:], connection_err_cache, args)
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

    if not invalid_err and not connection_err:
        try:
            epost_webenv = epost_result["WebEnv"]
            epost_query_key = epost_result["QueryKey"]
        except (TypeError, AttributeError) as err:
            logger.warning("Error occurred when query NCBI. Error:\n%s\nCaching this batches accessions", err)
            connection_err = True

    return epost_webenv, epost_query_key, invalid_err, connection_err


def get_protein_accession(record: SeqRecord, acc_to_retrieve: Optional[list] = None) -> Optional[str]:
    """Extract NCBI Protein accession from SeqRecord."""
    patterns = [
        r"\D{3}\d+\.\d+",
        r"\D{2}_\d+\.\d+",
        r"\D\d+\.\d+",
        r"\D\d{2}\D+\d+\.\d+",
        r"\D\d+\.\d+",
        r"\D\d{2}\D+\d+",
        r"\D\d+"
    ]

    for pattern in patterns:
        match = re.match(pattern, record.id)
        if match:
            return match.group()

    if acc_to_retrieve:
        for acc in acc_to_retrieve:
            if acc in record.id:
                return acc

    parts = record.id.split("|")
    if len(parts) == 3:
        if parts[0] == "sp" and parts[1]:
            return parts[1]
        return parts[2] if not parts[1] else parts[1]

    return None
