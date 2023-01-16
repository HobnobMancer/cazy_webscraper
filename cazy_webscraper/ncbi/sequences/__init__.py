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
import re

from http.client import IncompleteRead

from Bio import Entrez, SeqIO
from Bio.Entrez.Parser import NotXMLError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Entrez.Parser import NotXMLError
from saintBioutils.genbank import entrez_retry
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm
 

def post_accessions_to_entrez(batch, args):
    """Post NCBI protein accessions to NCBI via Entrez, and capture error message if one returned.

    :param batch: list of genbank accessions
    :param args: CLI args parser

    Return Entrez ePost web environment and query key.
    """
    logger = logging.getLogger(__name__)
    success, epost_result = None, None
    
    try:
        epost_result = Entrez.read(
            entrez_retry(
                args.retries,
                Entrez.epost,
                db="Protein",
                id=",".join(batch),
            ),
            validate=False,
        )

    except RuntimeError as err:
        if repr(err).startswith("RuntimeError('Some IDs have invalid value and were omitted.") or repr(err).startswith("RuntimeError('Empty ID list; Nothing to store')"):
            
            if len(batch) == 1:
                logger.warning(
                    f"Accession '{batch[0]}' not listed in NCBI. Querying NCBI returns the error:\n"
                    f"{repr(err)}\n"
                    f"Not retrieving seq for '{batch[0]}'"
                )
                success = "Invalid ID"
            
            else:
                logger.warning(
                    "Batch contains an accession no longer listed in NCBI\n"
                    f"Querying NCBI returns the error:\n{err}\n"
                    "Will identify invalid ID(s) later"
                )
                success = "Contains invalid ID"

        else:  # unrecognised Runtime error
            if len(batch) == 1:
                logger.warning(
                    f"Runtime error occured when querying NCBI by accession '{batch[0]}'\n"
                    f"Error returned:\n{err}\n"
                    "Interrupting error as recognition of an invalid ID. Therefore,\n"
                    f"not adding seq for '{batch[0]} to the local CAZyme db'"
                )
                success = "Invalid ID"

            else:
                logger.warning(
                    f"Runtime error occured when querying NCBI with batch of gbk accessions\n"
                    f"Error returned:\n{err}\n"
                    "Interrupting error as batch containing invalid ID.\n"
                    "Will identify invalid ID(s) later"
                )
                success = "Contains invalid ID"

    except IncompleteRead as err:
        logger.warning(
            "IncompleteRead error raised when querying NCBI:\n"
            f"{err}\n"
            "Will reattempt NCBI query later"
        )
        success = "Failed connection"

    except NotXMLError as err:
        logger.warning(
            "NotXMLError raised when querying NCBI:\n"
            f"{err}\n"
            "Will reattempt NCBI query later"
        )
        success = "Failed connection"

    except (TypeError, AttributeError) as err:  # if no record is returned from call to Entrez
        logger.warning(
            f"Error occurenced when batch posting IDs to NCBI\n"
            "Error retrieved:\n"
            f"{repr(err)}\n"
            "Will retry batch later"
        )
        success = "Failed connection"

    except Exception as err:  # if no record is returned from call to Entrez
        logger.warning(
            f"Error occurenced when posting IDs to NCBI\n"
            "Error retrieved:\n"
            f"{repr(err)}\n"
            "Will retry batch later"
        )
        success = "Failed connection"

    # retrieve the web environment and query key from the Entrez post
    try:
        epost_webenv = epost_result["WebEnv"]
        epost_query_key = epost_result["QueryKey"]
        success = "Complete"
    except (TypeError, AttributeError):
        epost_webenv, epost_query_key = None, None
        pass  # raised when could not epost failed

    return epost_webenv, epost_query_key, success


def fetch_ncbi_seqs(seq_records, epost_webenv, epost_query_key, acc_to_retrieve, args):
    """Retrieve protein sequences from NCBI from ePost of protein v.accs.

    :param seq_records: list of Bio.SeqRecords
    :param epost_websenv: Entrez ePost webenvironment
    :param epost_query_key: Entrez ePost query key
    :param acc_to_retrieve: set of NCBI protein version accessions to retrieve seqs for
    :param args: CLI args parser

    Return updated list of SeqRecords and string marking success/failure or seq retrieval.
    """
    logger = logging.getLogger(__name__)
    success, successful_accessions = None, set()

    try:
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            db="Protein",
            query_key=epost_query_key,
            WebEnv=epost_webenv,
            rettype="fasta",
            retmode="text",
        ) as seq_handle:
            for record in SeqIO.parse(seq_handle, "fasta"):
                retrieved_accession = get_protein_accession(record, acc_to_retrieve=acc_to_retrieve)

                if retrieved_accession not in acc_to_retrieve:
                    retrieved_accession = None

                if retrieved_accession is None:
                    logger.error(
                        "Could not retrieve a NCBI protein version accession matching\n"
                        f"an accession from the local database from the record id '{record.id}'\n"
                        "The sequence from this record will not be added to the db"
                    )
                    continue
    
                seq_records.append(record)
                successful_accessions.add(retrieved_accession)

    except RuntimeError as err:
        if repr(err).startswith("RuntimeError('Some IDs have invalid value and were omitted.") or repr(err).startswith("RuntimeError('Empty ID list; Nothing to store')"):
            
            if len(batch) == 1:
                logger.warning(
                    f"Accession '{batch[0]}' not listed in NCBI. Querying NCBI returns the error:\n"
                    f"{repr(err)}\n"
                    f"Not retrieving seq for '{batch[0]}'"
                )
                success = "Invalid ID"
            
            else:
                logger.warning(
                    "Batch contains an accession no longer listed in NCBI\n"
                    f"Querying NCBI returns the error:\n{err}\n"
                    "Will identify invalid ID(s) later"
                )
                success = "Contains invalid ID"

        else:  # unrecognised Runtime error
            if len(batch) == 1:
                logger.warning(
                    f"Runtime error occured when fetching record from NCBI for accession '{batch[0]}'\n"
                    f"Error returned:\n{err}\n"
                    "Interrupting error as recognition of an invalid ID. Therefore,\n"
                    f"not adding seq for '{batch[0]} to the local CAZyme db'"
                )
                success = "Invalid ID"

            else:
                logger.warning(
                    f"Runtime error occured when fetching records from NCBI\n"
                    f"Error returned:\n{err}\n"
                    "Interrupting error as batch containing invalid ID.\n"
                    "Will identify invalid ID(s) later"
                )
                success = "Contains invalid ID"

    except IncompleteRead as err:
        logger.warning(
            "IncompleteRead error raised when fetching record from NCBI:\n"
            f"{err}\n"
            "Will reattempt NCBI query later"
        )
        success = "Failed connection"

    except NotXMLError as err:
        logger.warning(
            "NotXMLError raised when fetching record(s) from NCBI:\n"
            f"{err}\n"
            "Will reattempt NCBI query later"
        )
        success = "Failed connection"

    except (TypeError, AttributeError) as err:  # if no record is returned from call to Entrez
        logger.warning(
            f"Error occurenced when fetching records from NCBI\n"
            "Error retrieved:\n"
            f"{repr(err)}\n"
            "Will retry batch later"
        )
        success = "Failed connection"

    except Exception as err:  # if no record is returned from call to Entrez
        logger.warning(
            f"Error occurenced when fetching sequences from NCBI\n"
            "Error retrieved:\n"
            f"{repr(err)}\n"
            "Will retry batch later"
        )
        success = "Failed connection"

    return seq_records, success, list(successful_accessions)


def get_protein_accession(record, acc_to_retrieve=None):
    """Extract NCBI Protein accession from SeqRecord.

    :param record: BioPython SeqRecord obj
    :param acc_to_retrieve: list of acc to retrieve seqs for

    Return str (accession) or None:
    """
    logger = logging.getLogger(__name__)
    retrieved_accession = None

    # check if multiple items returned in ID
    # try re search for accession in string
    try:
        retrieved_accession = re.match(r"\D{3}\d+\.\d+", record.id).group()

    except AttributeError:
        try:
            retrieved_accession = re.match(r"\D{2}_\d+\.\d+", record.id).group()

        except AttributeError:
            if acc_to_retrieve is not None:
                for acc in acc_to_retrieve:
                    if record.id.find(acc) != -1:
                        retrieved_accession = acc

    if retrieved_accession is None:
        # Accession may be a UniProt entry name or ID
        # e.g. 'prf||2109195A'or 'sp|B2FSW8.1|EALGL_STRMK'
        if len(record.id.split("|")) == 3:

            if record.id.startswith("sp"):
                if len(record.id.split("|")[1]) == 0:
                    retrieved_accession = record.id.split("|")[2]
                else:
                    retrieved_accession = record.id.split("|")[1]
            
            elif len(record.id.split("|")[1]) == 0:
                retrieved_accession = record.id.split("|")[2]

        if retrieved_accession is None:  # Accept Uniprot-style accessions which ar sometimes used
            try:
                retrieved_accession = re.match(r"\D\d+\.\d+", record.id).group()
            except AttributeError:
                try:
                    retrieved_accession = re.match(r"\D\d{2}\D+\d+\.\d+", record.id).group()
                except AttributeError:
                    try:
                        retrieved_accession = re.match(r"\D\d+\.\d+", record.id).group()
                    except AttributeError:
                        try:
                            if record.id == re.match(r"\D\d{2}\D+\d+", record.id).group():
                                retrieved_accession = re.match(r"\D\d{2}\D+\d+", record.id).group()
                        except AttributeError:
                            try:
                                if record.id == re.match(r"\D\d+", record.id).group():
                                    retrieved_accession = re.match(r"\D\d+", record.id).group()
                            except AttributeError:
                                retrieved_accession = None
    
    return retrieved_accession
