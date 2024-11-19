#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Parse taxonomy data from NCBI and CAZy to replace when multiple taxa are retrieved from CAZy"""


import logging

from argparse import Namespace
from http.client import IncompleteRead

from Bio.Entrez.Parser import NotXMLError, CorruptedXMLError
from saintBioutils.genbank import entrez_retry
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm

from Bio import Entrez


class NcbiProtein:
    """Represent the data extract from a NCBI protein record"""
    def __init__(self):
        self.protein_id = None
        self.genus = None
        self.species = None
        self.kingdom = None

    def get_ncbi_data(self, record):
        """extract the data from the record"""
        self.protein_id = record['GBSeq_accession-version']
        self.genus = record['GBSeq_organism'].split(maxsplit=1)[0]
        self.species = " ".join(record['GBSeq_organism'].split()[1:])
        self.kingdom = record['GBSeq_taxonomy'].split(';')[0]


logger = logging.getLogger(__name__)


def get_ncbi_tax(epost_results, args: Namespace) -> None | list[NcbiProtein]:
    """Parse the ePost output from Entrez

    :param epost_results: Entrez ePost output
    :param args: cmd-line args parser
    """

    # Retrieve web environment and query key from Entrez epost
    epost_webenv = epost_results["WebEnv"]
    epost_query_key = epost_results["QueryKey"]

    try:
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            db="Protein",
            query_key=epost_query_key,
            WebEnv=epost_webenv,
            retmode="xml",
        ) as record_handle:
            protein_records = Entrez.read(record_handle, validate=False)

    # if no record is returned from call to Entrez
    except (TypeError, AttributeError, RuntimeError,  NotXMLError, IncompleteRead, CorruptedXMLError) as error:
        logger.error(
            "Entrez failed to retireve accession numbers.\n%s\n"
            "Will use the last taxonomy to returned from CAZy", error
        )
        return

    ncbi_data = []
    for record in tqdm(protein_records, desc="Retrieving organism from NCBI"):
        protein = NcbiProtein()
        protein.get_ncbi_data(record)
        ncbi_data.append(protein)

    return ncbi_data
