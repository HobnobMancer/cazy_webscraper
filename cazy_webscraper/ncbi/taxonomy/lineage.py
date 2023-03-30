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

from saintBioutils.genbank import entrez_retry
from tqdm import tqdm

from Bio import Entrez


def fetch_lineages(tax_ids, query_key, web_env, args):
    """Fetch lineage data from NCBI Taxonomy database.

    :param tax_ids: list of NCBI tax ids
    :param query_key: str, query key from ncbi.entrez.epost
    :param web_env: str, web environment from ncbi.entrez.epost
    :param args: CLI args parser

    Return dict of lineage data {ncbi_tax_id: {rank: str/lineage}}
        or None if connection fails
    """
    lineage_dict = {}  # {ncbi_tax_id: {rank: str/lineage}}

    try:
        with entrez_retry(
            10,
            Entrez.efetch,
            db="Taxonomy",
            query_key=qk,
            WebEnv=we,
        ) as handle:
            batched_tax_records = Entrez.read(handle, validate=False)

    except (TypeError, AttributeError) as err:
        logger.warning(f"Failed to fetch tax record from NCBI tax for id '{tax_id}'':\n{err}")
        return

    for record in batched_tax_records:
        record_id = record['TaxId']
        if record_id not in tax_ids:
            continue

        # set all lineage data to None
        kingdom, phylum, tax_class, order, family, genus, species, strain = None, None, None, None, None, None, None, None

        for i in record['LineageEx']:
            rank = i['Rank']

            if rank == 'superkingdom':
                kingdom = i['ScientificName']

            elif rank == 'phylum':
                phylum = i['ScientificName']

            elif rank == 'class':
                tax_class = i['ScientificName']

            elif rank == 'order':
                order = i['ScientificName']

            elif rank == 'family':
                family = i['ScientificName']

            elif rank == 'genus':
                genus = i['ScientificName']

            elif rank == 'species' or 'species group':
                species = i['ScientificName']

            elif rank == 'serotype' or 'strain':
                strain = i['ScientificName']

        scientific_name = record['ScientificName']

        if genus is not None:
            # drop genus from species name
            if species is not None:
                species = species.replace(genus, "").strip()

            # extract strain from scientific name if not retrieved as rank
            if species is not None and strain is None:
                strain = scientific_name.replace(f"{genus} {species}", "").strip()

            # extract species from the scientific name if not retrieved as rank
            elif species is None:
                species = scientific_name.replace(genus, "").strip()

        lineage_dict[record_id] = {
            'kingdom': kingdom,
            'phylum': phylum,
            'class': tax_class,
            'order': order,
            'family': family,
            'genus': genus,
            'species': species,
            'strain': strain,
        }

    return lineage_dict
