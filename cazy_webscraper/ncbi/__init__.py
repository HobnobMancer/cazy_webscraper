#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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
"""Parse data retrieved from CAZy and build a dict of data matching the user's criteria."""


import logging

from saintBioutils.genbank import entrez_retry
from tqdm import tqdm

from Bio import Entrez


def replace_multiple_tax(cazy_data, args):
    """
    Identify GenBank accessions which have multiple source organisms listedi in CAZy. Replace with
    the latest source organism from NCBI.

    :param cazy_data: dict of CAZy data
    :param multiple_taxa_logger: Path, log file where proteins withe multiple taxa are listed
    :param args: cmd-line args parser

    Return dict updated cazy_data dict
    """
    logger = logging.getLogger(__name__)

    Entrez.email = args.email

    multi_taxa_gbk = {}

    for genbank_accession in tqdm(cazy_data, total=len(list(cazy_data.keys())), desc='Searching for multiple taxa annotations'):
        gbk_organisms = cazy_data[genbank_accession]["organism"]
        if len(gbk_organisms) > 1:
            multi_taxa_gbk[genbank_accession] = {"families": cazy_data[genbank_accession]["families"]}
            logger.warning(
                f"{genbank_accession} annotated with multiple taxa:\n"
                f"{cazy_data[genbank_accession]['organism']}"
            )
    
    logger.warning(
        f"{len(list(multi_taxa_gbk.keys()))} proteins found with multiple source organisms in CAZy\n"
        "Querying NCBI to retrieve the latest source organisms"
    )

    if len(list(multi_taxa_gbk.keys())) == 0:
        return cazy_data

    id_post_list = str(",".join(list(multi_taxa_gbk.keys())))
    try:
        epost_results = Entrez.read(
            entrez_retry(
                args.retries, Entrez.epost, "Protein", id=id_post_list
            )
        )
    # if no record is returned from call to Entrez
    except (TypeError, AttributeError):
        logger.error(
                f"Entrez failed to post assembly IDs.\n"
                "Exiting retrieval of accession numbers, and returning null value 'NA'"
        )
    # except RuntimeError:

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
    except (TypeError, AttributeError) as error:
        logger.error(
            f"Entrez failed to retireve accession numbers."
            "Exiting retrieval of accession numbers, and returning null value 'NA'"
        )
    
    for protein in tqdm(protein_records, desc="Retrieving organism from NCBI"):
        # retrieve NCBI taxonomy data
        accession = protein['GBSeq_accession-version']
        organism = protein['GBSeq_organism']
        kingdom = protein['GBSeq_taxonomy'].split(';')[0]

        # retrieve CAZy taxonomy data
        cazy_kingdom = cazy_data[accession]["kingdom"]
        cazy_organisms = cazy_data[accession]["organism"]
        
        try:
            cazy_data[accession] = {
                "kingdom": (kingdom,),
                "organism": (organism,),
                "families": multi_taxa_gbk[accession]["families"],
            }

            # log the difference
            cazy_organism_str = ','.join(cazy_organisms)
            multiple_taxa_logger.warning(
                f"{accession}\t{cazy_kingdom}: {cazy_organism_str}\t{kingdom}: {organism}"
            )

        except KeyError:
            err = f'GenBank accession {accession} retrieved from NCBI, but it is not present in CAZy'
            logger.error(err)

            file_content = f"{accession}\t{cazy_kingdom}: {cazy_organism_str}\t{err}"
            with open(multiple_taxa_logger, 'a') as fh:
                fh.write(file_content)

    return cazy_data


def get_ncbi_tax(epost_results, cazy_data, replaced_taxa_logger, args):
    """Parse the ePost output from Entrez and add the NCBI tax classifications to the CAZy data
    
    :param epost_results: Entrez ePost output
    :param cazy_data: dict, data retrieved from CAZy
    :param args: cmd-line args parser
    :param replaced_taxa_logger: logger, used for logging GenBank accessions with multiple taxa in CAZy
    
    Return cazy_data (dict)
    """
    logger = logging.getLogger(__name__)
    
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
    except (TypeError, AttributeError) as error:
        logger.error(
            f"Entrez failed to retireve accession numbers."
            "Exiting retrieval of accession numbers, and returning null value 'NA'"
        )
    
    for protein in tqdm(protein_records, desc="Retrieving organism from NCBI"):
        # retrieve NCBI taxonomy data
        accession = protein['GBSeq_accession-version']
        organism = protein['GBSeq_organism']
        kingdom = protein['GBSeq_taxonomy'].split(';')[0]

        # retrieve CAZy taxonomy data
        cazy_kingdom = cazy_data[accession]["kingdom"]
        cazy_organisms = cazy_data[accession]["organism"]
    
        cazy_kingdom_str = ",".join(cazy_kingdom)
        cazy_organism_str = ','.join(cazy_organisms)
        
        try:
            cazy_data[accession]['kingdom'] = kingdom
            cazy_data[accession]['organism'] = organism

            # log the difference
            replaced_taxa_logger.warning(
               f"{accession}\t{cazy_kingdom_str}: {cazy_organism_str}\t{kingdom}: {organism}"
            )

        except KeyError:
            err = f'GenBank accession {accession} retrieved from NCBI, but it is not present in CAZy'
            logger.error(err)
            
            replaced_taxa_logger.warning(
               f"{accession}\t{cazy_kingdom_str}: {cazy_organism_str}\t{err}"
            )

    return cazy_data


def select_first_organism(cazy_data, gbk_accessions, replaced_taxa_logger):
    """Select the first organism listed for each GenBank accession
    
    :param cazy_data: dict of data retrieved from CAZy
    :param gbk_accessions: list of GenBank accessions
    :param replaced_taxa_logger: logger, used for logging GenBank accessions with multiple taxa in CAZy
    
    Return cazy_data (dict)
    """
    for accession in tqdm(gbk_accessions, desc='Selecting the first retrieved organism'):
        cazy_kingdoms_str = ",".join(list(cazy_data[accession]["kingdom"]))
        cazy_organisms_str = ",".join(list(cazy_data[accession]["organism"]))

        cazy_kingdom = list(cazy_data[accession]["kingdom"])[0]
        cazy_organism = list(cazy_data[accession]["organism"])[0]

        cazy_data[accession]["kingdom"] = {cazy_kingdom}
        cazy_data[accession]["organism"] = {cazy_organism}

        # log the data that was replaced and the data it was replaced with
        replaced_taxa_logger.warning(
            f"{accession}\t{cazy_kingdoms_str}: {cazy_organisms_str}\t{cazy_kingdom}: {cazy_organism}"
        )
    
    return cazy_data
