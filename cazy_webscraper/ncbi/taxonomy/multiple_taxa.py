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
from saintBioutils.misc import get_chunks_list
from tqdm import tqdm

from Bio import Entrez


def identify_multiple_taxa(cazy_data, multiple_taxa_logger):
    """Identify GenBank accessions in the CAZy data

    :param cazy_data: dict of CAZy data, keyed by GenBank accessions, valued by dict with the keys
        kingdom, organism, families
    :param multiple_taxa_logger: logger, used for logging GenBank accessions
        with multiple taxa in CAZy

    Return list of GenBank accessions
    """
    multiple_taxa_gbk = []

    for genbank_accession in tqdm(
        cazy_data,
        total=len(list(cazy_data.keys())), desc='Searching for multiple taxa annotations',
    ):

        if len(cazy_data[genbank_accession]['taxonomy']) > 1:
            multiple_taxa_gbk.append(genbank_accession)

            for tax_tuple in cazy_data[genbank_accession]['taxonomy']:
                multiple_taxa_logger.warning(
                    f"{genbank_accession}\t{tax_tuple.kingdom}\t{tax_tuple.organism}"
                )
        
        else:
            cazy_data[genbank_accession]['organism'] = list(cazy_data[genbank_accession]['taxonomy'])[0]

    return multiple_taxa_gbk


def replace_multiple_tax(cazy_data, genbank_accessions, replaced_taxa_logger, args, invalid_ids):
    """Identify GenBank accessions which have multiple source organisms listedi in CAZy.
    Replace with the latest source organism from NCBI.

    :param cazy_data: dict of CAZy data
    :param genbank_accessions: list of genbank accessions with multiple taxa in CAZy
    :param replaced_taxa_logger: logger, used for logging GenBank accessions with 
        multiple taxa in CAZy, and the data this replaced and what is replaced by
    :param args: cmd-line args parser
    :param invalid_ids: boolean, potential presence of invalid GenBank accessions
        Set as true when func is called by replace_multiple_tax_with_invalid_ids()

    Return dict updated cazy_data dict and boolean, if multiple taxa were replaced by single taxa
    """
    logger = logging.getLogger(__name__)

    if args.skip_ncbi_tax:
        logger.warning(
            f"Skipping retrieving the latest taxonomy classification from the NCBI Taxonomy db\n"
            "Adding the first tax listed for each protein in the CAZy db"
        )
        cazy_data = select_first_organism(cazy_data, genbank_accessions, replaced_taxa_logger)
        success = True
        return cazy_data, success

    batches = get_chunks_list(genbank_accessions, args.ncbi_batch_size)

    for batch in tqdm(batches, desc=f"Batch retrieving tax info from NCBI. Batch size:{args.ncbi_batch_size}"):

        id_post_list = str(",".join(batch))

        success = False

        try:
            epost_results = Entrez.read(
                entrez_retry(
                    args.retries,
                    Entrez.epost,
                    "Protein",
                    id=id_post_list,
                )
            )
            success = True

        except (TypeError, AttributeError):  # if no record is returned from call to Entrez
            # error not due to the presence of invalid IDs
            logger.error(
                f"Entrez failed to post assembly IDs for this batch.\n"
                "Not retrieving tax data from NCBI for these proteins"
                "Selecting the first organism retrieved from CAZy as the source organism\nProtein accessions:\n"
                f"{batch}"
            )
            # cazy_data, gbk_accessions, replaced_taxa_logger
            cazy_data = select_first_organism(cazy_data, batch, replaced_taxa_logger)
            success = True
            continue

        except RuntimeError:
            logger.warning("Found GenBank accessions in CAZy data that are no longer in NCBI")

            if invalid_ids:
                # replace_multiple_tax was called by replace_multiple_tax_with_invalid_ids
                # return results, don't use recursive programming
                continue

            else:
                # first time replace_multiple_tax was called
                cazy_data, success = replace_multiple_tax_with_invalid_ids(
                    cazy_data,
                    genbank_accessions,
                    replaced_taxa_logger,
                    args,
                )

        if success is False:
            logger.error(
                "Could not retrieve taxonomy data from NCBI for this batch,\n"
                "Using the first source organism retrieved from CAZy for each GenBank accession\n"
                "Protein accessions:\n"
                f"{batch}"
            )

            cazy_data = select_first_organism(cazy_data, genbank_accessions, replaced_taxa_logger)
            success = True

        else:
            logger.info("Parsing data retrieved from NCBI")
            cazy_data = get_ncbi_tax(epost_results, cazy_data, replaced_taxa_logger, args)

    return cazy_data, success


def get_ncbi_tax(epost_results, cazy_data, replaced_taxa_logger, args):
    """Parse the ePost output from Entrez and add the NCBI tax classifications to the CAZy data

    :param epost_results: Entrez ePost output
    :param cazy_data: dict, data retrieved from CAZy
    :param args: cmd-line args parser
    :param replaced_taxa_logger: logger, used for logging GenBank
        accessions with multiple taxa in CAZy

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
        return cazy_data

    for protein in tqdm(protein_records, desc="Retrieving organism from NCBI"):
        # retrieve NCBI taxonomy data
        accession = protein['GBSeq_accession-version']
        organism = protein['GBSeq_organism']
        kingdom = protein['GBSeq_taxonomy'].split(';')[0]

        try:
            cazy_data[accession]['kingdom'] = kingdom
            cazy_data[accession]['organism'] = organism

        except KeyError:
            err = (
                f'GenBank accession {accession} retrieved from NCBI, but it is not present in CAZy'
            )
            logger.error(err)
            continue

        for tax_tuple in list(cazy_data[accession]['taxonomy'])[1:]:
            replaced_taxa_logger.warning(
                f"{accession}\t"
                f"SELECTED: {kingdom} -- {organism}"
                f"\tREPLACED: {tax_tuple.kingdom}: {tax_tuple.organism}"
            )

    return cazy_data


def replace_multiple_tax_with_invalid_ids(cazy_data, gbk_accessions, replaced_taxa_logger, args):
    """Retrieve the latest taxonomy classification for proteins with multiple source organisms.

    If a GenBank accession is not longer present in NCBI,
        select the first organism retrieved from CAZy.

    :param cazy_data: dict of CAZy data
    :param genbank_accessions: list of GenBank accessions
    :param replaced_taxa_logger: logger, log tax data that is replaced and what it is replaced with
    :param args: cmd-line args parser

    Return dict updated cazy_data dict
    """
    # retrieve the first half of the list
    mid_point = int((len(gbk_accessions)/2))

    half_gbk_accs = gbk_accessions[:mid_point]

    cazy_data, success = replace_multiple_tax(
        cazy_data,
        half_gbk_accs,
        replaced_taxa_logger,
        args,
        invalid_ids=True,
    )

    if success:
        # invalid IDs are stored in the second half of the accessions list
        half_gbk_accs = gbk_accessions[mid_point:]

        for accession in tqdm(half_gbk_accs, desc='Retrieving taxonomies individually'):
            cazy_data, success = replace_multiple_tax(
                cazy_data,
                [accession],
                replaced_taxa_logger,
                args,
                invalid_ids=True,
            )

            if success is False:
                cazy_data = select_first_organism(cazy_data, [accession], replaced_taxa_logger)

    else:
        # invalid gbk ID present in the first half of the accessions list
        for accession in tqdm(half_gbk_accs, desc='Retrieving taxonomies individually'):
            cazy_data, success = replace_multiple_tax(
                cazy_data,
                [accession],
                replaced_taxa_logger,
                args,
                invalid_ids=True,
            )

            if success is False:
                cazy_data = select_first_organism(cazy_data, [accession], replaced_taxa_logger)

        # parse the second half of the accessions list
        half_gbk_accs = gbk_accessions[mid_point:]

        cazy_data, success = replace_multiple_tax(cazy_data, half_gbk_accs, True)

        if success is False:
            # invalid gbk ID present in the second half of the accessions list
            for accession in tqdm(half_gbk_accs, desc='Retrieving taxonomies'):
                cazy_data, success = replace_multiple_tax(
                    cazy_data,
                    [accession],
                    replaced_taxa_logger,
                    args,
                    invalid_ids=True,
                )

                if success is False:
                    cazy_data = select_first_organism(
                        cazy_data,
                        [accession],
                        replaced_taxa_logger,
                    )

    return cazy_data, True


def select_first_organism(cazy_data, gbk_accessions, replaced_taxa_logger):
    """Select the first organism listed for each GenBank accession

    :param cazy_data: dict of data retrieved from CAZy
    :param gbk_accessions: list of GenBank accessions
    :param replaced_taxa_logger: logger, used for logging GenBank accessions with
        multiple taxa in CAZy

    Return cazy_data (dict)
    """
    for accession in tqdm(gbk_accessions, desc='Selecting the first retrieved organism'):
        selected_kingdom = list(cazy_data[accession]['taxonomy'])[0].kingdom
        selected_organism = list(cazy_data[accession]['taxonomy'])[0].organism

        for tax_tuple in list(cazy_data[accession]['taxonomy'])[1:]:
            replaced_taxa_logger.warning(
                f"{accession}\t"
                f"SELECTED: {selected_kingdom} -- {selected_organism}"
                f"\tREPLACED: {tax_tuple.kingdom}: {tax_tuple.organism}"
            )

        cazy_data[accession]["kingdom"] = selected_kingdom
        cazy_data[accession]["organism"] = selected_organism

    return cazy_data
