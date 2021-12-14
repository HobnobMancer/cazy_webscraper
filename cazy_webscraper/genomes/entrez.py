#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
"""Functions that call to and parse data from NCBI.Entrez."""


import logging

from Bio import Entrez
from saintBioutils.genbank import entrez_retry
from tqdm import tqdm


def get_linked_nucleotide_record_ids(batch_post, args):
    """Use Entrez.elink to get the ID of the NCBI.Nucleotide db record linked to the NCBI.Protein
    db record.
    
    In the resulting dict, the protein record ID is used to group nucleotide records that are 
    linked to the same protein record, and from which it will need to be determined which is
    the longest record, so that only one (and the longest) Nucleotide record is parsed for this protein
    not all retrieved Nucleotide records.

    :param batch_post: Entrez dict containing WebEnv and query_key from batch posting protein accessions
    :param args: cmd-line args parser

    Return dict {protein_accession: {nucleotide records IDs}}
    Return None if cannot retrieve data from NCBI
    """
    logger = logging.getLogger(__name__)

    nucleotide_ids = {}  # {protein_record_id: {nucleotide_ids}}
    # used for identifying which nucleotide record to retrieve protein accessions from
    with entrez_retry(
        args.retries,
        Entrez.elink,
        dbfrom="Protein",
        db="nuccore",
        query_key=batch_post['QueryKey'],
        WebEnv=batch_post['WebEnv'],
        linkname='protein_nuccore',
    ) as handle:
        try:
            batch_nuccore = Entrez.read(handle)
        except Exception as err:
            logger.warning(f"Failed Entrez connection: {err}\nNo nucletoide IDs retrieved for batch query")
            return

    for record in tqdm(batch_nuccore, desc="Retrieving Nucleotide db records IDs from nuccore"):
        protein_record_id = record['IdList'][0]
        record_nucleotide_ids = set()
    
        # Linked db records are contained in the 'LinkSetDb' field
        # multiple linked records may be retrieved
        
        for link_dict in record['LinkSetDb']:
            # record['LinkSetDb'] = [{'Link': [{'Id': '1127815138'}], 'DbTo': 'nuccore', 'LinkName': 'protein_nuccore'}]
            # link_list = {'Link': [{'Id': '1127815138'}], 'DbTo': 'nuccore', 'LinkName': 'protein_nuccore'}
            links = link_dict['Link']
            for link in links:
                # links = [{'Id': '1127815138'}]
                # link = {'Id': '1127815138'}
                nucleotide_id = link['Id']

                record_nucleotide_ids.add(nucleotide_id)

        try:
            existing_ids = nucleotide_ids[protein_record_id]
            all_ids = existing_ids.union(record_nucleotide_ids)
            nucleotide_ids[protein_record_id] = all_ids
        except KeyError:
            nucleotide_ids[protein_record_id] = record_nucleotide_ids
    
    return nucleotide_ids


def link_nucleotide_ids_individually(accessions_to_parse, no_accession_logger, args):
    """Use Entrez.elink to get the ID of the NCBI.Nucleotide db record linked to the NCBI.Protein
    db record.
    
    In the resulting dict, the protein record ID is used to group nucleotide records that are 
    linked to the same protein record, and from which it will need to be determined which is
    the longest record, so that only one (and the longest) Nucleotide record is parsed for this protein
    not all retrieved Nucleotide records.

    :param accessions_to_parse: list of protein accessions
    :param no_accession_logger: path to log file for logging protein accessions for which no 
        linked nucleotide record was retrieved
    :param args: cmd-line args parser

    Return dict {protein_accession: {nucleotide records IDs}}
    no_nucleotides: set of protein accessions for which no nucleotide record could be retrieved
    """
    logger = logging.getLogger(__name__)

    nucleotide_ids = {}  # {protein_record_id: {nucleotide_ids}}
    no_nucleotides = set()  # protein accessions for which no nucleotide record could be retrieved

    for protein_accession in tqdm(accessions_to_parse, "Retrieving Nucleotide IDs individually"):
    # used for identifying which nucleotide record to retrieve protein accessions from
        with entrez_retry(
            args.retries,
            Entrez.elink,
            dbfrom="Protein",
            db="nuccore",
            id=protein_accession,
        ) as handle:
            try:
                batch_nuccore = Entrez.read(handle)
            except Exception as err:
                logger.warning(
                    f"Could not retrieve linked Nucletoide record for {protein_accession}"
                )
                with open(no_accession_logger, "a") as fh:
                    fh.write(
                        f"Could not retrieve linked Nucletoide record for {protein_accession}\t"
                        f"Error: {err}"
                    )
                no_nucleotides.add(protein_accession)
                continue

        for record in tqdm(batch_nuccore, desc="Retrieving Nucleotide db records IDs from nuccore"):
            protein_record_id = record['IdList'][0]
            record_nucleotide_ids = set()
        
            # Linked db records are contained in the 'LinkSetDb' field
            # multiple linked records may be retrieved

            if len(record['LinkSetDb']) == 0:
                logger.warning(
                    f"Could not retrieve linked Nucletoide record for {protein_accession}"
                )
                err = "No data containiend in 'LinkSetDb' field from eLink output"
                with open(no_accession_logger, "a") as fh:
                    fh.write(
                        f"Could not retrieve linked Nucletoide record for {protein_accession}\t"
                        f"Error: {err}"
                    )
                no_nucleotides.add(protein_accession)
                continue
            
            for link_dict in record['LinkSetDb']:
                # record['LinkSetDb'] = [{'Link': [{'Id': '1127815138'}], 'DbTo': 'nuccore', 'LinkName': 'protein_nuccore'}]
                # link_list = {'Link': [{'Id': '1127815138'}], 'DbTo': 'nuccore', 'LinkName': 'protein_nuccore'}
                links = link_dict['Link']
                for link in links:
                    # links = [{'Id': '1127815138'}]
                    # link = {'Id': '1127815138'}
                    nucleotide_id = link['Id']

                    record_nucleotide_ids.add(nucleotide_id)

            try:
                existing_ids = nucleotide_ids[protein_record_id]
                all_ids = existing_ids.union(record_nucleotide_ids)
                nucleotide_ids[protein_record_id] = all_ids
            except KeyError:
                nucleotide_ids[protein_record_id] = record_nucleotide_ids
    
    return nucleotide_ids, no_nucleotides


def extract_protein_accessions(single_nucleotide_ids, retrieved_proteins, gbk_accessions, args):
    """Retrieve and parse Nucleotide db records, for NCBI.Protein records from which only
    one NCBI.Nucleotide db record ID was retrieved.
    
    :param retrieved_proteins: dict, {protein_accession: nucleotide record accession}
    :param single_nucleotide_ids: list of nucloetide record IDs
    :param gbk_accessions: list of protein GenBank accessions from the local CAZyme database
    :param args: cmd-line args parser
    
    Return retrieved_proteins (dict)
        newly_retrieved_proteins: set of CAZyme protein accessions retrieved from parsed records
        bool: True if successful Entrez connection, False is connection fails
    """
    logger = logging.getLogger(__name__)
    newly_retrieved_proteins = set()

    batch_query_ids = ",".join(list(single_nucleotide_ids))
    with entrez_retry(
        args.retries, Entrez.epost, "Nucleotide", id=batch_query_ids,
    ) as handle:
        batch_post = Entrez.read(handle)

    with entrez_retry(
        args.retries,
        Entrez.efetch,
        db="Nucleotide",
        query_key=batch_post['QueryKey'],
        WebEnv=batch_post['WebEnv'],
        retmode="xml",
    ) as handle:
        try:
            batch_nucleotide = Entrez.read(handle)
        except Exception as err:
            logger.warning(
                f"Failed Entrez connection for fetching Nucleotide records: {err}"
            )
            return retrieved_proteins, newly_retrieved_proteins, False
    
    for record in tqdm(batch_nucleotide, desc="Extracting data from Nucleotide records"):
        nucleotide_accession = record['GBSeq_accession-version']

        # retrieve protein accessions of proteins features in the nucletide record
        for feature_dict in record['GBSeq_feature-table']:
            # feature-table contains a list of features, one feature is one feature_dict

            for feature_qual in feature_dict['GBFeature_quals']:
                # feature_quals contains a list of dicts, one dict is feature_qual
                # looking for dict containing protein accession (protein_id)
                # e.g. {'GBQualifier_name': 'protein_id', 'GBQualifier_value': 'APS93952.1'}
                if feature_qual['GBQualifier_name'] == 'protein_id':
                    protein_accession = feature_qual['GBQualifier_value']

                    if protein_accession in gbk_accessions:
                        # protein is in the local CAZyme database
                        try:
                            retrieved_proteins[protein_accession].add(nucleotide_accession)
                        except KeyError:
                            retrieved_proteins[protein_accession] = {nucleotide_accession}
                        
                        newly_retrieved_proteins.add(protein_accession)

    return retrieved_proteins, newly_retrieved_proteins, True


def extract_protein_accessions_individually(single_nucleotide_ids, retrieved_proteins, gbk_accessions, args):
    """Retrieve and parse Nucleotide db records, for NCBI.Protein records from which only
    one NCBI.Nucleotide db record ID was retrieved.
    
    :param retrieved_proteins: dict, {protein_accession: nucleotide record accession}
    :param single_nucleotide_ids: list of nucloetide record IDs
    :param gbk_accessions: list of protein GenBank accessions from the local CAZyme database
    :param args: cmd-line args parser
    
    Return retrieved_proteins (dict)
        newly_retrieved_proteins: set of CAZyme protein accessions retrieved from parsed records
    """
    logger = logging.getLogger(__name__)

    newly_retrieved_proteins = set()

    for nucleotide_id in tqdm(single_nucleotide_ids, desc="Parsing nucelotide records individually"):
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            db="Nucleotide",
            id=nucleotide_id,
            retmode="xml",
        ) as handle:
            try:
                batch_nucleotide = Entrez.read(handle)
            except Exception as err:
                logger.warning(
                    f"Failed Entrez connection for fetching Nucleotide records: {err}"
                )
                pass
    
        for record in tqdm(batch_nucleotide, desc="Extracting data from Nucleotide records"):
            nucleotide_accession = record['GBSeq_accession-version']

            # retrieve protein accessions of proteins features in the nucletide record
            for feature_dict in record['GBSeq_feature-table']:
                # feature-table contains a list of features, one feature is one feature_dict

                for feature_qual in feature_dict['GBFeature_quals']:
                    # feature_quals contains a list of dicts, one dict is feature_qual
                    # looking for dict containing protein accession (protein_id)
                    # e.g. {'GBQualifier_name': 'protein_id', 'GBQualifier_value': 'APS93952.1'}
                    if feature_qual['GBQualifier_name'] == 'protein_id':
                        protein_accession = feature_qual['GBQualifier_value']

                        if protein_accession in gbk_accessions:
                            # protein is in the local CAZyme database
                            try:
                                retrieved_proteins[protein_accession].add(nucleotide_accession)
                            except KeyError:
                                retrieved_proteins[protein_accession] = {nucleotide_accession}
                            
                            newly_retrieved_proteins.add(protein_accession)

    return retrieved_proteins, newly_retrieved_proteins


def parse_longest_record(nucleotide_record_ids, retrieved_proteins, gbk_accessions, args):
    """Identify the longest NCBI.Nucleotide record, and extract Protein GenBank accessions
    
    :param nucleotide_record_ids: set, NCBI.Nucleotide records IDs retrieved for one Protein record
    :param retrieved_proteins: dict, {protein_accession: nucleotide record accession}
    :param gbk_accessions: list of protein GenBank accessions from the local CAZyme database
    :param args: cmd-line args parser
    
    Return retrieved_proteins (dict)
        newly_retrieved_proteins: set of CAZyme protein accessions retrieved from parsed records
        bool: True if successful Entrez connection, False is connection fails
    """
    logger = logging.getLogger(__name__)
    newly_retrieved_proteins = set()

    batch_query_ids = ",".join(list(nucleotide_record_ids))
    with entrez_retry(
        args.retries, Entrez.epost, "Nucleotide", id=batch_query_ids,
    ) as handle:
        batch_post = Entrez.read(handle)
    print("posted")
    with entrez_retry(
        args.retries,
        Entrez.efetch,
        db="Nucleotide",
        query_key=batch_post['QueryKey'],
        WebEnv=batch_post['WebEnv'],
        retmode="xml",
    ) as handle:
        try:
            batch_nucleotide = Entrez.read(handle)
        except Exception as err:
            logger.warning(
                f"Failed Entrez connection for fetching Nucleotide records: {err}"
            )
            return retrieved_proteins, newly_retrieved_proteins, False
    print("fetched")
    record_lengths = {}  # {Nucleotide record accession: {len: Number of features (int), record: record}
    # longest (most features) record interpretted as the most complete record
    
    for record in tqdm(batch_nucleotide, desc="Selecting longest Nucleotide record"):
        nucleotide_accession = record['GBSeq_accession-version']

        number_of_features = 0

        for feature_dict in record['GBSeq_feature-table']:
            for feature_qual in feature_dict['GBFeature_quals']:
                if feature_qual['GBQualifier_name'] == 'protein_id':
                    number_of_features += 1
        
        record_lengths[nucleotide_accession] = {'length': number_of_features, 'record': record}
    
    list_of_lengths = [acc['length'] for acc in list(record_lengths.keys())]
    list_of_lengths.sort(reverse=True)
    longest_length = list_of_lengths[0]

    for nucleotide_accession in record_lengths:
        if record_lengths[nucleotide_accession]['length'] == longest_length:
            # found the longest record
            # extract protein accessions for CAZymes in the local CAZyme database
            # method explained in extract_protein_accessions()
            for feature_qual in feature_dict['GBFeature_quals']:
                if feature_qual['GBQualifier_name'] == 'protein_id':
                    protein_accession = feature_qual['GBQualifier_value']
                    if protein_accession in gbk_accessions:
                        try:
                            retrieved_proteins[protein_accession].add(nucleotide_accession)
                        except KeyError:
                            retrieved_proteins[protein_accession] = {nucleotide_accession}

                        newly_retrieved_proteins.add(protein_accession)
            break
    
    return retrieved_proteins, newly_retrieved_proteins, True


def parse_longest_record_individually(nucleotide_record_ids, retrieved_proteins, gbk_accessions, args):
    """Identify the longest NCBI.Nucleotide record, and extract Protein GenBank accessions
    
    :param nucleotide_record_ids: set, NCBI.Nucleotide records IDs retrieved for one Protein record
    :param retrieved_proteins: dict, {protein_accession: nucleotide record accession}
    :param gbk_accessions: list of protein GenBank accessions from the local CAZyme database
    :param args: cmd-line args parser
    
    Return retrieved_proteins (dict)
        newly_retrieved_proteins: set of CAZyme protein accessions retrieved from parsed records
        bool: True if successful Entrez connection, False is connection fails
    """
    newly_retrieved_proteins = set()

    record_lengths = {}  # {Nucleotide record accession: {len: Number of features (int), record: record}
    # longest (most features) record interpretted as the most complete record

    for nucleotide_id in tqdm(nucleotide_record_ids, desc="Parsing nucleotide records individually"):
        with entrez_retry(
            args.retries,
            Entrez.efetch,
            db="Nucleotide",
            id=nucleotide_id,
            retmode="xml",
        ) as handle:
            try:
                nucleotide_records = Entrez.read(handle)
            except Exception:
                pass

            for record in nucleotide_records:
                nucleotide_accession = record['GBSeq_accession-version']

                number_of_features = 0

                for feature_dict in record['GBSeq_feature-table']:
                    for feature_qual in feature_dict['GBFeature_quals']:
                        if feature_qual['GBQualifier_name'] == 'protein_id':
                            number_of_features += 1
                
                record_lengths[nucleotide_accession] = {'length': number_of_features, 'record': record}
            
    list_of_lengths = [acc['length'] for acc in list(record_lengths.keys())]
    list_of_lengths.sort(reverse=True)
    longest_length = list_of_lengths[0]

    for nucleotide_accession in record_lengths:
        if record_lengths[nucleotide_accession]['length'] == longest_length:
            # found the longest record
            # extract protein accessions for CAZymes in the local CAZyme database
            # method explained in extract_protein_accessions()
            for feature_qual in feature_dict['GBFeature_quals']:
                if feature_qual['GBQualifier_name'] == 'protein_id':
                    protein_accession = feature_qual['GBQualifier_value']
                    if protein_accession in gbk_accessions:
                        try:
                            retrieved_proteins[protein_accession].add(nucleotide_accession)
                        except KeyError:
                            retrieved_proteins[protein_accession] = {nucleotide_accession}

                        newly_retrieved_proteins.add(protein_accession)
            break
    
    return retrieved_proteins, newly_retrieved_proteins, True

