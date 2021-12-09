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
"""Explore the number of GenBank genomes annotated by CAZy."""


import json
import logging
import os
import time

import pandas as pd

from datetime import datetime
from matplotlib import pyplot as plt
from typing import List, Optional

from Bio import Entrez
from saintBioutils.genbank import entrez_retry
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from cazy_webscraper import cazy_scraper
from cazy_webscraper.expand import get_chunks_list
from cazy_webscraper.sql.sql_interface import get_table_dicts
from cazy_webscraper.utilities.parsers import genbank_cov_parser


KINGDOMS = ['Bacteria', 'Eukaryota', 'Archaea', 'Viruses', 'unclassified']


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # Program preparation
    if argv is None:
        parser = genbank_cov_parser.build_parser()
        args = parser.parse_args()
    else:
        parser = genbank_cov_parser.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)

    Entrez.email = args.email

    # check if need to build output dir
    if os.getcwd() != args.output_dir:
        make_output_directory(args.output_dir, args.force, args.nodelete)

    # connect to the local CAZyme database
    connection, logger_name, cache_dir = cazy_scraper.connect_existing_db(args, time_stamp)
    # make cache_dir
    make_output_directory(cache_dir, args.force_cache, args.nodelete_cache)

    no_accession_logger = cache_dir / "logs"
    # make logs dir
    make_output_directory(no_accession_logger, args.force_cache, args.nodelete_cache)
    no_accession_logger = no_accession_logger / f"no_genomic_accession_retrieved_{time_stamp}.log"

    # load Genbank and Kingdom records from the db
    logger.warning("Retrieving Genbanks, Taxs and Kingdoms records from the local CAZyme db")
    genbank_kingdom_dict = get_table_dicts.get_gbk_kingdom_dict(connection)
    logger.warning("Retrieved Genbanks, Taxs and Kingdoms records from the local CAZyme db")

    add_bioproject_id, genomic_assembly_names = get_assebmly_names(genbank_kingdom_dict, no_accession_logger, args)

    output_path = cache_dir / f"genomic_bioproject_ids_{time_stamp}.json"
    with open(output_path, 'w') as fh:
        json.dump(add_bioproject_id, fh)

    output_path = cache_dir / f"genomic_assembly_names_{time_stamp}.json"
    with open(output_path, 'w') as fh:
        json.dump(genomic_assembly_names, fh)
    
    genomic_accession_dict = get_genomic_accessions(genomic_assembly_names, no_accession_logger, args)
    output_path = cache_dir / f"genomic_accession_numbers_{time_stamp}.json"
    with open(output_path, 'w') as fh:
        json.dump(genomic_accession_dict, fh)

    write_out_genomic_accessions(genomic_accession_dict, time_stamp, args)

    ncbi_genomes_totals = get_ncbi_counts(args)

    write_out_genome_coverage(ncbi_genomes_totals, genomic_accession_dict, time_stamp, args)

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    if args.verbose:
        logger.info(
            "Finished calculting the local CAZyme db's coverage of GenBank\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_scraper.VERSION_INFO}\n"
            f"Citation: {cazy_scraper.CITATION_INFO}"
        )
    else:
        print(
            "=====================cazy_webscraper=====================\n"
            "Finished calculting the local CAZyme db's coverage of GenBank\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_scraper.VERSION_INFO}\n"
            f"Citation: {cazy_scraper.CITATION_INFO}"
        )


def get_nucleotide_ids(genbank_kingdom_dict, no_accession_logger, args):
    """Retrieve the NCBI Nucleotide db records ID's containing the GenBank protein accessions.
    
    :param genbank_kingdom_dict: dict of Genbank and Kingdom records from db
        {kingdom: {genus: {species: {protein_accessions}}}
    :param no_accession_logger: Path, path to log file to save protein accessions for which no 
        genomic accession was retrieved
    :param args: cmd-line args parser
    
    Return dict {kingdom: {genus: {species: {nucleotide_id: {protein_accessions},},},},}
    """
    logger = logging.getLogger(__name__)

    nucleotide_accessions_dict = {}
    # {kingdom: {genus: {species: {nucleotide record accessio: {protein_accessions},},},},}

    retrieved_proteins = {}  # {protein_accession: nucleotide record accession}

    # store all retrieved ncucleotide accessions to prevent retrieval of the same assembly record multiple times
    
    for kingdom in tqdm(genbank_kingdom_dict, desc="Retrieving nucleotide record accessions per kingdom"):
        genera = genbank_kingdom_dict[kingdom]

        gbk_organism_dict = {}  # protein_accession: {species: str, genus: str}

        for genus in genera:
            organisms = genera[genus]

            for species in organisms:
                for gbk_acc in organisms[species]:
                    gbk_organism_dict[gbk_acc] = {'species': species, 'genus': genus}
        
        # retrieve all Gbk protein accessions for the given genera
        gbk_accessions = list(gbk_organism_dict.keys())

        # break up the list into a series of smaller lists that can be batched querried
        batch_queries = get_chunks_list(gbk_accessions, args.batch_size)

        failed_gbk_batches = []
        failed_nuc_batches = []

        for batch_query in tqdm(batch_queries, desc=f"Batch querying NCBI for {kingdom}"):
            # removed proteins for which a nucleotide accession has already been retrieved
            cleaned_batch = [_ for _ in batch_query is _ not in list(retrieved_proteins.keys())]
            if len(cleaned_batch) == 0:
                continue

            batch_query_ids = ",".join(cleaned_batch)
            with entrez_retry(
                args.retries, Entrez.epost, "Protein", id=batch_query_ids,
            ) as handle:
                batch_post = Entrez.read(handle)

            # eLink Protein to Nuccore db and retrieve IDs of linked nucleotide records
            # {protein record ID: {nucleotide records IDs}}
            nucleotide_ids = get_linked_nucleotide_record_ids(batch_post)
            if nucleotide_ids is None:
                failed_gbk_batches.append(cleaned_batch)
                continue  # onto the next batch of protein accessions

            # retrieves the nucleotide records IDs for protein records for whcih only one
            # nuclotide ID was retrieved
            single_nucleotide_ids = set()

            # retrieve the protein_record_ids of protein records from which multiple nucletoide 
            # record IDs were retrieved
            protein_records_multi_nuc = set()

            for protein_record_id in nucleotide_ids:
                if len(nucleotide_ids[protein_record_id]) == 1:
                    single_nucleotide_ids.add(list(nucleotide_ids[protein_record_id])[0])
                else:
                    protein_records_multi_nuc.add(protein_record_id)

            # batch query to fetch nucletoide records for protein records
            # from which only a sinlge nucleotide ID was retrieved
            retrieved_proteins, newly_retrieved_proteins, succcess = extract_protein_accessions(
                single_nucleotide_ids,
                retrieved_proteins,
                gbk_accessions,
                args,
            )
            if succcess is False:
               failed_nuc_batches.append(single_nucleotide_ids)
                continue  # on to the next batch of protein accessions

            # add the nucleotide accessions to the nucleotide_accessions_dict
            # {kingdom: {genus: {species: {nucleotide record accessio: {protein_accessions},},},},}
            for protein_accession in newly_retrieved_proteins:
                nucleotide_accession = retrieved_proteins[protein_accession]

                try:
                    species = gbk_organism_dict[protein_accession]['species']
                    genus = gbk_organism_dict[protein_accession]['genus']
                except KeyError:
                    logger.error(
                        f"Retrieved protein accession from NCBI that is not in the local CAZyme db"
                    )
                    with open(no_accession_logger, 'a') as fh:
                        fh.write(
                            f"{protein_accession}\t"
                            "Retrieved protein accession from NCBI that is not in the local CAZyme db\t"
                            f"{genus} {species}\n"
                        )
                    continue

            bioproject_acc = None
            try:
                xrefs = record['GBSeq_xrefs']
                for xref_dict in xrefs:
                    try:
                        if xref_dict['GBXref_dbname'].strip() == 'BioProject':
                            bioproject_acc = xref_dict['GBXref_id']
                    except KeyError:
                        pass
            except KeyError:
                pass

            if bioproject_acc is not None:
                genomic_bioproject_ids = add_bioproject_id(
                    genomic_bioproject_ids,
                    bioproject_acc,
                    kingdom,
                    genus,
                    species,
                    protein_accession,
                )
                continue

            else:
                # attempt to retrieve assembly name
                assembly_name = None

                try:
                    for item in record['GBSeq_comment'].split(";"):
                        if item.strip().startswith("Assembly Name ::"):
                            assembly_name = item.strip()[len("Assembly Name :: "):]
                except KeyError as err:
                    pass

                if assembly_name is not None:
                    genomic_assembly_names = add_assembly_name(
                        genomic_assembly_names,
                        assembly_name,
                        kingdom,
                        genus,
                        species,
                        protein_accession,
                    )
                    continue

                else:
                    message = (
                        f"Could not retrieve BioProject ID or assembly name for {protein_accession}\t"
                        f"{kingdom}: {genus} {species}"
                    )
                    logger.warning(message)
                    with open(no_accession_logger, 'a') as fh:
                        fh.write(message)
                    continue

    return genomic_bioproject_ids, genomic_assembly_names


def get_linked_nucleotide_record_ids(batch_post, all_parsed_nucleotide_ids, args):
    """Use Entrez.elink to get the ID of the NCBI.Nucleotide db record linked to the NCBI.Protein
    db record.
    
    In the resulting dict, the protein record ID is used to group nucleotide records that are 
    linked to the same protein record, and from which it will need to be determined which is
    the longest record, so that only one (and the longest) Nucleotide record is parsed for this protein
    not all retrieved Nucleotide records.

    :param batch_post: Entrez dict containing WebEnv and query_key from batch posting protein accessions
    :param all_parsed_nucleotide_ids: set of IDs of previusly downloaded and parsed Nucleotide records
    :param args: cmd-line args parser

    Return dict {protein_accession: {nucleotide records IDs}}
    Return None if cannot retrieve data from NCBI
    """
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
        except Exception:
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


def extract_protein_accessions(single_nucleotide_ids, retrieved_proteins, gbk_accessions, args):
    """Retrieve and parse Nucleotide db records, for NCBI.Protein records from which only
    one NCBI.Nucleotide db record ID was retrieved.
    
    :param retrieved_proteins: dict, {protein_accession: nucleotide record accession}
    :param single_nucleotide_ids: list of nucloetide record IDs
    :param gbk_accessions: list of protein GenBank accessions from the local CAZyme database
    :param args: cmd-line args parser
    
    Return retrieved_proteins (dict) and bool
        True if successful Entrez connection, False is connection fails
    """
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
        except Exception:
            return retrieved_proteins, False
    
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

    return retrieved_proteins, True


def add_assembly_name(genomic_assembly_names, assembly_name, kingdom, genus, species, protein_accession):
    """Add assembly name to dict.
    
    :param:
    
    Return dict
    """
    try:
        genomic_assembly_names[kingdom]

        try:
            genomic_assembly_names[kingdom][genus]

            try:
                genomic_assembly_names[kingdom][genus][species]

                try:
                    genomic_assembly_names[kingdom][genus][species][assembly_name].add(protein_accession)
                
                except KeyError:
                    genomic_assembly_names[kingdom][genus][species][assembly_name] = {protein_accession}

            except KeyError:
                genomic_assembly_names[kingdom][genus][species] = {
                    assembly_name: {protein_accession},
                }
            
        except KeyError:
            genomic_assembly_names[kingdom][genus] = {
                species: {
                    assembly_name: {protein_accession},
                },
            }

    except KeyError:
        genomic_assembly_names[kingdom] = {
            genus: {
                species: {
                    assembly_name: {protein_accession},
                },
            },
        }
    
    return genomic_assembly_names


def add_bioproject_id(genomic_bioproject_ids, bioproject_acc, kingdom, genus, species, protein_accession):
    """Add assembly name to dict.
    
    :param:
    
    Return dict
    """
    try:
        genomic_bioproject_ids[kingdom]

        try:
            genomic_bioproject_ids[kingdom][genus]

            try:
                genomic_bioproject_ids[kingdom][genus][species]

                try:
                    genomic_bioproject_ids[kingdom][genus][species][bioproject_acc].add(protein_accession)
                
                except KeyError:
                    genomic_bioproject_ids[kingdom][genus][species][bioproject_acc] = {protein_accession}

            except KeyError:
                genomic_bioproject_ids[kingdom][genus][species] = {
                    bioproject_acc: {protein_accession},
                }
            
        except KeyError:
            genomic_bioproject_ids[kingdom][genus] = {
                species: {
                    bioproject_acc: {protein_accession},
                },
            }

    except KeyError:
        genomic_bioproject_ids[kingdom] = {
            genus: {
                species: {
                    bioproject_acc: {protein_accession},
                },
            },
        }
    
    return genomic_bioproject_ids


def get_genomic_accessions(genomic_assembly_names, no_accession_logger, args):
    """Retrieve genomic accessions for the genomic assemblies

    :param genomic_assembly_names: dict
        {kingdom: {genus: {species: {assembly_name: {protein_accessions},},},},}
    :param no_accession_logger: Path, path to log file to write out assembly names for which no
        genomic accession was retrieved
    :param args: cmd-line args parser
    
    Return dict,
    {kingdom: {genus: {species: {genomic_accession: {proteins: {protein_accessions}, count=int},},},},}
    """
    logger = logging.getLogger(__name__)

    genomic_accession_dict = {}
    # {kingdom: {genus: {species: {genomic_accession: {proteins: {protein_accessions}, count=int},},},},}

    for kingdom in tqdm(genomic_assembly_names, desc='Retrieving genomic accessions per kingdom'):
        genera = genomic_assembly_names[kingdom]
        for genus in genera:
            organisms = genera[genus]
            for species in organisms:
                # retrieve all genomic assembly names for the given species
                assembly_names = list(organisms[species].keys()) 

                # break up the list into a series of smaller lists that can be batched querried
                batch_queries = get_chunks_list(args.batch_size, assembly_names)

                for batch_query in tqdm(batch_queries, desc=f"Batch querying for {genus} {species}"):
                    batch_query_ids = ",".join(batch_query)

                # retrieve the records IDs for the assembly names
                with entrez_retry(
                    args.retries,
                    Entrez.esearch,
                    "Assembly",
                    id=batch_query_ids,
                ) as handle:
                    batch_post = Entrez.read(handle)
                
                with entrez_retry(
                    args.retries,
                    Entrez.efetch,
                    db="Assembly",
                    query_key=batch_post['QueryKey'],
                    WebEnv=batch_post['WebEnv'],
                    retmode="xml",
                ) as handle:
                    batch_fetch = Entrez.read(handle)

                genomic_accessions = {}

                for genome_record in tqdm(batch_fetch, desc="Retrieving assembly IDs"):
                    index = 0
                    accessions = set()

                    for index in range(len(genome_record['IdList'])):
                        with entrez_retry(
                            10, Entrez.efetch, db="Assembly", id=genome_record['IdList'][index], retmode="xml", rettype="docsum",
                        ) as handle:
                            result = Entrez.read(handle)
                        genomic_accession = result['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
                        assembly_name = result['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']

                        genomic_accessions[genomic_accession] = assembly_name

                    accessions = list(genomic_accessions.keys)
                    accessions.sort(reverse=True)
                    latest_accession = accessions[0]
                    latest_assembly_name = genomic_accessions[latest_accession]

                    # replace assemlby name for genomic accession
                    try:
                        protein_accessions = genomic_assembly_names[kingdom][genus][species][latest_assembly_name]

                    except KeyError:
                        logger.warning(f"Retrieved assembly name {latest_assembly_name}, but not retrieved previously")
                        with open(no_accession_logger, 'a') as fh:
                            fh.write(
                                f"{latest_assembly_name}\tRetrieved assembly name, but not retrieved previously\t"
                                f"{latest_accession}\t{genus} {species}\n"
                            )
                        continue

                    try:
                        genomic_accession_dict[kingdom]

                        try:
                            genomic_accession_dict[kingdom][genus]

                            try:
                                genomic_accession_dict[kingdom][genus][species]

                            except KeyError:
                                genomic_accession_dict[kingdom][genus][species] = {
                                    genomic_accession: {
                                        'proteins': protein_accessions,
                                        'count': len(protein_accessions),
                                    },
                                }

                        except KeyError:
                            genomic_accession_dict[kingdom][genus] = {
                                species: {
                                    genomic_accession: {
                                        'proteins': protein_accessions,
                                        'count': len(protein_accessions),
                                    },
                                },
                            }

                    except KeyError:
                        genomic_accession_dict[kingdom] = {
                            genus: {
                                species: {
                                    genomic_accession: {
                                        'proteins': protein_accessions,
                                        'count': len(protein_accessions),
                                    },
                                },
                            },
                        }

    return genomic_accession_dict


def write_out_genomic_accessions(genomic_accession_dict, time_stamp, args):
    """Compile and write output listing the genomic accessions.
    
    :param genomic_accession_dict:
    :param time_stamp: str, date and time script was invoked
    :param args: cmd-line args parser
    
    Return nothing
    """
    raw_json_path = args.output_dir / f"cazy_genomic_accessions_{time_stamp}.json"
    with open(raw_json_path, 'w') as fh:
        json.dump(genomic_accession_dict, fh)

    genomic_df_column_names = ['Kingdom','Genus','Species','Genomic_accession','#ofProteins']
    genomic_df = pd.DataFrame(columns=genomic_df_column_names)

    protein_df_column_names = ['Kingdom', 'Genus', 'Species', 'Genomic_accession', 'Protein_accession']
    protein_df = pd.DataFrame(columns=protein_df_column_names)

    for kingdom in genomic_accession_dict:
        genera = genomic_accession_dict[kingdom]

        for genus in genera:
            organisms = genera[genus]
            
            for organism in organisms:
                genomes = organisms[organism]
            
                for genomic_accession in genomes:
                    protein_accessions = organisms[genomic_accession]['Proteins']
                    number_of_proteins = organisms[genomic_accession]['Count']

                    g_new_row_data = [kingdom, genus, organism, genomic_accession, number_of_proteins]
                    g_new_row = pd.DataFrame([g_new_row_data], columns=genomic_df_column_names)
                    genomic_df = genomic_df.append(g_new_row)

                    for protein in protein_accessions:
                        p_new_row_data = [kingdom, genus, organism, genomic_accession, protein]
                        p_new_row = pd.DataFrame([p_new_row_data], columns=protein_df_column_names)
                        protein_df = protein_df.append(p_new_row)
    
    genomic_csv = args.output_dir / f"genomic_accessions_{time_stamp}.csv"
    genomic_df.to_csv(genomic_csv)

    protein_csv = args.output_dir / f"protein_genomic_accessions_{time_stamp}.csv"
    protein_df.to_csv(protein_csv)

    return


def get_ncbi_counts(args):
    """Retrieve the number of genomic assemblies per Kingdom from NCBI.
    
    :param args: cmd-line args parser
    
    Return dict {kingdom: count}
    """
    counts = {}
    
    for kingdom in KINGDOMS:
        with entrez_retry(
            args.retries,
            Entrez.esearch,
            db="Assembly",
            term=kingdom,
            retmode="xml",
        ) as record_handle:
            record = Entrez.read(record_handle, validate=False)
    
        number_of_assemblies = record['Count']
        counts[kingdom] = number_of_assemblies

    return counts


def write_out_genome_coverage(ncbi_genomes_totals, genomic_accession_dict, time_stamp, args):
    """Write out the genome coverage of NCBI GenBank database by the local CAZyme database
    
    :param ncbi_genomes_totals: dict {kingdom: number of NCBI GenBank genomes
    :param genomic_accession_dict: dict
        {kingdom: {genus: {species: {accession: {proteins: set(), counts: int}}}}
    :param time_stamp: str, date and time script was invoked
    :param args: cmd-line args parser
    
    Return nothing
    """
    column_names = ['Kingdom', 'NCBI_genomes', 'CAZy_genomes', 'Coverage_percent']
    coverage_df = pd.DataFrame(columns=column_names)
    graph_columns = ['Kingdom', 'NCBI', 'CAZy']
    graph_df = pd.DateFrame(columns=graph_columns)

    for kingdom in KINGDOMS:
        ncbi = ncbi_genomes_totals[kingdom]

        cazy = 0
        genera = genomic_accession_dict[kingdom]
        for genus in genera:
            organisms = genera[genus]
            for species in organisms:
                species_genome_accessions = len(list(organisms[species].keys()))
                cazy += species_genome_accessions

        coverage = (cazy / ncbi) * 100

        row_data = [kingdom, ncbi, cazy, coverage]
        new_row = pd.DataFrame([row_data], columns=column_names)
        coverage_df = coverage_df.append(new_row)

        row_data = [kingdom, int(ncbi), int(cazy)]
        new_row = pd.DataFrame([row_data], columns=graph_columns)
        graph_df = graph_df.append(new_row)

    output_path = args.output_dir / f"cazy_genbank_genome_coverage_{time_stamp}.csv"
    coverage_df.to_csv(output_path)

    fig, ax = plt.subplots()
    # plot CAZy bars
    ax.bar(
        graph_df['Kingdom'],
        graph_df['CAZy'],
        label='CAZy',
        color='orange',
    )
    # add NCBI bars (the higher bars)
    ax.bar(
        graph_df['Kingdom'],
        graph_df['NCBI'],
        bottom=graph_df['CAZy'],
        label='NCBI',
        color='dodgerblue',
    )
    ax.set_ylabel('Kingdom')
    ax.set_xlabel('Number of genomes in the database')
    ax.set_title('GenBank genomes included in CAZy')
    
    ax.legend()

    output_path = args.output_dir / f"gbk_cazy_genomes_plot_{time_stamp}.png"
    fig.savefig(output_path, bbox_inches='tight', dpi=360)

    return


if __name__ == "__main__":
    main()
