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
from cazy_webscraper.genomes import entrez
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

    nucleotide_accessions_dict = get_nucleotide_accessions(
        genbank_kingdom_dict,
        no_accession_logger,
        args,
    )

    output_path = cache_dir / f"nucleotide_accessions_{time_stamp}.json"
    with open(output_path, 'w') as fh:
        json.dump(nucleotide_accessions_dict, fh)
    
    genomic_accession_dict = get_genomic_accessions(nucleotide_accessions_dict, no_accession_logger, args)
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


def get_nucleotide_accessions(genbank_kingdom_dict, no_accession_logger, args):
    """Retrieve the NCBI Nucleotide db records ID's containing the GenBank protein accessions.
    
    :param genbank_kingdom_dict: dict of Genbank and Kingdom records from db
        {kingdom: {genus: {species: {protein_accessions}}}
    :param no_accession_logger: Path, path to log file to save protein accessions for which no 
        genomic accession was retrieved
    :param args: cmd-line args parser
    
    Return dict {kingdom: {genus: {species: {nucleotide_accession: {protein_accessions},},},},}
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

        # gbk accessions waiting for a nucleotide id to be retrieved
        remaining_accessions = list(gbk_organism_dict.keys())

        while len(remaining_accessions) != 0:
            starting_loop_length = len(remaining_accessions)
            logger.warning(f"{len(remaining_accessions)} Gbk accessions remaining to be parsed for {kingdom}")
            
            accessions_to_parse = remaining_accessions[:args.batch_size]

            # eLink Protein to Nuccore db and retrieve IDs of linked nucleotide records
            batch_query_ids = ",".join(accessions_to_parse)
            with entrez_retry(
                args.retries, Entrez.epost, "Protein", id=batch_query_ids,
            ) as handle:
                batch_post = Entrez.read(handle)

            # {protein record ID: {nucleotide records IDs}}
            nucleotide_ids = entrez.get_linked_nucleotide_record_ids(batch_post, args)

            if nucleotide_ids is None: 
                # issue with at least one accession in the batch
                # e.g. it is not longer stored in NCBI
                # pass individually to find/parse the bad accession(s)
                nucleotide_ids, no_nucleotides = entrez.link_nucleotide_ids_individually(
                    accessions_to_parse,
                    no_accession_logger,
                    args,
                )

                if len(no_nucleotides) != 0:
                    logger.warning(f"{len(no_nucleotides)} proteins found with no linked Nucleotide records")
                    for protein_accession in no_nucleotides:
                        # already logged in link_nucleotide_ids_individually()
                        try:
                            remaining_accessions.remove(protein_accession)
                        except ValueError:
                            pass
            
            if len(list(nucleotide_ids.values())) == 0:
                # no linked Nucleotide records retrieved for the current batch of protein accessions
                for protein_accession in accessions_to_parse:
                    logger.warning(
                        f"Could not reitreve linked Nucleotide record for {protein_accession}"
                    )
                    # do not try to retrieve record again
                    try:
                        remaining_accessions.remove(protein_accession)
                    except ValueError:
                        pass
                continue

            # retrieves the nucleotide records IDs for protein records for whcih only one
            # nuclotide ID was retrieved
            single_nucleotide_ids = set()

            # retrieve the protein_record_ids of protein records from which multiple nucletoide 
            # record IDs were retrieved
            protein_records_multi_nuc = set()

            for protein_record_id in tqdm(nucleotide_ids, desc="Idenitfying proteins with multiple linked nucleotide records"):
                if len(nucleotide_ids[protein_record_id]) == 1:
                    single_nucleotide_ids.add(list(nucleotide_ids[protein_record_id])[0])
                else:
                    logger.warning(
                        f"Found {len(nucleotide_ids[protein_record_id])} linked nucletoide records "
                        f"for protein record {protein_record_id}"
                        )
                    protein_records_multi_nuc.add(protein_record_id)

            # batch query to fetch nucletoide records for protein records
            # from which only a sinlge nucleotide ID was retrieved
            if len(single_nucleotide_ids) != 0:
                retrieved_proteins, newly_retrieved_proteins, succcess = entrez.extract_protein_accessions(
                    single_nucleotide_ids,
                    retrieved_proteins,
                    gbk_accessions,
                    args,
                )
                if succcess is False:
                    # issue with at least one accession in the batch
                    # e.g. it is not longer stored in NCBI
                    # pass individually to find/parse the bad accession(s)
                    retrieved_proteins, newly_retrieved_proteins = entrez.extract_protein_accessions_individually(
                        single_nucleotide_ids,
                        retrieved_proteins,
                        gbk_accessions,
                        args,
                    )

                # add the nucleotide accessions to the nucleotide_accessions_dict
                # {kingdom: {genus: {species: {nucleotide record accessio: {protein_accessions},},},},}
                nucleotide_accessions_dict = add_nucleotide_accessions(
                    nucleotide_accessions_dict,
                    gbk_organism_dict,
                    retrieved_proteins,
                    newly_retrieved_proteins,
                    kingdom,
                    no_accession_logger,
                )

            # for Protein records for which multiple Nucleotide record IDs were retrieved
            # Identify the longest Nucleotide record and retrieve protein accessions from it
            # The longest record is most likely to be the most complete record
            for protein_record_id in tqdm(protein_records_multi_nuc, "Parsing protein records with multiple linked Nucletide records"):
                nucleotide_record_ids = nucleotide_ids[protein_record_id]

                retrieved_proteins, newly_retrieved_proteins, success = entrez.parse_longest_record(
                    nucleotide_record_ids,
                    retrieved_proteins,
                    gbk_accessions,
                    args,
                )

                if success is False:
                    # issue with at least one accession in the batch
                    # e.g. it is not longer stored in NCBI
                    # pass individually to find/parse the bad accession(s)
                    retrieved_proteins, newly_retrieved_proteins, success = entrez.parse_longest_record_individually(
                        nucleotide_record_ids,
                        retrieved_proteins,
                        gbk_accessions,
                        args,
                    )

                nucleotide_accessions_dict = add_nucleotide_accessions(
                    nucleotide_accessions_dict,
                    gbk_organism_dict,
                    retrieved_proteins,
                    newly_retrieved_proteins,
                    kingdom,
                    no_accession_logger,
                )

            # remove protein accessions from remaining_accessions because the linked Nucleotide 
            # record ID has already been retrieved
            for protein_accession in retrieved_proteins:
                try:
                    remaining_accessions.remove(protein_accession)
                except ValueError:
                    pass

            if starting_loop_length == len(remaining_accessions) and len(remaining_accessions) != 0:
                # failing to retrieve data for protein accessions
                nucleotide_ids, no_nucleotides = entrez.link_nucleotide_ids_individually(
                    accessions_to_parse,
                    no_accession_logger, 
                    args,
                )

                if nucleotide_ids is None:
                    for protein_accession in remaining_accessions:
                        logger.warning(
                                f"Could not retrieve  Nucletoide record for {protein_accession}"
                            )
                        try:
                            remaining_accessions.remove(protein_accession)
                        except ValueError:
                            pass
                    continue
                
            if len(list(nucleotide_ids.values())) == 0:
                # no linked Nucleotide records retrieved for the current batch of protein accessions
                for protein_accession in accessions_to_parse:
                    logger.warning(
                        f"Could not reitreve linked Nucleotide record for {protein_accession}"
                    )
                    # do not try to retrieve record again
                    try:
                        remaining_accessions.remove(protein_accession)
                    except ValueError:
                        pass
        
            if len(no_nucleotides) != 0:
                for protein_accession in no_nucleotides:
                    # already logged in link_nucleotide_ids_individually()
                    try:
                        remaining_accessions.remove(protein_accession)
                    except ValueError:
                        pass
            
    return nucleotide_accessions_dict


def add_nucleotide_accessions(
    nucleotide_accessions_dict,
    gbk_organism_dict,
    retrieved_proteins,
    newly_retrieved_proteins,
    kingdom,
    no_accession_logger,
):
    """Add Nucleotide accession and protein acccession to the nucleotide_accessions_dict.
    
    :param nucleotide_accessions_dict: dict
        {kingdom: {genus: {species: {nucleotide record accessio: {protein_accessions},},},},}
    :param gbk_organism_dict: dict,
        {protein_accession: {species: str, genus: str}}
    :param retrieved_proteins: dict
        {protein record id: {nucleotide records ids}}
    :param newly_retrieved_proteins: list of protein's retrieved from Nucleotide records and 
        not yet added to nucleotide_accessions_dict
    :param kingdom: str, taxonomic kingdom of protein's source organism
    :param no_accession_logger: Path, path to log file to save protein accessions for which no 
            genomic accession was retrieved
    
    Return dict
    """
    logger = logging.getLogger(__name__)

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
    try:
        nucleotide_accessions_dict[kingdom]

        try:
            nucleotide_accessions_dict[kingdom][genus]

            try:
                nucleotide_accessions_dict[kingdom][genus][species]

                try:
                    nucleotide_accessions_dict[kingdom][genus][species][nucleotide_accession].add(protein_accession)
                
                except KeyError:
                    nucleotide_accessions_dict[kingdom][genus][species][nucleotide_accession] = {protein_accession}

            except KeyError:
                nucleotide_accessions_dict[kingdom][genus][species] = {
                    nucleotide_accession: {protein_accession},
                }
            
        except KeyError:
            nucleotide_accessions_dict[kingdom][genus] = {
                species: {
                    nucleotide_accession: {protein_accession},
                },
            }

    except KeyError:
        nucleotide_accessions_dict[kingdom] = {
            genus: {
                species: {
                    nucleotide_accession: {protein_accession},
                },
            },
        }
    
    return nucleotide_accessions_dict


def get_genomic_accessions(nucleotide_accessions_dict, no_accession_logger, args):
    """Retrieve genomic accessions for the genomic assemblies

    :param nucleotide_accessions_dict: dict
        {kingdom: {genus: {species: {nucleotide record accession: {protein_accessions},},},},}
    :param no_accession_logger: Path, path to log file to write out assembly names for which no
        genomic accession was retrieved
    :param args: cmd-line args parser
    
    Return dict,
    {kingdom: {genus: {species: {genomic_accession: {proteins: {protein_accessions}, count=int},},},},}
    """
    logger = logging.getLogger(__name__)

    genomic_accession_dict = {}
    # {kingdom: {genus: {species: {genomic_accession: {proteins: {protein_accessions}, count=int},},},},}

    for kingdom in tqdm(nucleotide_accessions_dict, desc='Retrieving genomic accessions per kingdom'):
        genera = nucleotide_accessions_dict[kingdom]
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
                        protein_accessions = nucleotide_accessions_dict[kingdom][genus][species][latest_assembly_name]

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
