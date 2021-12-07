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
    make_output_directory(cache_dir, args.force_cache, args.nodelete_cache)

    no_accession_logger = cache_dir / "logs"
    make_output_directory(cache_dir, args.force_cache, args.nodelete_cache)
    no_accession_logger = no_accession_logger / f"no_genomic_accession_retrieved_{time_stamp}.log"

    # load Genbank and Kingdom records from the db
    logger.warning("Retrieving Genbanks, Taxs and Kingdoms records from the local CAZyme db")
    genbank_kingdom_dict = get_table_dicts.get_gbk_kingdom_dict(connection)
    logger.warning("Retrieved Genbanks, Taxs and Kingdoms records from the local CAZyme db")

    genomic_assembly_names = get_assebmly_names(genbank_kingdom_dict, no_accession_logger, args)

    genomic_accession_dict = get_genomic_accessions(genomic_assembly_names, no_accession_logger, args)

    write_out_genomic_accessions(genomic_accession_dict, time_stamp, args)

    ncbi_genomes_totals = get_ncbi_counts(args)

    write_out_genome_coverage(ncbi_genomes_totals, genomic_accession_dict, time_stamp, args)

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    if args.verbose:
        logger.info(
            "Finished getting data from UniProt\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_scraper.VERSION_INFO}\n"
            f"Citation: {cazy_scraper.CITATION_INFO}"
        )
    else:
        print(
            "=====================cazy_webscraper=====================\n"
            "Finished getting data from UniProt\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {cazy_scraper.VERSION_INFO}\n"
            f"Citation: {cazy_scraper.CITATION_INFO}"
        )


def get_assebmly_names(genbank_kingdom_dict, no_accession_logger, args):
    """Retrieve assembly names of the source genomic accessions for the protein accessions in the local db.
    
    :param genbank_kingdom_dict: dict of Genbank and Kingdom records from db
        {kingdom: {genus: {species: {protein_accessions}}}
    :param no_accession_logger: Path, path to log file to save protein accessions for which no 
        genomic accession was retrieved
    :param args: cmd-line args parser
    
    Return dict {kingdom: {genus: {species: {assembly_name: {protein_accessions},},},},}
    """
    logger = logging.getLogger(__name__)

    genomic_assembly_names = {}
    # {kingdom: {genus: {species: {assembly_name: {protein_accessions},},},},}
    
    for kingdom in tqdm(genbank_kingdom_dict, desc="Retrieving genomic assembly names per kingdom"):
        genera = genbank_kingdom_dict[kingdom]
        for genus in genera:
            organisms = genera[genus]
            for species in organisms:
                # retrieve all Gbk protein accessions for the given organism
                gbk_accessions = list(organisms[species])

                # break up the list into a series of smaller lists that can be batched querried
                batch_queries = get_chunks_list(gbk_accessions, args.batch_size)

                for batch_query in tqdm(batch_queries, desc=f"Batch querying NCBI for {genus} {species}"):
                    batch_query_ids = ",".join(batch_query)
                    with entrez_retry(
                        args.retries, Entrez.epost, "Protein", id=batch_query_ids,
                    ) as handle:
                        batch_post = Entrez.read(handle)

                    # eFetch against the Protein database, retrieve in xml retmode
                    with entrez_retry(
                        args.retries,
                        Entrez.efetch,
                        db="Protein",
                        query_key=batch_post['QueryKey'],
                        WebEnv=batch_post['WebEnv'],
                        retmode="xml",
                    ) as handle:
                        batch_fetch = Entrez.read(handle)
                    
                    for record in tqdm(batch_fetch, desc="Retrieving assembly name from protein records"):
                        assembly_name = None
                        protein_accession = record['GBSeq_accession-version']

                        for item in record['GBSeq_comment'].split(";"):
                            if item.strip().startswith("Assembly Name ::"):
                                assembly_name = item.strip()[len("Assembly Name :: "):]

                        if assembly_name is None:
                            logger.warning(
                                f"Could not retrieve genome assembly name for {protein_accession}"
                            )
                            with open(no_accession_logger, 'a') as fh:
                                fh.write(
                                    f"{protein_accession}\tNo genome assembly name retrieved\t"
                                    f"Protein record comment: {record['GBSeq_comment']}\t"
                                    f"{genus} {species}\n"
                                )
                            continue
                        
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
