#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
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
#
# Bio.PDB reference:
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
# implemented in Python. Bioinformatics 19: 2308–2310
"""Extract protein sequences from local CAZyme database, write to FASTA files and/or BLAST db"""


import logging
import sys

import pandas as pd

from datetime import datetime
from typing import List, Optional

from tqdm import tqdm
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger

from cazy_webscraper import closing_message, connect_existing_db
from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks, get_table_dicts
from cazy_webscraper.sql.sql_interface.get_data.get_records import (
    get_user_genbank_sequences,
    get_user_uniprot_sequences
)
from cazy_webscraper.utilities import parse_configuration
from cazy_webscraper.utilities.parsers.extract_seq_parser import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)
    # parse cmd-line arguments
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        args = build_parser(argv).parse_args()

    if logger is None:
        logger = logging.getLogger(__package__)
        config_logger(args)

    validate_user_options(args)

    # make output directories
    if args.fasta_file:
        target_dir = args.fasta_file.parent
        make_output_directory(target_dir, args.force, args.nodelete)
    if args.fasta_dir:
        make_output_directory(args.fasta_dir, args.force, args.nodelete)

    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "sequence_extraction"
        make_output_directory(cache_dir, args.force, args.nodelete_cache)

    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = parse_configuration.get_expansion_configuration(args)

    gbk_table_dict = get_table_dicts.get_gbk_table_dict(connection)
    # {genbank_accession: 'taxa_id': str, 'gbk_id': int}

    # check what additional tabled needed to be loaded
    if any( ((args.genbank_accessions is not None), (args.uniprot_accessions is not None), ('genbank' in args.source)) ):
        logger.info("Loading the GenBank table")
        gbk_seq_dict = get_table_dicts.get_gbk_table_seq_dict(connection)
        # {genbank_accession: 'sequence': str, 'seq_date': str}

    if any( ((args.uniprot_accessions is not None), ('uniprot' in args.source)) ):
        logger.info("Loading the UniProt table")
        uniprot_table_dict = get_table_dicts.get_uniprot_table_dict(connection)
        # {acc: {name: str, gbk_id: int, seq: str, seq_date:str } }

    # build dick {gbk_acc: db_id} matching the users specified criteria
    # either via a list in a file or parameters provided via config file and/or command line

    gbk_dict = {}  # {gbk_acc: gbk_id}

    if args.genbank_accessions is not None:
        logger.warning(f"Extracting protein sequences for GenBank accessions listed in {args.genbank_accessions}")
        gbk_dict.update(get_user_genbank_sequences(gbk_table_dict, args))
    
    if args.uniprot_accessions is not None:
        logger.warning(f"Extracting protein sequences for UniProt accessions listed in {args.uniprot_accessions}")
        gbk_dict.update(get_user_uniprot_sequences(gbk_table_dict, uniprot_table_dict, args))

    if len(list(gbk_dict.keys())) == 0:
        gbk_dict = get_selected_gbks.get_genbank_accessions(
            class_filters,
            family_filters,
            taxonomy_filter_dict,
            kingdom_filters,
            ec_filters,
            connection,
        )

    # extract protein sequences from the database

    extracted_sequences = {}  # {accession: {'db': str, 'seq': str}}

    if 'genbank' in args.source:
        extracted_sequences.update(get_genbank_sequences(gbk_seq_dict, gbk_dict))

    if 'uniprot' in args.source:
        extracted_sequences.update(get_uniprot_sequences(uniprot_table_dict, gbk_dict))

    protein_records = []
    
    for protein_accession in tqdm(extracted_sequences, "Compiling SeqRecords"):
        try:
            new_record = SeqRecord(
                Seq(extracted_sequences[protein_accession]['seq']),
                id=protein_accession,
                description=extracted_sequences[protein_accession]['db']
            )
            protein_records.append(new_record)
        except TypeError:
            if extracted_sequences[protein_accession]['seq'] is not None:
                logger.warning(
                    f"Seq for {protein_accession} retrieved as type {type(extracted_sequences[protein_accession]['seq'])}\n"
                    "Not adding to FASTA file"
                )
            pass  # passed when sequence is None

    # write out the sequences to the specified outputs

    logger.warning(f"Extracted {len(protein_records)}")

    write_output(protein_records, cache_dir, args)

    closing_message("extract_sequences", start_time, args)
    

def validate_user_options(args):
    """Check the user has provided suitable operational options
    
    :param args: cmd-line args parser
    
    Return nothing
    """
    logger = logging.getLogger(__name__)

    if 'genbank' in args.source:
        logger.info("Extract GenBank protein sequences")
    
    if 'uniprot' in args.source:
        logger.info("Extracting UniProt protein sequences")

    if args.blastdb is None and args.fasta_dir is None and args.fasta_file is None:
        logger.error(
            "No output option provided.\n"
            "No path to create a BLAST db, output dir or output file provided\n"
            "Call at least one of --blastdb, --fasta_dir, --fasta_file\n"
            "Terminating program."
        )
        sys.exit(1)

    if args.blastdb:
        logger.info(f"Enabled building a BLAST db at:\n{args.blastdb}")
    
    if args.fasta_dir:
        logger.info(
            f"Enabled writing a unqiue sequence to a separate FASTA file in:\n{args.fasta_dir}"
        )
    
    if args.fasta_file:
        logger.info(
            f"Enabled writing all extracted sequences to a single FASTA file at:{args.fasta_file}"
        )

    return


def get_genbank_sequences(gbk_seq_dict, gbk_dict):
    """Extract protein seqeunces from the database

    :param gbk_seq_dict: dict {gbk_acc: {'sequence': str}}
    :param gbk_dict: dict of selected GenBank record {acc: id}
    
    Return dict {gbk_acc: {'db': 'genbank', 'seq': str}}
    """
    extracted_sequences = {}  # {accession: {'db': str, 'seq': str}}

    for gbk_accession in tqdm(gbk_dict, desc="Getting GenBank sequences"):
        extracted_sequences[gbk_accession] = {
            'db': 'GenBank',
            'seq': gbk_seq_dict[gbk_accession]['sequence'],
        }

    return extracted_sequences


def get_uniprot_sequences(uniprot_table_dict, gbk_dict):
    """Extract protein seqeunces from the database
    
    :param uniprot_table_dict: {acc: {name: str, gbk_id: int, seq: str, seq_date:str } }
    :param gbk_dict: dict of selected GenBank record {acc: id}
    
    Return dict {gbk_acc: {'db': 'genbank', 'seq': str}}
    """
    # get the db ids of selected gbks
    selected_genbanks = [gbk_dict[gbk_acc] for gbk_acc in gbk_dict]

    extracted_sequences = {}  # {accession: {'db': str, 'seq': str}}

    for uniprot_accession in tqdm(uniprot_table_dict, desc="Getting UniProt sequencse from db"):
        gbk_id = uniprot_table_dict[uniprot_accession]['genbank_id']
        if gbk_id in selected_genbanks:
            extracted_sequences[uniprot_accession] = {
                'db': 'UniProt',
                'seq': uniprot_table_dict[uniprot_accession]['seq'],
            }

    return extracted_sequences


def write_output(protein_records, cache_dir, args):
    """Write out the extract protein sequences to the specified output(s)
    
    :param protein_records: list of SeqRecords
    :param cache_dir: Path, cache directory
    :param args: cmd-line args parser
    
    Return nothing.
    """
    logger = logging.getLogger(__name__)

    if args.fasta_file:
        SeqIO.write(protein_records, args.fasta_file, "fasta")

    if args.fasta_dir:
        for record in protein_records:
            accession = record.id
            target_path = args.fasta_dir / f'{accession}.fasta'

            SeqIO.write([record], target_path, "fasta")
    
    if args.blastdb:
        # check if building a dir
        if str(args.blastdb).find("/") != -1:
            # build output dir
            make_output_directory(args.blastdb.parent, args.force, args.nodelete)

        fasta_name = f"{args.blastdb}_blastdb.fasta"
        SeqIO.write(protein_records, fasta_name, "fasta")

        cmd_makedb = NcbimakeblastdbCommandline(
            cmd='makeblastdb',
            dbtype='prot',
            input_file=fasta_name,
            title=args.blastdb,
        )
        stdout, stderr = cmd_makedb()

        # check the command was successfully exectured
        if len(stderr) != 0:
            logger.warning()
            print(f"Could not build non-CAZyme db.\nstdout={stdout}\nstderr={stderr}")

    cache_path = cache_dir / 'extracted_sequences.txt'
    with open(cache_path, 'a') as fh:
        for record in protein_records:
            fh.write(f"{record.id}\n")

    return


if __name__ == "__main__":
    main()
