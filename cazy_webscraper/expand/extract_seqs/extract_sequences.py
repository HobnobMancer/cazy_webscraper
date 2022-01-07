#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
from Bio.SeqRecord import SeqRecord


from cazy_webscraper import cazy_webscraper
from cazy_webscraper.expand import get_chunks_gen
from cazy_webscraper.sql.sql_interface import get_selected_gbks, get_table_dicts
from cazy_webscraper.utilities import config_logger, file_io, parse_configuration
from cazy_webscraper.utilities.parsers import pdb_strctre_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)
    # parse cmd-line arguments
    if argv is None:
        parser = pdb_strctre_parser.build_parser()
        args = parser.parse_args()
    else:
        args = pdb_strctre_parser.build_parser(argv).parse_args()

    if logger is None:
        logger = logging.getLogger(__package__)
        config_logger(args)

    validate_user_options(args)

    # make output directories
    if args.fasta_file:
        target_dir = args.fasta_file.parent
        file_io.make_output_directory(target_dir, args.force, args.nodelete)
    if args.fasta_dir:
        file_io.make_output_directory(args.fasta_dir, args.force, args.nodelete)


    connection, logger_name, cache_dir = cazy_webscraper.connect_existing_db(args, time_stamp)

    if args.cache_dir is not None:  # use user defined cache dir
        cache_dir = args.cache_dir
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)
    else:
        cache_dir = cache_dir / "sequence_extraction"
        file_io.make_output_directory(cache_dir, args.force, args.nodelete_cache)

    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = parse_configuration.get_expansion_configuration(args)

    gbk_dict = get_selected_gbks.get_genbank_accessions(
        class_filters,
        family_filters,
        taxonomy_filter_dict,
        kingdom_filters,
        ec_filters,
        connection,
    )
    genbank_accessions = list(gbk_dict.keys())

    extracted_sequences = {}  # {accession: {'db': str, 'seq': str}}
    if args.genbank:
        extracted_sequences.update(get_genbank_sequences(gbk_dict, connection))
    if args.uniprot:
        extracted_sequences.update(get_uniprot_sequences(gbk_dict, connection))

    protein_records = []
    
    for protein_accession in tqdm(extracted_sequences, "Compiling SeqRecords"):
        new_record = SeqRecord(
            extracted_sequences[protein_accession]['seq'],
            id=protein_accession,
            description=extracted_sequences[protein_accession]['db']
        )
        protein_records.append(new_record)

    write_output(protein_records, args)
    

def validate_user_options(args):
    """Check the user has provided suitable operational options
    
    :param args: cmd-line args parser
    
    Return nothing
    """
    logger = logging.getLogger(__name__)

    if args.genbank is False and args.uniprot is False:
        logger.error(
            "No external sequence sources provided\n"
            "Call at least one of --genbank and --uniprot\n"
            "Terminating program."
        )
        sys.exit(1)

    if args.genbank:
        logger.info("Extract GenBank protein sequences")
    
    if args.uniprot:
        logger.info("Extracting UniProt protein sequences")

    if args.blastdb is False and args.fasta_dir is False and args.fasta_file is False:
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


def get_genbank_sequences(gbk_dict, connection):
    """Extract protein seqeunces from the database
    
    :param gbk_dict: dict of selected GenBank record {acc: id}
    :param connection: open sqlalchemy connection to an SQLite db
    
    Return dict {gbk_acc: {'db': 'genbank', 'seq': str}}
    """
    gbk_seq_dict = get_table_dicts.get_gbk_table_seq_dict(connection)

    extracted_sequences = {}  # {accession: {'db': str, 'seq': str}}

    for gbk_accession in gbk_dict:
        extracted_sequences[gbk_accession] = {
            'db': 'GenBank',
            'seq': gbk_seq_dict[gbk_accession]['sequence'],
        }

    return extracted_sequences


def get_uniprot_sequences(gbk_dict, connection):
    """Extract protein seqeunces from the database
    
    :param gbk_dict: dict of selected GenBank record {acc: id}
    :param connection: open sqlalchemy connection to an SQLite db
    
    Return dict {gbk_acc: {'db': 'genbank', 'seq': str}}
    """
    uniprot_table_dict = get_table_dicts.get_uniprot_table_dict(connection)
    # dict = {acc: {name: str, gbk_id: int, seq: str, seq_date:str } }

    gbk_uniprot_dict = {}  # {gbk_id: {uniprot_acc: str, seq: str}}
    
    for uniprot_accession in uniprot_table_dict:
        gbk_id = uniprot_table_dict[uniprot_accession]['gbk_id']
        seq = uniprot_table_dict[uniprot_accession]['seq']
        gbk_uniprot_dict[gbk_id] = {'uniprot_accession': uniprot_accession, 'seq': seq}

    extracted_sequences = {}  # {accession: {'db': str, 'seq': str}}

    for gbk_accession in gbk_dict:
        gbk_id = gbk_dict[gbk_accession]['gbk_id']

        extracted_sequences[gbk_accession] = {
            'db': 'UniProt',
            'seq': gbk_uniprot_dict[gbk_accession][gbk_id]['seq'],
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
        fasta_name = args.blastdb / 'blastdb.fasta'
        SeqIO.write(protein_records, target_path, "fasta")

        cmd_makedb = NcbimakeblastdbCommandline(
            cmd='makeblastdb',
            dbtype='prot',
            input_file=fasta_name,
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
