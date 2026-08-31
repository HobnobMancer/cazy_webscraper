#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
#
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
"""Coordinate extracting data from a local CAZyme database out to files."""


import logging

from argparse import Namespace
from pathlib import Path

from saintBioutils.utilities.file_io import make_output_directory

from src.export import check_existing_outputs, compile_output_path
from src.export.sequences import build_blast_db, write_sequences
from src.export.tabular import write_csv, write_json
from src.sql.interface.connect import connect_existing_db
from src.sql.interface.get_data.get_extract_data import (
    iter_protein_records,
    iter_sequences,
)
from src.sql.interface.get_data.get_proteins import (
    get_ncbi_acc_for_uniprot_acc,
    get_ncbi_prot_accessions,
)
from src.utilities.parse_configuration import get_expansion_configuration
from src.utilities.parse_configuration.accession_files import get_acc_from_file


logger = logging.getLogger(__name__)


# file types that are built from protein sequences rather than tabular records
SEQUENCE_FILE_TYPES = {"fasta", "fasta_dir", "blastdb"}


def main(args: Namespace, time_stamp: str, start_time):
    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)

    cache_dir = args.cache_dir if args.cache_dir else (cache_dir / "extract_data")
    make_output_directory(cache_dir, args.force, args.nodelete_cache)

    # NOTE: deliberately not make_output_directory() here. That helper deletes the whole
    # directory when --force is used, which is acceptable for our own cache dir but not for a
    # directory the user chose to write results into. Clobbering is guarded per file by
    # check_existing_outputs()/--overwrite instead, which is the right granularity.
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)

    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = get_expansion_configuration(args)

    file_types = set(args.file_types)
    include = set(args.include) if args.include else set()

    # work out every path this run will write, so existing files are caught before any work
    output_stem = compile_output_path(args)
    output_paths = {
        "csv": Path(f"{output_stem}.csv"),
        "json": Path(f"{output_stem}.json"),
        "fasta": Path(f"{output_stem}.fasta"),
        "blastdb": Path(f"{output_stem}_blastdb.fasta"),
    }
    fasta_dir = Path(f"{output_stem}_fasta") if "fasta_dir" in file_types else None

    to_check = [output_paths[ft] for ft in file_types if ft in output_paths]
    if not check_existing_outputs(to_check, args):
        return "extract_data"

    # select the proteins to export
    if args.genbank_accessions or args.uniprot_accessions:
        selected = set()
        if args.genbank_accessions:
            selected.update(get_acc_from_file(args.genbank_accessions, args.database))
        if args.uniprot_accessions:
            uniprot_accessions = get_acc_from_file(
                args.uniprot_accessions, args.database, table="UniProt"
            )
            selected.update(
                get_ncbi_acc_for_uniprot_acc(uniprot_accessions, args.database)
            )
        accessions = sorted(selected)
    else:
        accessions = get_ncbi_prot_accessions(
            class_filters,
            family_filters,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            args.database,
        )

    if not accessions:
        logger.warning("No proteins in the local database matched the criteria provided")
        return "extract_data"

    logger.warning("Extracting data for %s proteins", len(accessions))

    # tabular outputs
    if "csv" in file_types:
        rows = write_csv(
            iter_protein_records(accessions, args.database, include, args.batch_size),
            output_paths["csv"], include,
        )
        logger.warning("Wrote %s rows to %s", rows, output_paths["csv"])

    if "json" in file_types:
        records = write_json(
            iter_protein_records(accessions, args.database, include, args.batch_size),
            output_paths["json"], include,
        )
        logger.warning("Wrote %s records to %s", records, output_paths["json"])

    # sequence outputs
    if file_types & SEQUENCE_FILE_TYPES:
        if fasta_dir:
            fasta_dir.mkdir(parents=True, exist_ok=True)

        written = write_sequences(
            iter_sequences(accessions, args.database, set(args.source), args.batch_size),
            args,
            fasta_path=output_paths["fasta"] if "fasta" in file_types else None,
            fasta_dir=fasta_dir,
            blastdb_fasta=output_paths["blastdb"] if "blastdb" in file_types else None,
        )
        logger.warning("Wrote %s sequences", written)

        if not written:
            logger.warning(
                "No sequences were found for the selected proteins. Retrieve them first with "
                "'cazy_webscraper get_ncbi_seqs' (or 'get_uniprot_data --sequence')"
            )
        elif "blastdb" in file_types:
            build_blast_db(output_paths["blastdb"], output_stem.name)

    return "extract_data"
