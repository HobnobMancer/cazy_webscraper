#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2025
# (c) University of Strathclyde 2025
# (c) James Hutton Institute 2025
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# ehobbs@ebi.ac.uk
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

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from src.scripts import extract_data


def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    parser = subps.add_parser(
        "extract_data",
        description="Extract data from a local CAZyme database to CSV, JSON, FASTA or a BLAST database",
        help="Write data held in the local CAZyme database out to files",
        formatter_class=ArgumentDefaultsHelpFormatter
    )

    # Add positional/required arguments
    parser.add_argument(
        "database",
        type=Path,
        help="Path to local CAZy database",
    )

    # Add optional arguments to parser

    filters_group = parser.add_argument_group("Filtering arguments")
    data_group = parser.add_argument_group("Output arguments")
    operational_group = parser.add_argument_group("Operational arguments")
    utilities_group = parser.add_argument_group("Utility arguments")

    filters_group.add_argument(
        "--classes",
        type=str,
        default=None,
        help="Classes from which all families are to be scraped. Separate classes by ','"
    )
    filters_group.add_argument(
        "--families",
        type=str,
        default=None,
        help="Families to scrape. Separate families by commas 'GH1,GH2'",
    )
    filters_group.add_argument(
        "--kingdoms",
        type=str,
        default=None,
        help=(
            "Kingdoms to scrape. Separate by a single comma.\n"
            "Options= archaea, bacteria, eukaryota, viruses, unclassified (not case sensitive)"
        ),
    )
    filters_group.add_argument(
        "--genera",
        type=str,
        default=None,
        help="Genera to restrict the scrape to"
    )
    filters_group.add_argument(
        "--species",
        type=str,
        default=None,
        help="Species (written as Genus Species) to restrict the scrape to"
    )
    filters_group.add_argument(
        "--strains",
        type=str,
        default=None,
        help=(
            "Specific strains of organisms to restrict the scrape to "
            "(written as Genus Species Strain)"
        ),
    )
    filters_group.add_argument(
        "--ec_filter",
        type=str,
        default=None,
        help="Limit retrieval to proteins annotated with the provided EC numbers. Separate EC numbers with single commas"
    )

    data_group.add_argument(
        "--file_types",
        nargs="+",
        choices=["csv", "json", "fasta", "fasta_dir", "blastdb"],
        default=["csv"],
        help=(
            "Output file type(s) to write. 'fasta' writes one multi-record file, 'fasta_dir' "
            "writes one file per protein, and 'blastdb' builds a BLAST protein database "
            "(requires makeblastdb from BLAST+ on your PATH)"
        ),
    )
    data_group.add_argument(
        "--include",
        nargs="+",
        choices=[
            "class", "family", "subfamily", "kingdom", "genus", "organism",
            "ec", "pdb", "uniprot_acc", "uniprot_name", "genbank_seq", "uniprot_seq",
        ],
        default=None,
        help=(
            "Additional data to include as columns in the csv/json output. "
            "Has no effect on the sequence output types"
        ),
    )
    data_group.add_argument(
        "--source",
        nargs="+",
        choices=["genbank", "uniprot"],
        default=["genbank"],
        help=(
            "Which protein sequences to write to the fasta/fasta_dir/blastdb outputs. "
            "Has no effect on the csv/json output types"
        ),
    )
    data_group.add_argument(
        "-o",
        "--output_dir",
        type=Path,
        default=None,
        help="Directory to write output files to. Default: the current working directory",
    )
    data_group.add_argument(
        "--prefix",
        type=str,
        default=None,
        help="String to prefix every output file name with",
    )
    data_group.add_argument(
        "--overwrite",
        dest="overwrite",
        action="store_true",
        default=False,
        help="Overwrite existing output files",
    )

    operational_group.add_argument(
        "-c",
        "--config",
        type=Path,
        metavar="Config file",
        default=None,
        help="Path to configuration file. Default: None, exports the entire database",
    )
    operational_group.add_argument(
        "--genbank_accessions",
        type=Path,
        default=None,
        help="Path to a text file containing a list of GenBank accessions to extract data for",
    )
    operational_group.add_argument(
        "--uniprot_accessions",
        type=Path,
        default=None,
        help="Path to a text file containing a list of UniProt accessions to extract data for",
    )

    utilities_group.add_argument(
        "--batch_size",
        type=int,
        default=1000,
        help="Number of proteins to assemble and write out at a time",
    )
    utilities_group.add_argument(
        "--cache_dir",
        type=Path,
        default=None,
        help="Path to cache directory",
    )
    utilities_group.add_argument(
        "--cazy_synonyms",
        type=Path,
        default=None,
        help="Path to JSON file containing CAZy class synoymn names",
    )
    utilities_group.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force writing in existing cache directory",
    )
    utilities_group.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help="Do not delete content in existing output directory",
    )
    utilities_group.add_argument(
        "--nodelete_cache",
        dest="nodelete_cache",
        action="store_true",
        default=False,
        help="Do not delete content in existing cache dir",
    )
    utilities_group.add_argument(
        "-r",
        "--retries",
        type=int,
        default=10,
        help="Number of times to retry a failed request to RCSB PDB",
    )
    utilities_group.add_argument(
        "-t",
        "--timeout",
        type=int,
        default=45,
        help="Time in seconds before a connection to RCSB PDB times out",
    )

    parser.set_defaults(func=extract_data.main)
