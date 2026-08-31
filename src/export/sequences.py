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
# Bio.PDB / Biopython reference:
# Cock, P. J. A, Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A. et al. (2009)
# 'Biopython: freely available Python tools for computational molecular biology and
# bioinformatics', Bioinformatics, 25(11), pp. 1422-3.
"""Write protein sequences from a local CAZyme database out to FASTA files and BLAST databases."""


import logging
import subprocess

from argparse import Namespace
from pathlib import Path

from tqdm import tqdm


logger = logging.getLogger(__name__)


def _fasta_record(accession: str, source: str, sequence: str) -> str:
    """Format one sequence as a FASTA record, wrapped at 60 characters."""
    wrapped = "\n".join(sequence[i:i + 60] for i in range(0, len(sequence), 60))
    return f">{accession} {source}\n{wrapped}\n"


def write_sequences(
    sequence_iter,
    args: Namespace,
    fasta_path: Path = None,
    fasta_dir: Path = None,
    blastdb_fasta: Path = None,
) -> int:
    """Write sequences to the requested outputs as they are yielded.

    Records are streamed straight to the open file handles rather than collected into a list
    first, so exporting a whole database does not hold every sequence in memory.

    Returns the number of sequences written.
    """
    written = 0

    handles = []
    fasta_handle = open(fasta_path, "w") if fasta_path else None
    if fasta_handle:
        handles.append(fasta_handle)

    # the BLAST database is built from a FASTA file, so it is written the same way
    blastdb_handle = open(blastdb_fasta, "w") if blastdb_fasta else None
    if blastdb_handle:
        handles.append(blastdb_handle)

    try:
        for accession, source, sequence in tqdm(sequence_iter, desc="Writing sequences"):
            record = _fasta_record(accession, source, sequence)

            if fasta_handle:
                fasta_handle.write(record)

            if blastdb_handle:
                blastdb_handle.write(record)

            if fasta_dir:
                with open(fasta_dir / f"{accession}.fasta", "w") as fh:
                    fh.write(record)

            written += 1
    finally:
        for handle in handles:
            handle.close()

    return written


def build_blast_db(fasta_path: Path, db_title: str) -> bool:
    """Build a protein BLAST database from a FASTA file.

    Returns True on success. makeblastdb is invoked directly rather than through Biopython's
    command line wrappers: Bio.Blast.Applications (which the v2 code used, via
    NcbimakeblastdbCommandline) is deprecated in biopython 1.83, still present in 1.85, and
    gone by 1.88 - so under this project's biopython>=1.85 requirement it cannot be relied on.
    """
    command = [
        "makeblastdb",
        "-dbtype", "prot",
        "-in", str(fasta_path),
        "-title", str(db_title),
        "-out", str(fasta_path.with_suffix("")),
    ]

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        logger.error(
            "Could not find 'makeblastdb'. Install BLAST+ and make sure makeblastdb is on "
            "your PATH to build a BLAST database. The FASTA file at %s was still written.",
            fasta_path,
        )
        return False

    if result.returncode != 0:
        logger.error(
            "makeblastdb failed (exit code %s):\n%s\n%s",
            result.returncode, result.stdout.strip(), result.stderr.strip(),
        )
        return False

    logger.warning("Built BLAST database from %s", fasta_path)
    return True
