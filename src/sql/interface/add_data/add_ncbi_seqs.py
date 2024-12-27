#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
# (c) James Hutton Institute 2024
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
# SOFTWARE.
"""Add NCBI seqs to a local SQLite database"""


import logging
import sqlite3

from pathlib import Path

from Bio.Seq import Seq

from src.sql.interface.get_data.get_table_dicts import get_protein_table_dict

logger = logging.getLogger(__name__)


def update_ncbi_seqs(seq_dict: dict[str, Seq], db: Path, update: bool):
    """Add seqs from NCBI to a local CAZyme db.

    Update existing seqs in the db if enabled"""
    conn = sqlite3.connect(db)
    protein_table = get_protein_table_dict(conn)  # {protein_acc: dict['sequence': seq]}
    records_to_update = set()  # (str(seq), protein acc)

    for protein_acc, seq in seq_dict.items():
        if protein_acc not in protein_table:
            logger.warning("Retrieved seq for %s but this protein acc is not in the local db", protein_acc)
            continue
        if not seq_dict[protein_acc]:
            records_to_update.add( (str(seq), protein_acc) )
        elif seq_dict[protein_acc] and update:
            records_to_update.add( (str(seq), protein_acc) )

    if records_to_update:
        logger.warning("Updating %s seq records in the local CAZyme db", len(records_to_update))
        for record in records_to_update:
            conn.execute(
                """UPDATE Proteins SET sequence = ? WHERE protein_accession = ?""",
                record
            )

    conn.commit()
    conn.close()
