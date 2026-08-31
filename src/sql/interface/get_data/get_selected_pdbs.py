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
"""Retrieve PDB accessions matching user criteria from the local CAZyme db"""


import sqlite3

from pathlib import Path

from src.sql.interface.get_data.get_proteins import (
    PROTEIN_FILTER_JOINS,
    build_protein_filters,
)


def get_pdb_accessions(
    class_filters: set[str],
    family_filters: set[str],
    kingdom_filters: set[str],
    taxonomy_filter_dict: dict[str, set[str]],
    ec_filters: set[str],
    db_path: Path,
    protein_accessions: list[str] = None,
) -> list[str]:
    """Retrieve PDB accessions for proteins matching the user's criteria.

    Uses the same filter clause as the protein queries, joined through to the PDB tables, so
    the whole selection is done in one query rather than by loading the Proteins, Pdbs and
    Proteins_Pdbs tables into memory as dicts (the v2 approach).

    If protein_accessions is given, the selection is further restricted to the PDB accessions
    of those proteins (used by the --genbank_accessions/--uniprot_accessions options).

    Returns a sorted list of unique parent PDB accessions (any ``[chain]`` annotation removed).
    """
    query = "SELECT DISTINCT Pdbs.pdb_accession" + PROTEIN_FILTER_JOINS + """
    JOIN Proteins_Pdbs ON P.protein_id = Proteins_Pdbs.protein_id
    JOIN Pdbs ON Proteins_Pdbs.pdb_id = Pdbs.pdb_id
    """

    clause, params = build_protein_filters(
        class_filters, family_filters, kingdom_filters, taxonomy_filter_dict, ec_filters
    )
    query += clause

    if protein_accessions is not None:
        placeholders = ','.join('?' for _ in protein_accessions)
        query += f" AND P.protein_accession IN ({placeholders})"
        params.extend(protein_accessions)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(query, params)
    # CAZy/UniProt list some structures with a chain suffix, e.g. 1ABC[A]; the structure file
    # is per-entry, so the chain annotation is dropped and duplicates collapsed
    accessions = {row[0].split("[")[0] for row in cursor if row[0]}
    conn.close()

    return sorted(accessions)
