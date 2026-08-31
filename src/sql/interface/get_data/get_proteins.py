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
# SOFTWARE.
"""Funcs for retrieving GenBank accessions of interest, matching user filter criteria"""


import logging
import sqlite3

from pathlib import Path

logger = logging.getLogger(__name__)


# The Proteins table is aliased as P, and the class/family/taxonomy/EC joins below are
# shared by every "which proteins match the user's filters" query, so they live in one place.
PROTEIN_FILTER_JOINS = """
    FROM Proteins P
    LEFT JOIN Proteins_CazyFamilies ON P.protein_id = Proteins_CazyFamilies.protein_id
    LEFT JOIN CazyFamilies ON Proteins_CazyFamilies.family_id = CazyFamilies.family_id
    LEFT JOIN Taxs ON P.taxonomy_id = Taxs.taxonomy_id
    LEFT JOIN Kingdoms ON Taxs.kingdom_id = Kingdoms.kingdom_id
    LEFT JOIN Proteins_Ecs ON P.protein_id = Proteins_Ecs.protein_id
    LEFT JOIN ECs ON Proteins_Ecs.ec_id = ECs.ec_id
"""


def build_protein_filters(
    class_filters: set[str],
    family_filters: set[str],
    kingdom_filters: set[str],
    taxonomy_filter_dict: dict[str, set[str]],
    ec_filters: set[str],
) -> tuple[str, list]:
    """Build the shared WHERE clause restricting proteins to the user's filter criteria.

    Returns the SQL fragment (starting with the WHERE) and its bound parameters.
    """
    query = " WHERE P.source = 'ncbi'"
    params = []

    if class_filters:
        query += " AND (" + " OR ".join("CazyFamilies.family LIKE ?" for _ in class_filters) + ")"
        params.extend(f'{cls}%' for cls in class_filters)

    if family_filters:
        query += " AND (" + " OR ".join("CazyFamilies.family = ?" for _ in family_filters) + ")"
        params.extend(family_filters)

    if kingdom_filters:
        query += " AND (" + " OR ".join("Kingdoms.kingdom = ?" for _ in kingdom_filters) + ")"
        params.extend(kingdom_filters)

    if taxonomy_filter_dict["genera"]:
        query += " AND (" + " OR ".join("Taxs.genus = ?" for _ in taxonomy_filter_dict["genera"]) + ")"
        params.extend(taxonomy_filter_dict["genera"])

    if taxonomy_filter_dict["species"]:
        for species in taxonomy_filter_dict["species"]:
            genus, *species_parts = species.split()
            query += " AND (Taxs.genus = ? AND Taxs.species LIKE ?)"
            params.append(genus)
            params.append(f"{' '.join(species_parts)}%")

    if taxonomy_filter_dict["strains"]:
        for strain in taxonomy_filter_dict["strains"]:
            genus, *strain_parts = strain.split()
            query += " AND (Taxs.genus = ? AND Taxs.species = ?)"
            params.append(genus)
            params.append(' '.join(strain_parts))

    if ec_filters:
        query += " AND (" + " OR ".join("ECs.ec_number = ?" for _ in ec_filters) + ")"
        params.extend(ec_filters)

    return query, params


def get_ncbi_prot_accessions(
    class_filters: set[str],
    family_filters: set[str],
    kingdom_filters: set[str],
    taxonomy_filter_dict: dict[str, set[str]],
    ec_filters: set[str],
    db_path: Path,
    additional_join: str = None,
    additional_filter: str = None
) -> list[str]:
    query = """
    SELECT P.protein_accession
    FROM Proteins P
    LEFT JOIN Proteins_CazyFamilies ON P.protein_id = Proteins_CazyFamilies.protein_id
    LEFT JOIN CazyFamilies ON Proteins_CazyFamilies.family_id = CazyFamilies.family_id
    LEFT JOIN Taxs ON P.taxonomy_id = Taxs.taxonomy_id
    LEFT JOIN Kingdoms ON Taxs.kingdom_id = Kingdoms.kingdom_id
    LEFT JOIN Proteins_Ecs ON P.protein_id = Proteins_Ecs.protein_id
    LEFT JOIN ECs ON Proteins_Ecs.ec_id = ECs.ec_id
    """

    if additional_join:
        query += f" {additional_join}"

    clause, params = build_protein_filters(
        class_filters, family_filters, kingdom_filters, taxonomy_filter_dict, ec_filters
    )
    query += clause

    if additional_filter:
        query += f" AND ({additional_filter})"

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(query, params)
    seq_acc_to_retrieve = [row[0] for row in cursor]
    conn.close()

    return seq_acc_to_retrieve


def get_uniprot_accessions(
    class_filters: set[str],
    family_filters: set[str],
    kingdom_filters: set[str],
    taxonomy_filter_dict: dict[str, set[str]],
    ec_filters: set[str],
    db_path: Path,
    protein_accessions: list[str] = None,
    additional_filter: str = None,
) -> list[tuple[int, str]]:
    """Retrieve (protein_id, uniprot_accession) pairs for proteins matching the user's criteria.

    Only proteins that already have UniProt data in the local db are returned, since Pfam
    retrieval from InterPro is keyed by UniProt accession, not GenBank accession - run
    get_uniprot_data first to populate these.

    If protein_accessions is given, the selection is further restricted to those proteins
    (used by the --genbank_accessions/--uniprot_accessions options).
    """
    query = "SELECT DISTINCT P.protein_id, Uniprots.uniprot_accession" + PROTEIN_FILTER_JOINS + """
    JOIN Uniprots ON P.uniprot_id = Uniprots.uniprot_id
    """

    clause, params = build_protein_filters(
        class_filters, family_filters, kingdom_filters, taxonomy_filter_dict, ec_filters
    )
    query += clause

    if protein_accessions is not None:
        placeholders = ','.join('?' for _ in protein_accessions)
        query += f" AND P.protein_accession IN ({placeholders})"
        params.extend(protein_accessions)

    if additional_filter:
        query += f" AND ({additional_filter})"

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(query, params)
    uniprot_pairs = [(row[0], row[1]) for row in cursor]
    conn.close()

    return uniprot_pairs


def get_ncbi_acc_for_uniprot_acc(uniprot_accs: list[str], db: Path) -> set[str]:
    ncbi_acc = set()
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    placeholders = ','.join('?' for _ in uniprot_accs)
    query = f"""SELECT protein_accession
        FROM Proteins
        LEFT JOIN Uniprots ON Proteins.uniprot_id = Uniprots.uniprot_id
        WHERE Uniprots.uniprot_accession IN ({placeholders})"""
    cur.execute(query, tuple(uniprot_accs))
    for row in cur:
        ncbi_acc.add(row[0])

    if not ncbi_acc:
        logger.warning("Did not retrieve any NCBI accessions for the %s UniProt accessions", len(uniprot_accs))

    return ncbi_acc


def get_ncbi_acc_to_id(conn: sqlite3.Connection, ncbi_accs: set[str]) -> dict[str, int]:
    cur = conn.cursor()
    cur.execute("""SELECT protein_accession, protein_id FROM Proteins""")
    acc2id = {}
    for row in cur:
        if row[0] in ncbi_accs:
            acc2id[row[0]] = row[1]
    cur.close()
    return acc2id


def get_protein_table_dict(connection: sqlite3.Connection) -> dict:
    """Compile a dict of the data in the Proteins table
    Return dict {genbank_accession: 'taxa_id': int, 'protein_id': int}
    """
    prot_cur = connection.cursor()
    prot_cur.execute("""SELECT * FROM Proteins""")
    db_protein_dict = {}  # {genbank_accession: 'taxa_id': str, 'id': int}
    for row in prot_cur:
        # [0] protein_id, [1] protein_accession
        # [2] sequence, [3] sequence_update
        # [4] taxonomy_id, [5] ncbi_tax_id
        # [6] uniprot_id, [7] source
        db_protein_dict[f"{row[1]}"] = {
            'taxa_id': row[4],
            'protein_id': row[0],
            'sequence': row[2]
        }
    prot_cur.close()
    return db_protein_dict
