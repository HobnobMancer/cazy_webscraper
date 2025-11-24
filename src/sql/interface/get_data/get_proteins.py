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


import sqlite3

from pathlib import Path


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
    SELECT Proteins.protein_accession
    FROM Proteins
    LEFT JOIN Proteins_CazyFamilies ON Proteins.protein_id = Proteins_CazyFamilies.protein_id
    LEFT JOIN CazyFamilies ON Proteins_CazyFamilies.family_id = CazyFamilies.family_id
    LEFT JOIN Taxs ON Proteins.taxonomy_id = Taxs.taxonomy_id
    LEFT JOIN Kingdoms ON Taxs.kingdom_id = Kingdoms.kingdom_id
    LEFT JOIN Proteins_Ecs ON Proteins.protein_id = Proteins_Ecs.protein_id
    LEFT JOIN ECs ON Proteins_Ecs.ec_id = ECs.ec_id
    """

    if additional_join:
        query += f" {additional_join}"

    query += " WHERE Proteins.source = 'ncbi'"

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

    if additional_filter:
        query += f" AND ({additional_filter})"

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(query, params)
    seq_acc_to_retrieve = [row[0] for row in cursor]
    conn.close()

    return seq_acc_to_retrieve
