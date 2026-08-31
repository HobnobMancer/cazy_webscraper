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
"""Assemble the protein records requested by the extract_data subcommand.

Data is yielded a batch of proteins at a time rather than assembled in one dict, so the
memory used does not scale with the size of the database being exported.
"""


import sqlite3

from pathlib import Path

from saintBioutils.misc import get_chunks_list


# what each --include option adds to a protein record
INCLUDE_CHOICES = (
    "class",
    "family",
    "subfamily",
    "kingdom",
    "genus",
    "organism",
    "ec",
    "pdb",
    "uniprot_acc",
    "uniprot_name",
    "genbank_seq",
    "uniprot_seq",
)


def _fetch_many_to_many(conn, accessions: list[str], query: str) -> dict[str, set]:
    """Run a query returning (protein_accession, value) pairs and group by accession."""
    placeholders = ','.join('?' for _ in accessions)
    grouped = {}
    for accession, value in conn.execute(query.format(placeholders=placeholders), accessions):
        if value is None:
            continue
        grouped.setdefault(accession, set()).add(value)
    return grouped


def _assemble_batch(conn, accessions: list[str], include: set[str]) -> dict[str, dict]:
    """Build the record dicts for one batch of protein accessions."""
    placeholders = ','.join('?' for _ in accessions)
    records = {accession: {} for accession in accessions}

    # CAZy class/family/subfamily all come from the same join
    if {"class", "family", "subfamily"} & include:
        families = {}
        subfamilies = {}
        for accession, family, subfamily in conn.execute(
            f"""SELECT P.protein_accession, F.family, F.subfamily
                FROM Proteins P
                JOIN Proteins_CazyFamilies PF ON P.protein_id = PF.protein_id
                JOIN CazyFamilies F ON PF.family_id = F.family_id
                WHERE P.protein_accession IN ({placeholders})""",
            accessions,
        ):
            if family:
                families.setdefault(accession, set()).add(family)
            if subfamily:
                subfamilies.setdefault(accession, set()).add(subfamily)

        for accession, record in records.items():
            fams = families.get(accession, set())
            if "family" in include:
                record["family"] = fams
            if "subfamily" in include:
                record["subfamily"] = subfamilies.get(accession, set())
            if "class" in include:
                # the CAZy class is the leading alphabetic part of the family, e.g. GH13 -> GH
                record["class"] = {
                    ''.join(c for c in fam if c.isalpha()) for fam in fams
                }

    if {"kingdom", "genus", "organism"} & include:
        for accession, kingdom, genus, species in conn.execute(
            f"""SELECT P.protein_accession, K.kingdom, T.genus, T.species
                FROM Proteins P
                LEFT JOIN Taxs T ON P.taxonomy_id = T.taxonomy_id
                LEFT JOIN Kingdoms K ON T.kingdom_id = K.kingdom_id
                WHERE P.protein_accession IN ({placeholders})""",
            accessions,
        ):
            record = records[accession]
            if "kingdom" in include:
                record["kingdom"] = kingdom
            if "genus" in include:
                record["genus"] = genus
            if "organism" in include:
                record["organism"] = " ".join(p for p in (genus, species) if p) or None

    if "ec" in include:
        ecs = _fetch_many_to_many(conn, accessions, """
            SELECT P.protein_accession, E.ec_number
            FROM Proteins P
            JOIN Proteins_Ecs PE ON P.protein_id = PE.protein_id
            JOIN Ecs E ON PE.ec_id = E.ec_id
            WHERE P.protein_accession IN ({placeholders})""")
        for accession, record in records.items():
            record["ec"] = ecs.get(accession, set())

    if "pdb" in include:
        pdbs = _fetch_many_to_many(conn, accessions, """
            SELECT P.protein_accession, D.pdb_accession
            FROM Proteins P
            JOIN Proteins_Pdbs PP ON P.protein_id = PP.protein_id
            JOIN Pdbs D ON PP.pdb_id = D.pdb_id
            WHERE P.protein_accession IN ({placeholders})""")
        for accession, record in records.items():
            record["pdb"] = pdbs.get(accession, set())

    if {"uniprot_acc", "uniprot_name", "uniprot_seq"} & include:
        for accession, uniprot_acc, name, sequence in conn.execute(
            f"""SELECT P.protein_accession, U.uniprot_accession, U.name, U.sequence
                FROM Proteins P
                LEFT JOIN Uniprots U ON P.uniprot_id = U.uniprot_id
                WHERE P.protein_accession IN ({placeholders})""",
            accessions,
        ):
            record = records[accession]
            if "uniprot_acc" in include:
                record["uniprot_acc"] = uniprot_acc
            if "uniprot_name" in include:
                record["uniprot_name"] = name
            if "uniprot_seq" in include:
                record["uniprot_seq"] = sequence

    if "genbank_seq" in include:
        for accession, sequence in conn.execute(
            f"""SELECT protein_accession, sequence FROM Proteins
                WHERE protein_accession IN ({placeholders})""",
            accessions,
        ):
            records[accession]["genbank_seq"] = sequence

    return records


def iter_protein_records(
    accessions: list[str],
    db_path: Path,
    include: set[str],
    batch_size: int = 1000,
):
    """Yield {accession: record} dicts, one batch of proteins at a time."""
    conn = sqlite3.connect(db_path)
    try:
        for batch in get_chunks_list(accessions, batch_size):
            yield _assemble_batch(conn, batch, set(include))
    finally:
        conn.close()


def iter_sequences(
    accessions: list[str],
    db_path: Path,
    sources: set[str],
    batch_size: int = 1000,
):
    """Yield (accession, source, sequence) tuples for the requested sequence sources.

    'genbank' yields the NCBI sequence stored against the protein; 'uniprot' yields the
    sequence of the UniProt record linked to it, keyed by its UniProt accession.
    """
    conn = sqlite3.connect(db_path)
    try:
        for batch in get_chunks_list(accessions, batch_size):
            placeholders = ','.join('?' for _ in batch)

            if "genbank" in sources:
                for accession, sequence in conn.execute(
                    f"""SELECT protein_accession, sequence FROM Proteins
                        WHERE protein_accession IN ({placeholders})""",
                    batch,
                ):
                    if sequence:
                        yield accession, "GenBank", sequence

            if "uniprot" in sources:
                for uniprot_acc, sequence in conn.execute(
                    f"""SELECT U.uniprot_accession, U.sequence
                        FROM Proteins P
                        JOIN Uniprots U ON P.uniprot_id = U.uniprot_id
                        WHERE P.protein_accession IN ({placeholders})""",
                    batch,
                ):
                    if sequence:
                        yield uniprot_acc, "UniProt", sequence
    finally:
        conn.close()
