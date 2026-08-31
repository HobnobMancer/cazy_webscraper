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
"""Retrieve proteins with no genome data in the local database."""


import sqlite3

from src.sql.interface.get_data.get_proteins import (
    PROTEIN_FILTER_JOINS,
    build_protein_filters,
)

from tqdm import tqdm


def get_cached_proteins(protein_accs: list[str], db_connection: sqlite3.Connection) -> set[str]:
    """Get protein accessions that are already cached in temp tables

    Args:
        protein_accs: List of protein accessions to check
        db_connection: Database connection

    Returns:
        Set of protein accessions that are already cached
    """
    cursor = db_connection.cursor()
    placeholders = ','.join('?' * len(protein_accs))
    cursor.execute(f"""
        SELECT DISTINCT p.protein_accession 
        FROM Proteins p 
        JOIN TEMP_GENOME2PROTEIN tg2p ON p.protein_id = tg2p.protein_id
        WHERE p.protein_accession IN ({placeholders})
    """, protein_accs)

    return set(row[0] for row in cursor.fetchall())


def get_genome_table(genomes_of_interest, connection):
    """Load genome table into a dict

    :param genomes_of_interest: list of genome names
    :param connection: open sql db connection

    Return dict {genome name: db id}
    """
    if not genomes_of_interest:
        return {}

    # Create placeholders for the IN clause
    placeholders = ','.join('?' * len(genomes_of_interest))

    cursor = connection.execute(
        f"""SELECT genome_name, genome_id 
            FROM Genomes 
            WHERE assembly_name IN ({placeholders})""",
        genomes_of_interest
    )

    db_genome_dict = {}
    for row in cursor:
        genome_name, genome_id = row
        db_genome_dict[genome_name] = genome_id

    return db_genome_dict


def get_protein_genome_table_data(connection):
    """Parse the Genbanks_Genomes table into a set of tuples, one row per tuple.

    :param connection: open sql db connection

    Return set of tuples [(protein_id, genome_id), ...]
    """
    cursor = connection.execute("SELECT protein_id, genome_id FROM Genbanks_Genomes")

    prot_gnm_records = set()
    for row in cursor.fetchall():
        prot_gnm_records.add((row[0], row[1]))

    return prot_gnm_records


def get_genomes(gbk_dict, args, connection):
    """Retrieve genomic version accessions for proteins in gbk_dict

    I TIHNK THIS IS NEEDED FOR GTDB TAX RETRIEVAL

    :param gbk_dict: dict, gbk_ver_acc: local db ID
    :param args: CLI argument parser
    :param connection: open connection to a SQLite db engine

    Return dict {local db genome id: {
        'gbk_genome': str-version acc,
        'refseq_genome': str
    }}
    """
    genome_dict = {}

    for gbk in tqdm(gbk_dict, desc="Getting genomes for proteins of interest"):
        cursor = connection.execute(
            """SELECT Gn.gbk_version_accession, Gn.refseq_version_accession, Gn.genome_id, Gn.gtdb_tax_id 
               FROM Genomes AS Gn 
               INNER JOIN Genbanks_Genomes AS GG ON Gn.genome_id = GG.genome_id 
               INNER JOIN Genbanks AS Gb ON GG.genbank_id = Gb.genbank_id 
               WHERE Gb.genbank_accession = ?""",
            (gbk,)
        )
        
        results = cursor.fetchall()
        
        if results:
            for result in results:
                gbk_version_acc, refseq_version_acc, genome_id, gtdb_tax_id = result
                
                if gtdb_tax_id is not None:
                    if args.update_genome_lineage is False:
                        continue

                if genome_id not in genome_dict:
                    genome_dict[genome_id] = {}

                if gbk_version_acc is not None:
                    if 'gkb_genomes' not in genome_dict[genome_id]:
                        genome_dict[genome_id]['gkb_genomes'] = set()
                    genome_dict[genome_id]['gkb_genomes'].add(gbk_version_acc)

                if refseq_version_acc is not None:
                    if 'ref_genomes' not in genome_dict[genome_id]:
                        genome_dict[genome_id]['ref_genomes'] = set()
                    genome_dict[genome_id]['ref_genomes'].add(refseq_version_acc)

    selected_genomes = set()
    for genome_db_id in genome_dict:
        if 'gkb_genomes' in genome_dict[genome_db_id]:
            selected_genomes.update(genome_dict[genome_db_id]['gkb_genomes'])
        if 'ref_genomes' in genome_dict[genome_db_id]:
            selected_genomes.update(genome_dict[genome_db_id]['ref_genomes'])

    return genome_dict, selected_genomes


def get_selected_genomes(
    class_filters: set,
    family_filters: set,
    kingdom_filters: set,
    taxonomy_filter_dict: dict,
    ec_filters: set,
    db_path,
    protein_accessions: list = None,
    only_missing_gtdb: bool = False,
) -> dict:
    """Retrieve the genomes of proteins matching the user's filter criteria.

    A genome carries both a GenBank (GCA_) and a RefSeq (GCF_) accession, and GTDB indexes its
    lineages under whichever of the two it built the entry from, so both are returned and the
    caller matches on either.

    Returns {assembly_accession: genome_id}, with an entry for each accession a genome has.
    """
    query = "SELECT DISTINCT G.genome_id, G.gbk_version_accession, G.refseq_version_accession" \
        + PROTEIN_FILTER_JOINS + """
    JOIN Proteins_Genomes PG ON P.protein_id = PG.protein_id
    JOIN Genomes G ON PG.genome_id = G.genome_id
    """

    clause, params = build_protein_filters(
        class_filters, family_filters, kingdom_filters, taxonomy_filter_dict, ec_filters
    )
    query += clause

    if protein_accessions is not None:
        placeholders = ','.join('?' for _ in protein_accessions)
        query += f" AND P.protein_accession IN ({placeholders})"
        params.extend(protein_accessions)

    if only_missing_gtdb:
        query += " AND G.gtdb_tax_id IS NULL"

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(query, params)

    accession_to_genome_id = {}
    for genome_id, gbk_accession, refseq_accession in cursor:
        for accession in (gbk_accession, refseq_accession):
            if accession:
                accession_to_genome_id[accession] = genome_id

    conn.close()

    return accession_to_genome_id
