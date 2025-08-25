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
"""Retrieve proteins with no assembly data in the local database."""

import sqlite3
from tqdm import tqdm


def get_no_assembly_proteins(gbk_dict, connection):
    """Filter a gbk_dict to retain only those proteins with no assembly data in the local db

    :param gbk_dict: dict, {protein gbk acc: db id}
    :param connection: open sqlite db connection

    Return gbk_dict"""
    filtered_gbk_dict = {}

    for gbk_acc in tqdm(gbk_dict, desc="Filtering for proteins with no assembly data in the db"):
        cursor = connection.execute(
            """SELECT Genbanks.genbank_id, Genomes.genome_id 
               FROM Genbanks 
               INNER JOIN Genbanks_Genomes ON Genbanks.genbank_id = Genbanks_Genomes.genbank_id
               INNER JOIN Genomes ON Genbanks_Genomes.genome_id = Genomes.genome_id
               WHERE Genbanks.genbank_accession = ?""",
            (gbk_acc,)
        )
        
        results = cursor.fetchall()
        if results:
            filtered_gbk_dict[gbk_acc] = gbk_dict[gbk_acc]

    return filtered_gbk_dict


def get_records_to_update(gbk_dict, connection):
    """Filter a gbk_dict to retain only those proteins with no assembly data in the local db

    :param gbk_dict: dict, {protein gbk acc: db id}
    :param connection: open sqlite db connection

    Return gbk_dict"""
    update_gbk_dict = {}  # proteins to update the new genome data
    add_gbk_dict = {}  # proteins to add new genome data

    for gbk_acc in tqdm(gbk_dict, desc="Filtering for proteins with no assembly data in the db"):
        cursor = connection.execute(
            """SELECT Genbanks.genbank_id, Genomes.genome_id 
               FROM Genbanks 
               INNER JOIN Genbanks_Genomes ON Genbanks.genbank_id = Genbanks_Genomes.genbank_id
               INNER JOIN Genomes ON Genbanks_Genomes.genome_id = Genomes.genome_id
               WHERE Genbanks.genbank_accession = ?""",
            (gbk_acc,)
        )
        
        results = cursor.fetchall()
        if len(results) == 0:
            add_gbk_dict[gbk_acc] = gbk_dict[gbk_acc]
        else:
            update_gbk_dict[gbk_acc] = gbk_dict[gbk_acc]

    return update_gbk_dict, add_gbk_dict


def get_assembly_table(genomes_of_interest, connection):
    """Load assembly table into a dict

    :param genomes_of_interest: list of assembly names
    :param connection: open sql db connection

    Return dict {assembly name: db id}
    """
    if not genomes_of_interest:
        return {}
    
    # Create placeholders for the IN clause
    placeholders = ','.join('?' * len(genomes_of_interest))
    
    cursor = connection.execute(
        f"""SELECT assembly_name, genome_id 
            FROM Genomes 
            WHERE assembly_name IN ({placeholders})""",
        genomes_of_interest
    )
    
    db_genome_dict = {}
    for row in tqdm(cursor.fetchall(), desc="Retrieving genome records from the local db"):
        assembly_name, genome_id = row
        db_genome_dict[assembly_name] = genome_id

    return db_genome_dict


def get_gbk_genome_table_data(connection):
    """Parse the Genbanks_Genomes table into a set of tuples, one row per tuple.

    :param connection: open sql db connection

    Return set of tuples
    """
    cursor = connection.execute("SELECT genbank_id, genome_id FROM Genbanks_Genomes")
    
    prot_gnm_records = set()
    for row in cursor.fetchall():
        prot_gnm_records.add((row[0], row[1]))

    return prot_gnm_records


def get_genomes(gbk_dict, args, connection):
    """Retrieve genomic version accessions for proteins in gbk_dict

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
