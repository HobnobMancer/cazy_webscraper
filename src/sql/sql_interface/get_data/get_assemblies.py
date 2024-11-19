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


from sqlalchemy import text
from tqdm import tqdm

from cazy_webscraper.sql.sql_orm import (
    Genbank,
    Genome,
    Session,
    genbanks_genomes,
)


def get_no_assembly_proteins(gbk_dict, connection):
    """Filter a gbk_dict to retain only those proteins with no assembly data in the local db

    :param gbk_dict: dict, {protein gbk acc: db id}
    :param connection: open sqlite db connection

    Return gbk_dict"""
    filtered_gbk_dict = {}

    for gbk_acc in tqdm(gbk_dict, desc="Filtering for proteins with no assembly data in the db"):
        with Session(bind=connection) as session:
            query_result = session.query(Genbank, Genome).\
                join(Genome, (Genome.genbank_id == Genbank.genbank_id)).\
                    filter(Genbank.genbank_accession == gbk_acc).\
                        all()

        for result in query_result:
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
        with Session(bind=connection) as session:
            query_result = session.query(Genbank, Genome).\
                join(Genome, (Genome.genbank_id == Genbank.genbank_id)).\
                    filter(Genbank.genbank_accession == gbk_acc).\
                        all()

        if len(query_result) == 0:
            add_gbk_dict[gbk_acc] = gbk_dict[gbk_acc]
        else:
            update_gbk_dict[gbk_acc] = gbk_dict[gbk_acc]

    return update_gbk_dict, add_gbk_dict


def get_assembly_table(genomes_of_interest, connection):
    """Load assembly table into a dict

    :param genomes_of_interest: list of assmebly names
    :param connection: open sql db connection

    Return dict {assembly name: db id}
    """
    with Session(bind=connection) as session:
        genome_records = session.query(Genome).all()

    db_genome_dict = {}  # {assembly name: db id}

    for record in tqdm(genome_records, desc="Retrieving genome records from the local db"):
        assembly_name = record.assembly_name
        db_id = record.genome_id

        if assembly_name in genomes_of_interest:
            db_genome_dict[assembly_name] = db_id

    return db_genome_dict


def get_gbk_genome_table_data(connection):
    """Parse the Genbanks_Genomes table into a set of tuples, one row per tuple.

    :param connection: opwn sql db connection

    Return set of tuples
    """
    with Session(bind=connection) as session:
        table_objs = session.query(genbanks_genomes).all()

    prot_gnm_records = set()

    for record in table_objs:
        prot_gnm_records.add((record.genbank_id, record.genome_id))

    return prot_gnm_records


def get_genomes(gbk_dict, args, connection):
    """Retrieve genomic version accessions for proteins in gbk_dict

    :param gbk_dict: dict, gbk_ver_acc: local db ID
    :param args: CLI argument parser
    :param connection: open connectoin to a SQLite db engine

    Return dict {local db genome id: {
        'gbk_genome': str-version acc,
        'refseq_genome': str
    }}
    """
    genome_dict = {}

    for gbk in tqdm(gbk_dict, desc="Getting genomes for proteins of interest"):
        with connection.begin():
            cmd = text(
                "SELECT Gn.gbk_version_accession, Gn.refseq_version_accession, Gn.genome_id, Gn.gtdb_tax_id "
                "FROM Genomes AS Gn "
                "INNER JOIN Genbanks_Genomes AS GG ON Gn.genome_id = GG.genome_id "
                "INNER JOIN Genbanks AS Gb ON GG.genbank_id = Gb.genbank_id "
                f"WHERE Gb.genbank_accession = '{gbk}'"
            )
            query_result = connection.execute(cmd).fetchall()

        if len(query_result) != 0:
            for result in query_result:
                if result[3] is not None:
                    if args.update_genome_lineage is False:
                        continue

                genome_db_id = result[2]
                try:
                    genome_dict[genome_db_id]
                except KeyError:
                    genome_dict[genome_db_id] = {}

                if result[0] is not None:
                    try:
                        genome_dict[genome_db_id]['gkb_genomes'].add(result[0])
                    except KeyError:
                        genome_dict[genome_db_id]['gkb_genomes'] = {result[0]}

                if result[1] is not None:
                    try:
                        genome_dict[genome_db_id]['ref_genomes'].add(result[1])
                    except KeyError:
                        genome_dict[genome_db_id]['ref_genomes'] = {result[1]}

    selected_genomes = set()
    for genome_db_id in genome_dict:
        try:
            for genome in genome_dict[genome_db_id]['gkb_genomes']:
                selected_genomes.add(genome)
        except KeyError:
            pass
        try:
            for genome in genome_dict[genome_db_id]['ref_genomes']:
                selected_genomes.add(genome)
        except KeyError:
            pass

    return genome_dict, selected_genomes
