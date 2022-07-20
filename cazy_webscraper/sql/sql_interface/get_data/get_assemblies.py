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
