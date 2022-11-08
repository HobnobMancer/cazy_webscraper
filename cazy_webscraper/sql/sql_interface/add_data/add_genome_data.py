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
"""Add genomic assembly information to the local CAZyme database."""


import logging

from sqlalchemy import text
from tqdm import tqdm

from cazy_webscraper.sql.sql_interface import insert_data
from cazy_webscraper.sql.sql_interface.get_data.get_assemblies import (
    get_records_to_update,
    get_assembly_table,
    get_gbk_genome_table_data,
)


def add_assembly_data(assembly_prot_dict, ncbi_genome_dict, gbk_dict, connection, args):
    """Coordinate adding genomic data and associting proteins with genomes in the local database

    :param assembly_prot_dict: dict, {assembly_name: {protein_accessions}}
    :param ncbi_genome_dict: dict, listing assembly meta data -- data retrieved from ncbi
    :param gbk_dict: dict, {gbk_accession: local db id}
    :param connection: open sqlite db connection
    :param args: cmd-line args parser

    Return nothing
    """
    logger = logging.getLogger(__name__)

    db_genome_table_dict = get_assembly_table(ncbi_genome_dict, connection)

    genomes_to_add = []
    genomes_to_update = []

    for assembly_name in ncbi_genome_dict:
        try:
            db_genome_table_dict[assembly_name]
            genomes_to_update.append(assembly_name)

        except KeyError:  # assembly name not in the local db, so need to add new record
            genomes_to_add.append(assembly_name)

    genomes_of_interst = genomes_to_add + genomes_to_update  # list of assembly names

    if args.update and len(genomes_to_update) > 0:
        update_genomic_data(genomes_to_update, db_genome_table_dict, ncbi_genome_dict, connection)
        # load updated table
        db_genome_table_dict = get_assembly_table(genomes_of_interst, connection)

    if len(genomes_to_add) > 0:
        add_genomic_data(ncbi_genome_dict, genomes_to_add, connection)
        # load updated table
        db_genome_table_dict = get_assembly_table(genomes_of_interst, connection)

    # add genbank protein ids and genome ids to Genbanks_Genomes table to link proteins and genomes
    add_prot_genome_relationships(assembly_prot_dict, gbk_dict, db_genome_table_dict, connection)

    return


def add_genomic_data(ncbi_genome_dict, genomes_of_interst, connection):
    """Add new genomic assembly data to the local CAZyme database

    :param ncbi_genome_dict: dict, listing assembly meta data
    :param genomes_of_interst: list, assembly names of new assemblies to add to the db
    :param connection: open sqlite db connection

    Return nothing
    """
    logger = logging.getLogger(__name__)

    records_to_add = set()

    for assembly_name in tqdm(ncbi_genome_dict, desc="Compiling data to add to db"):
        if assembly_name not in genomes_of_interst:
            continue
        gbk_ver_acc = ncbi_genome_dict[assembly_name]['gbk_acc']
        gbk_ncbi_id = ncbi_genome_dict[assembly_name]['gbk_uid']
        ref_ver_acc = ncbi_genome_dict[assembly_name]['refseq_acc']
        ref_ncbi_id = ncbi_genome_dict[assembly_name]['refseq_uid']

        records_to_add.add((assembly_name, gbk_ver_acc, gbk_ncbi_id, ref_ver_acc, ref_ncbi_id))

    if len(records_to_add) != 0:
        insert_data(
            connection,
            'Genomes',
            ['assembly_name', 'gbk_version_accession', 'gbk_ncbi_id', 'refseq_version_accession', 'refseq_ncbi_id'],
            list(records_to_add),
        )

    return


def update_genomic_data(genomes_of_interest, genome_table_dict, ncbi_genome_dict, connection, unit_test=False):
    """Update existing genomic assembly data in the local CAZyme database

    :param genomes_of_interest: list, assembly names of records to update in the local db
    :param genome_table_dict: dict, assembly_name: db_dib
    :param ncbi_genome_dict: dict, listing assembly meta data
    :param connection: open sqlite db connection

    Return nothing
    """
    logger = logging.getLogger(__name__)

    for assembly_name in tqdm(ncbi_genome_dict, desc="Updating data in the db"):
        if assembly_name not in genomes_of_interest:
            continue

        db_id = genome_table_dict[assembly_name]

        gbk_ver_acc = ncbi_genome_dict[assembly_name]['gbk_acc']
        gbk_ncbi_id = ncbi_genome_dict[assembly_name]['gbk_uid']
        ref_ver_acc = ncbi_genome_dict[assembly_name]['refseq_acc']
        ref_ncbi_id = ncbi_genome_dict[assembly_name]['refseq_uid']

        with connection.begin():
            connection.execute(
                text(
                    "UPDATE Genomes "
                    f"SET gbk_version_accession = {gbk_ver_acc} AND "
                    f"gbk_ncbi_id = {gbk_ncbi_id} AND "
                    f"refseq_version_accession = {ref_ver_acc} AND "
                    f"refseq_ncbi_id = {ref_ncbi_id} AND "
                    f"WHERE genome_id = '{db_id}'"
                )
            )
            if unit_test:
                connection.rollback()

    return


def add_prot_genome_relationships(assembly_prot_dict, gbk_dict, db_genome_table_dict, connection):
    """Add protein (Gbk) and genome links to the Genbanks_Genomes table

    :param assembly_prot_dict: dict, {assembly_name: {protein_accessions}}
    :param db_genome_table_dict: Genomes table loaded into dict
    :param gbk_dict: dict, {gbk_accession: local db id}
    :param connection: open sqlite db connection

    Return nothing
    """
    # load Genbanks_Genomes table --> set of tuples, each tuple is one row (gbk_id, genome_id)
    gbk_genome_table_data = get_gbk_genome_table_data(connection)  
    relationships_to_add = set()

    for assembly_name in tqdm(assembly_prot_dict, desc="Linking proteins and genomes"):
        assembly_proteins = assembly_prot_dict[assembly_name]
        for protein in assembly_proteins:
            protein_id = gbk_dict[protein]
            assembly_id = db_genome_table_dict[assembly_name]

            relationship = (protein_id, assembly_id)

            if relationship not in gbk_genome_table_data:
                relationships_to_add.add(relationship)

    if len(relationships_to_add) > 0:
        insert_data(
            connection,
            'Genbanks_Genomes',
            ['genbank_id', 'genome_id'],
            list(relationships_to_add),
        )
