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
from cazy_webscraper.sql.sql_interface.get_data.get_assemblies import get_records_to_update, get_assembly_table


def add_assembly_data(assembly_dict, genome_dict, gbk_dict, connection, args):
    """Coordinate adding genomic data and associting proteins with genomes in the local database
    
    :param assembly_dict: dict, {assembly_name: {protein_accessions}}
    :param genome_dict: dict, listing assembly meta data
    :param gbk_dict: dict, {gbk_accession: local db id}
    :param connection: open sqlite db connection
    :param args: cmd-line args parser
    
    Return nothing
    """
    logger = logging.getLogger(__name__)

    db_genome_table_dict = get_assembly_table(genome_dict, connection)

    genomes_to_add = {}
    new_genome_table = {}
    genomes_to_update = {}
    updated_genome_table = {}

    for assembly_name in genome_dict:
        try:
            db_genome_table_dict[assembly_name]
            genomes_to_update[assembly_name] = genome_dict[assembly_name]

        except KeyError:
            genomes_to_add[assembly_name] = genome_dict[assembly_name]

    if args.update:
        update_genomic_data(genomes_to_update, connection)
        updated_genome_table = get_assembly_table(genomes_to_update, connection)

    add_genomic_data(genomes_to_add, connection)
    new_genome_table = get_assembly_table(genomes_to_add, connection)

    # combine dicts
    new_genome_table.update(updated_genome_table)

    db_genome_table_dict = get_assembly_table(new_genome_table, connection)

    # identify protein records to update
    gbk_genome_dict = {}  # assembly db id: protein db id
    
    for assembly_name in tqdm(assembly_dict, desc="Identifying protein records to update"):
        assembly_proteins = assembly_dict[assembly_name]

        for protein in assembly_proteins:
            protein_id = gbk_dict[protein]
            assembly_id = db_genome_table_dict[assembly_name]
            gbk_genome_dict[protein_id] = assembly_id

    # update protein records
    update_protein_genomic_data(gbk_genome_dict, connection)

    return


def add_genomic_data(genome_dict, connection):
    """Add new genomic assembly data to the local CAZyme database
    
    :param genome_dict: dict, listing assembly meta data
    :param connection: open sqlite db connection

    Return nothing
    """
    logger = logging.getLogger(__name__)

    records_to_add = set()

    for assembly_name in tqdm(genome_dict, desc="Compiling data to add to db"):
        gbk_ver_acc = genome_dict[assembly_name]['gbk_acc']
        gbk_ncbi_id = genome_dict[assembly_name]['gbk_uid']
        ref_ver_acc = genome_dict[assembly_name]['refseq_acc']
        ref_ncbi_id = genome_dict[assembly_name]['refseq_uid']
        
        records_to_add.add((assembly_name, gbk_ver_acc, gbk_ncbi_id, ref_ver_acc, ref_ncbi_id))

    if len(records_to_add) != 0:
        insert_data(
            connection,
            'Genomes',
            ['assembly_name', 'gbk_version_accession', 'gbk_ncbi_id', 'refseq_version_accession', 'refseq_ncbi_id'],
            list(records_to_add),
        )
    
    return


def update_genomic_data(update_gbk_dict, assembly_dict, genome_dict, connection):
    """Update existing genomic assembly data in the local CAZyme database

    :param update_gbk_dict: dict, {gbk_accession: local db id}
    :param assembly_dict: dict, {assembly_name: {protein_accessions}}
    :param genome_dict: dict, listing assembly meta data
    :param connection: open sqlite db connection
    
    Return nothing
    """
    logger = logging.getLogger(__name__)

    for assembly_name in tqdm(assembly_dict, desc="Updating data in the db"):
        assembly_proteins = assembly_dict[assembly_name]
        gbk_ver_acc = genome_dict[assembly_name]['gbk_acc']
        gbk_ncbi_id = genome_dict[assembly_name]['gbk_uid']
        ref_ver_acc = genome_dict[assembly_name]['refseq_acc']
        ref_ncbi_id = genome_dict[assembly_name]['refseq_uid']

        for protein_acc in assembly_proteins:
            try:
                db_id = update_gbk_dict[protein_acc]
            except KeyError:
                continue  # protein not to be updated

            with connection.begin():
                connection.execute(
                    text(
                        "UPDATE Genomes "
                        f"SET assembly_name = {assembly_name} AND "
                        f"gbk_version_accession = {gbk_ver_acc} AND "
                        f"gbk_ncbi_id = {gbk_ncbi_id} AND "
                        f"refseq_version_accession = {ref_ver_acc} AND "
                        f"refseq_ncbi_id = {ref_ncbi_id} AND "
                        f"WHERE genbank_id = '{db_id}'"
                    )
                )
    
    return


def update_protein_genomic_data(gbk_genome_dict, connection):
    """Update the genomic assembly data in the protein (gbk) table
    
    :param gbk_genome_dict: dict, db protein id: assembl db id
    :param conneciton: open sql db connection
    
    Return nothing
    """
    for protein_id in tqdm(gbk_genome_dict, desc="Updating db protein records"):
        assembly_id = gbk_genome_dict[protein_id]

        with connection.begin():
            connection.execute(
                text(
                    "UPDATE Genbanks "
                    f"SET genome_id = {assembly_id} "
                    f"WHERE genbank_id = '{protein_id}'"
                )
            )

    return
