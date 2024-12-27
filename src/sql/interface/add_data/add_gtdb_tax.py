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
"""Add GTDB Taxonomy data to a local SQLite database"""


import logging

from sqlalchemy import text
from tqdm import tqdm

from cazy_webscraper.sql.sql_interface import insert_data
from cazy_webscraper.sql.sql_interface.get_data.get_table_dicts import (
    get_gtdb_table_dict,
    get_genome_table,
)
from cazy_webscraper.sql.sql_interface.get_data.get_assemblies import get_assembly_table


def add_gtdb_taxs(gtdb_lineages, connection):
    """Add GTDB lineages to GtdbTaxs table

    :param gtdb_lineags: dict
    {genome_version_accession: {'lineage': raw str lineage from gtdb. 'release': release-str}}
    :param connection: open connection to an SQlite db engine

    Return nothing
    """
    # loading existing table data
    existing_gtdb_table = get_gtdb_table_dict(connection)

    lineages_to_add = set()

    for genome in tqdm(gtdb_lineages, desc='Adding GTDB lineages to db'):
        lineage = [" ".join(_.split("__")[1:]) for _ in gtdb_lineages[genome]['lineage'].split(";")]
        if tuple(lineage) not in list(existing_gtdb_table.values()):
            lineage.append(gtdb_lineages[genome]['release'])
            lineages_to_add.add(tuple(lineage))

    if len(lineages_to_add) != 0:
        insert_data(
            connection,
            'GtdbTaxs',
            [
                'kingdom',
                'phylum',
                'tax_class',
                'tax_order',
                'family',
                'genus',
                'species',
                'release',
            ],
            list(lineages_to_add),
        )

    return


def add_genome_gtdb_relations(gtdb_lineages, args, connection):
    """Add genome - gtdb tax relationships
    
    :param gtdb_lineags: dict
    {genome_version_accession: {'lineage': raw str lineage from gtdb. 'release': release-str}}
    :param CLI argument parser
    :param connection: open connection to an SQlite db engine

    Return nothing
    """
    logger = logging.getLogger(__name__)

    genome_table = get_genome_table(connection)  # {genomic ver acc: {'db_id': db_id, 'gtdb_id': gtdb_id}}
    gtdb_table = get_gtdb_table_dict(connection)  # {db id: (tuple of lineage data)}
    # flip gtdb_table dict key/value pairs
    gtdb_lineage_table_dict = {}
    for db_id in gtdb_table:
        gtdb_lineage_table_dict[gtdb_table[db_id]] = db_id

    relationships_to_update = set()  # tuples (genome_id, gtdb_id)

    for genome_ver_acc in tqdm(gtdb_lineages, desc="Adding genome-gtdb relations to db"):
        try:
            genome_id = genome_table[genome_ver_acc]['db_id']
            current_gtdb_id = genome_table[genome_ver_acc]['gtdb_id']
        except KeyError:
            logger.error(
                f"Data retrieved for {genome_ver_acc}\n"
                "but genome v.accession not found in local db"
                "Not linking genome to GTDB lineages in local db"
            )
            continue

        try:
            gtdb_id = gtdb_lineage_table_dict[gtdb_lineages[genome_ver_acc]['lineage']]
        except KeyError:
            logger.error(
                "Could not find local db id for the folllowing lienage:\n"
                f"{gtdb_lineages[genome_ver_acc]['lineage']}\n"
                "Not linking lineage to genomes in local db"
            )
            continue
    
        if current_gtdb_id is None:
            relationships_to_update.add((genome_id, gtdb_id))
        else:
            if args.update_genome_lineage:
                relationships_to_update.add((genome_id, gtdb_id))

    if len(relationships_to_update):
        for record in tqdm(relationships_to_update, desc="Update genonme GTDB lineage links"):
            with connection.begin():
                connection.execute(
                    text(
                        "UPDATE Genomes "
                        f"SET gtdb_tax_id = {record[1]} "
                        f"WHERE genome_id = '{record[0]}'"
                    )
                )
    
    return

#
#
# Fix unique constrain issue when adding data to the GTDB table