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
from cazy_webscraper.sql.sql_interface.get_data.get_table_dicts import get_gtdb_table_dict


def add_gtdb_taxs(gtdb_lineages, connection):
    """Add GTDB lineages to GtdbTaxs table

    :param gtdb_lineags:
    :param connection: open connection to an SQlite db engine

    Return nothing
    """
    # loading existing table data
    existing_gtdb_table = get_gtdb_table_dict(connection)

    lineages_to_add = set()

    for genome in tqdm(gtdb_lineages, desc='Adding GTDB lineages to db'):
        lineage = [" ".join(_.split("__")[1:]) for _ in gtdb_lineages[genome].split(";")]
        if lineage not in list(existing_gtdb_table.values()):
            lineages_to_add.add(lineage)

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
                'strain',
            ],
            list(lineages_to_add),
        )

    return
