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
"""Add Taxonomy data to a local SQLite database"""


import sqlite3

from cazy_webscraper.sql.sql_interface import insert_data
from cazy_webscraper.sql.sql_interface.get_data.get_table_dicts import (
    get_kingdom_table_dict,
)


def add_kingdoms(kingdoms: list, conn: sqlite3.connect) -> None:
    """Add new Kingdoms objects to database.

    Check existing kingdom objects in the db against kingdoms retrieved from the 
    CAZy txt file, so as to only add new kingdoms.

    :param cazy_taxa_dict: dict of kingdoms and organisms from the cazy_data dict
        {kingdom: {organism}}
    :param connection: open sqlalchemy connection to a local SQLite db engine

    Return nothing
    """
    kingdom_table_dict = get_kingdom_table_dict(conn)
    # dict {kingdom: {organisms}}
    existing_kingdom_records = list(kingdom_table_dict.keys())

    # create list of tuples for db insert
    kingdoms_db_insert_values = [
        (kngdm,) for kngdm in kingdoms if kngdm not in existing_kingdom_records
    ]

    if len(kingdoms_db_insert_values) != 0:
        insert_data(conn, 'Kingdoms', ['kingdom'], kingdoms_db_insert_values)

    return