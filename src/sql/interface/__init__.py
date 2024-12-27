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
"""Submodule to interact with local SQLite database, and adding data other than CAZyme records."""


import logging
import sqlite3


logger = logging.getLogger(__name__)


class SqlInterfaceException(Exception):
    """General exception for SQL interface"""

    def __init__(self, message):
        self.message = message


def insert_data(
    connection: sqlite3.Connection,
    table_name: str,
    column_names: list[str],
    insert_values: list[tuple[str]]
):
    """Insert values into one or multiple rows in the database.

    :param connection: open connection to SQLite db engine
    :param table_name: str, name of table to be inserted into
    :param column_names: list of columns (str) to insert data into
    :param insert_values: list of tuples, one tuple per inserted row in the db

    Return nothing.
    """
    logger.info("Bulk inserting data into db")

    # Set up placeholders for the VALUES statement
    placeholders = ', '.join(['?' for _ in column_names])

    query = f"INSERT INTO {table_name} ({', '.join(column_names)}) VALUES ({placeholders})"

    insert_cur = connection.cursor()
    try:
        insert_cur.executemany(query, insert_values)
        connection.commit()
    except Exception as db_error:
        connection.rollback()
        raise SqlInterfaceException(f"Database error: {str(db_error)}") from db_error
    finally:
        insert_cur.close()


# # def get_gbk_table_dict(connection):
# #     """Compile a dict of the data in the Genbanks table
    
# #     :param connection: open connection to an SQLite3 database
    
# #     Return dict {gbk accession : gbk id}
# #     """
# #     logger = logging.getLogger(__name__)

# #     logger.info("Compiling Genbank protein table into dict")

# #     with sql_orm.Session(bind=connection) as session:
# #         all_genbank = session.query(sql_orm.Genbank).all()

# #     db_gbk_dict = {}  # {genbank_accession: db genbank id number}
# #     for gbk in all_genbank:
# #         db_gbk_dict[f"{gbk.genbank_accession}"] = gbk.genbank_id
    
# #     return db_gbk_dict
