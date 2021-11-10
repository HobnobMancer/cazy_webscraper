#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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

from tqdm import tqdm

from cazy_webscraper.sql import sql_orm


class SqlInterfaceException(Exception):
    """General exception for SQL interface"""

    def __init__(self, message):
        self.message = message


def log_scrape_in_db(
    time_stamp,
    config_dict,
    taxonomy_filters,
    kingdoms,
    ec_filter,
    db,
    retrieved_annotations,
    session,
    args,
):
    """Add a log of scraping CAZy to the local database.

    :param time_stamp: str, date and time cazy_webscraper was invoked
    :param config_dict: dict of CAZy classes and families to be scraped
    :param taxonomy_filters: dict of genera, species and strains to restrict the scrape to
    :param kingdoms: list of taxonomy Kingdoms to restrict scrape to
    :param ec_filter: set of EC numbers to restrict retrieval of data to
    :param db: str, name of the external database from which data is retrieved
    :param retrieved_annotations: str, types of annotations retrieved (e.g. UniProt accessions)
    :param session: open SQL database session
    :param args: cmd arguments

    Return nothing."""
    logger = logging.getLogger(__name__)

    date = time_stamp[:time_stamp.find("--")]
    time = time_stamp[((time_stamp.find("--")) + 2):].replace("-", ":")

    new_log = sql_orm.Log(
        date=date,
        time=time,
        database=db,
        retrieved_annotations=retrieved_annotations,
    )

    if config_dict is not None:
        # get classes that user named to be scraped
        try:
            classes = config_dict["classes"]
            if classes is not None:
                classes = str(classes).replace("[", "").replace("]", "").replace("'", "")
                new_log.classes = classes
        except KeyError:
            pass

        # create a list of families instructed to be scraped
        families = []
        for key in config_dict:
            if (key != "classes") and (config_dict[key] is not None):
                families.append(config_dict[key])

        if len(families) != 0:
            families = str(families).replace("[", "").replace("]", "").replace("'", "")
            new_log.families = families

    # get taxonomy filters defined by user, and separate into genera, species and strains
    try:
        if len(taxonomy_filters["genera"]) != 0:
            genera = str(taxonomy_filters["genera"]).replace("[", "").\
                replace("]", "").replace("'", "")
            new_log.genera = genera
    except TypeError:
        pass

    try:
        if len(taxonomy_filters["species"]) != 0:
            species = str(taxonomy_filters["species"])
            species = species.replace("[", "").replace("]", "").replace("'", "")
            new_log.species = species
    except TypeError:
        pass

    try:
        if len(taxonomy_filters["strains"]) != 0:
            strains = str(taxonomy_filters["strains"])
            strains = strains.replace("[", "").replace("]", "").replace("'", "")
            new_log.strains = strains
    except TypeError:
        pass

    # get Taxonomy Kingdoms defined by user to be scraped
    if kingdoms is not None:
        new_log.kingdoms = str(kingdoms).replace("[", "").replace("]", "").replace("'", "")
    else:
        new_log.kingdoms = "ALL (Archaea, Bacteria, Eukaryota, Viruses, Unclassified"
    
    # retrieve commands from the command line
    cmd_line = ""
    for cmd in [
        [args.classes, " --classes '"],
        [args.families, " --families '"],
        [args.kingdoms, " --kingdoms"],
        [args.genera, " --genera '"],
        [args.species, " --species '"],
        [args.strains, " --strains '"],
    ]:
        try:
            cmd_line = cmd_line + cmd[1] + cmd[0] + "'"
        except TypeError:
            pass

    if len(ec_filter) != 0:
        cmd_line = cmd_line + "--ec_filter" + (args.ec) + "'"
        new_log.ec_filter = ','.join(list(ec_filter))

    if len(cmd_line) != 0:
        cmd_line = cmd_line.strip()
        new_log.cmd_line = cmd_line

    session.add(new_log)
    session.commit()

    return


def insert_data(connection, table_name, column_names, insert_values):
    """Insert values into one or multiple rows in the database.
    
    :param connection: open connection to SQLite db engine
    :param table_name: str, name of table to be inserted into
    :param column_names: list of columns (str) to insert data into
    :param insert_values: list of tuples, one tuple per inserted row in the db
    
    Return nothing.
    """
    # set up series of ? to fill in the VALUES statement
    value_stmt = ''
    for name in range((len(column_names)) - 1):
        value_stmt += '?, '
    value_stmt += '?'  # statement should not end with a comma
    
    with connection.begin():
        try:
            connection.exec_driver_sql(
                f"INSERT INTO {table_name} ({', '.join(column_names)}) VALUES ({value_stmt})",
                insert_values,
            )
        except Exception as db_error:
            raise SqlInterfaceException(db_error)

    return


def get_gbk_table_dict(connection):
    """Compile a dict of the data in the Genbanks table
    
    :param connection: open connection to an SQLite3 database
    
    Return dict {gbk accession : gbk id}
    """
    with sql_orm.Session(bind=connection) as session:
        all_genbank = session.query(sql_orm.Genbank).all()

    db_gbk_dict = {}  # {genbank_accession: db genbank id number}
    for gbk in all_genbank:
        db_gbk_dict[f"{gbk.genbank_accession}"] = gbk.genbank_id
    
    return db_gbk_dict
