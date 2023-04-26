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
"""Add NCBI Taxonomy data to a local SQLite database"""


import logging
from sqlalchemy import text
from tqdm import tqdm

from cazy_webscraper.sql.sql_interface import insert_data
from cazy_webscraper.sql.sql_interface.get_data.get_table_dicts import get_ncbi_tax_table, get_no_tax_gbk_table_dict


def add_ncbi_taxonomies(tax_dict, connection, args):
    """Add NCBI taxonomy data to the local db

    :param tax_dict: dict, keyed by ncbi tax db id, and valued by linegae data
    :param connection: open connection to an sqlite db
    :param args: cmd-line args parser

    Return nothing
    """
    logger = logging.getLogger(__name__)

    # load ncbiTax table into dict
    db_ncbi_tax_table = get_ncbi_tax_table(connection)  # {ncbi_tax_id: local db id}

    records_to_add = set()
    records_to_update = set()

    for ncbi_tax_id in tax_dict:
        tax_data = (
            ncbi_tax_id,
            tax_dict[ncbi_tax_id]['kingdom'],
            tax_dict[ncbi_tax_id]['phylum'],
            tax_dict[ncbi_tax_id]['class'],
            tax_dict[ncbi_tax_id]['order'],
            tax_dict[ncbi_tax_id]['family'],
            tax_dict[ncbi_tax_id]['genus'],
            tax_dict[ncbi_tax_id]['species'],
            tax_dict[ncbi_tax_id]['strain'],
        )

        try:
            db_ncbi_tax_table[int(ncbi_tax_id)]
            if args.update_taxs:
                records_to_update.add(tax_data)

        except KeyError:
            records_to_add.add(tax_data)

    if args.update_taxs and len(records_to_update) != 0:
        # update tax records
        update_ncbi_tax_records(records_to_update, connection)

    if len(records_to_add) != 0:
        # add records
        insert_data(
            connection,
            'NcbiTaxs',
            [
                'ncbi_tax_id',
                'kingdom',
                'phylum',
                'tax_class',
                'tax_order',
                'family',
                'genus',
                'species',
                'strain',
            ],
            list(records_to_add),
        )

    return


def update_ncbi_tax_records(records_to_update, connection, unit_test=False):
    """Update ncbi tax records in the local CAZyme database

    :param records_to_update: set of tuples, one tuple per record to update
    :param connection: open connection to an sqlite db

    Return nothing
    """
    for record in tqdm(records_to_update, desc="Updating NCBI tax records in the local db"):
        with connection.begin():
            connection.execute(
                text(
                    "UPDATE NcbiTaxs "
                    f"SET kingdom = '{record[1]}' AND "
                    f"phylum = '{record[2]}' AND "
                    f"tax_class = '{record[3]}' AND "
                    f"tax_order = '{record[4]}' AND "
                    f"family = '{record[5]}' AND "
                    f"genus = '{record[6]}' AND "
                    f"species = '{record[7]}' AND "
                    f"strain = '{record[8]}' "
                    f"WHERE ncbi_tax_id = '{record[0]}'"
                )
            )
            if unit_test:
                connection.rollback()


def update_genbank_ncbi_tax(tax_prot_dict, connection, args, unit_test=False):
    """Update tax information in the Genbanks table

    :param tax_prot_dict: dict {tax_id: lineage and protein data}
    :param connection: open connection to a sql db
    :param args: cmd-line args parser

    Return nothing
    """
    logger = logging.getLogger(__name__)
    
    db_ncbi_tax_table = get_ncbi_tax_table(connection)  # {ncbi_tax_id: local db id}

    if args.update_gbk:
        with connection.begin():
            for tax_id in tqdm(tax_prot_dict, desc="Updating Genbanks table"):
                try:
                    proteins = tax_prot_dict[int(tax_id)]['proteins']
                except KeyError:
                    try:
                        proteins = tax_prot_dict[tax_id]['proteins']
                    except KeyError:
                        logger.warning(
                            f"Could not retrieve proteins linked to tax id {tax_id}\n"
                            "Will not update the respective records in the Genbanks table"
                        )
                        continue
                try:    
                    tax_db_id = db_ncbi_tax_table[int(tax_id)]
                except KeyError:
                    try:
                        tax_db_id = db_ncbi_tax_table[str(tax_id)]
                    except KeyError:
                        logger.warning(
                            f"Could not retrieve the local db ID for the NCBI tax id {tax_id}\n"
                            "Will not update the respective records in the Genbanks table"
                        )
                        continue
                for prot_db_id in proteins:
                    connection.execute(
                        text(
                            "UPDATE Genbanks "
                            f"SET ncbi_tax_id = {tax_db_id} "
                            f"WHERE genbank_id = '{prot_db_id}'"
                        )
                    )
                if unit_test:
                    connection.rollback()

    else:
        # filter for protein records with no tax data, and only add new tax data, do no overwrite existing data
        gbk_db_ids = get_no_tax_gbk_table_dict(connection)

        for tax_id in tqdm(tax_prot_dict, desc="Updating Genbanks table"):
            tax_db_id = db_ncbi_tax_table[int(tax_id)]
            proteins = tax_prot_dict[tax_id]['proteins']
            for prot_db_id in proteins:
                if prot_db_id in gbk_db_ids:
                    with connection.begin():
                        connection.execute(
                            text(
                                "UPDATE Genbanks "
                                f"SET ncbi_tax_id = {tax_db_id} "
                                f"WHERE genbank_id = {prot_db_id}"
                            )
                        )
                        if unit_test:
                            connection.rollback()

    return
