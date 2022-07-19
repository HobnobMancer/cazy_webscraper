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
"""Add protein sequences retrieved from GenBank to a local SQLite database"""


import logging
from sqlalchemy import text
from tqdm import tqdm

from cazy_webscraper.sql.sql_interface.get_data import get_table_dicts


def add_gbk_seqs_to_db(seq_dict, retrieval_date, gbk_dict, connection, args):
    """Add sequences retrieved from GenBank to the local CAZyme db
    
    :param seq_dict: dict of retrieved seq {gbk_acc: seq}
    :param retrieval_date: str, date seqs were retrieved in ISO format
    :param gbk_dict: dict {gbk_acc: db_id}
    :param connection: open sqlalchemy connection to a SQLite db
    :param args: cmd-line args parser
    
    Return nothing
    """
    logger = logging.getLogger(__name__)
    
    # load the current Genbanks table into a dict
    gbk_table_dict = get_table_dicts.get_gbk_table_seq_dict(connection)

    records = set()  # records to update, contains tuples (gbk_acc, seq, )

    for gbk_accession in tqdm(seq_dict, desc="Compiling records for db insertion"):
        existing_record = gbk_table_dict[gbk_accession]
        seq = seq_dict[gbk_accession]
        db_id = gbk_dict[gbk_accession]

        if existing_record['sequence'] is None:
            records.add( (db_id, seq) )
        else:
            if (args.seq_update) and (retrieval_date < existing_record['seq_date']):
                records.add( (db_id, seq) )

    if len(records) != 0:
        with connection.begin():
            for record in tqdm(records, desc="Adding seqs to db"): 
                    connection.execute(
                        text(
                            "UPDATE Genbanks "
                            f"SET sequence = '{record[1]}'"
                            f"WHERE genbank_id = '{record[0]}'"
                        )
                    )
                    #, seq_update_date = {retrieval_date} 
                    connection.execute(
                        text(
                            "UPDATE Genbanks "
                            f"SET seq_update_date = '{retrieval_date}'"
                            f"WHERE genbank_id = '{record[0]}'"
                        )
                    )
    else:
        logger.warning("No records to update")

    return
