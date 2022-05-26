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
"""Retrieve proteins with no assembly data in the local database."""


from tqdm import tqdm

from cazy_webscraper.sql.sql_orm import (
    Genbank,
    Genome,
    Session,
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
