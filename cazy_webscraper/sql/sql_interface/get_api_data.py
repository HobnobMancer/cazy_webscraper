#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2021
# (c) University of Strathclyde 2021
# (c) James Hutton Institute 20201
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Retrieve protein data for provided GenBank accessions"""


import logging
import re

from tqdm import tqdm

from cazy_webscraper.sql.sql_orm import (
    CazyFamily,
    Ec,
    Genbank,
    Kingdom,
    Taxonomy,
    Session,
)


def get_class_fam_annotations(gbk_dict, query_data, connection, args):
    """Retrieve CAZy class and/or family annotations for the provided Gbks.
    
    :param gbk_dict: dict of selected GenBank accessions {acc: id}
    :param query_data: dict containing all data retrieved from the db
    :param connection: open sqlaclchemy connection for an SQLite db
    :param args: cmd-line args parser

    Return query_data: dict containing all data retrieved from the db
    """
    logger = logging.getLogger(__name__)
    
    gbk_accessions = list(gbk_dict.keys())

    # retrieve the data from the CAZy Family table
    with Session(bind=connection) as session:
        fam_table_query = session.query(Genbank.genbank_id).\
            join(CazyFamily, Genbank.families).\
            filter(Genbank.genbank_accession.in_(gbk_accessions)).\
            all()

    if len(fam_table_query) == 0:
        logger.warning(
            "No CAZy class/family annotations retrieved for any of the selected "
            "GenBank accessions."
        )
        return query_data

    for record in tqdm(fam_table_query, desc="Getting CAZy class/family annotations"):
        gbk_acc = record[0].genbank_accession

        if args.cazy_class:
            fam = record[1].family
            cazy_class = re.match(r"\D{2,3}\d", fam).group()[:-1]

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['class'].add(cazy_class)
                except KeyError:
                    query_data[gbk_acc]['class'] = {cazy_class}

            except KeyError:
                query_data[gbk_acc] = {'class': {cazy_class}}
        
        if args.cazy_family:
            fam = record[1].family

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['family'].add(fam)
                except KeyError:
                    query_data[gbk_acc]['family'] = {fam}

            except KeyError:
                query_data[gbk_acc] = {'family': {fam}}

        if args.cazy_subfamily:
            subfam = record[1].subfamily

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['subfamily'].add(subfam)
                except KeyError:
                    query_data[gbk_acc]['subfamily'] = {subfam}

            except KeyError:
                query_data[gbk_acc] = {'subfamily': {subfam}}

    return query_data



def get_tax_annotations(gbk_dict, query_data, connection, args):
    """Retrieve kingdom, genus and/or scientific name of the source organism for the provided Gbks.
    
    :param gbk_dict: dict of selected GenBank accessions {acc: id}
    :param query_data: dict containing all data retrieved from the db
    :param connection: open sqlaclchemy connection for an SQLite db
    :param args: cmd-line args parser

    Return query_data: dict containing all data retrieved from the db
    """
    logger = logging.getLogger(__name__)
    
    gbk_accessions = list(gbk_dict.keys())

    # retrieve the data from the Taxonomy and Kingdom tables
    with Session(bind=connection) as session:
        tax_query = session.query(Genbank, Taxonomy, Kingdom).\
            join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
            join(Genbank, (Genbank.taxonomy_id == Taxonomy.taxonomy_id)).\
            filter(Genbank.genbank_accession.in_(gbk_accessions)).\
            all()

    if len(tax_query) == 0:
        logger.warning("No taxonomy data retrieved for any of the selected GenBank accessions.")
        return query_data

    for record in tqdm(tax_query, desc="Getting taxonomy data"):
        gbk_acc = record[0].genbank_accession

        if args.kingdom:
            kingdom = record[2].kingdom

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['kingdom'].add(kingdom)
                except KeyError:
                    query_data[gbk_acc]['kingdom'] = {kingdom}

            except KeyError:
                query_data[gbk_acc] = {'kingdom': {kingdom}}
        
        if args.genus:
            genus = record[1].genus

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['genus'].add(genus)
                except KeyError:
                    query_data[gbk_acc]['genus'] = {genus}

            except KeyError:
                query_data[gbk_acc] = {'genus': {genus}}

        if args.organism:
            genus = record[1].genus
            species = record[1].species
            organism = f"{genus} {species}"

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['organism'].add(organism)
                except KeyError:
                    query_data[gbk_acc]['organism'] = {organism}

            except KeyError:
                query_data[gbk_acc] = {'organism': {organism}}

    return query_data
