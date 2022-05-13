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
    Pdb,
    Taxonomy,
    Session,
    Uniprot,
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
        fam_table_query = session.query(Genbank, CazyFamily).\
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

        if 'class' in args.include:
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
        
        if 'family' in args.include:
            fam = record[1].family

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['family'].add(fam)
                except KeyError:
                    query_data[gbk_acc]['family'] = {fam}

            except KeyError:
                query_data[gbk_acc] = {'family': {fam}}

        if 'subfamily' in args.include:
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

        if 'kingdom' in args.include:
            kingdom = record[2].kingdom

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['kingdom']
                    logger.warning(
                        f"Multiple taxa found for {gbk_acc}\n"
                        "Retreiving only one record."
                    )
                    query_data[gbk_acc]['kingdom'] = kingdom
                except KeyError:
                    query_data[gbk_acc]['kingdom'] = kingdom

            except KeyError:
                query_data[gbk_acc] = {'kingdom': kingdom}
        
        if 'genus' in args.include:
            genus = record[1].genus

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['genus']
                    logger.warning(
                        f"Multiple taxa found for {gbk_acc}\n"
                        "Retreiving only one record."
                    )
                    query_data[gbk_acc]['genus'] = genus
                except KeyError:
                    query_data[gbk_acc]['genus'] = genus

            except KeyError:
                query_data[gbk_acc] = {'genus': genus}

        if 'organism' in args.include:
            genus = record[1].genus
            species = record[1].species
            organism = f"{genus} {species}"

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['organism']
                    logger.warning(
                        f"Multiple taxa found for {gbk_acc}\n"
                        "Retreiving only one record."
                    )
                    query_data[gbk_acc]['organism'] = organism
                except KeyError:
                    query_data[gbk_acc]['organism'] = organism

            except KeyError:
                query_data[gbk_acc] = {'organism': organism}

    return query_data


def get_ec_annotations(gbk_dict, query_data, connection):
    """Retrieve EC number annotations for the provided Gbks.
    
    :param gbk_dict: dict of selected GenBank accessions {acc: id}
    :param query_data: dict containing all data retrieved from the db
    :param connection: open sqlaclchemy connection for an SQLite db

    Return query_data: dict containing all data retrieved from the db
    """
    logger = logging.getLogger(__name__)
    
    gbk_accessions = list(gbk_dict.keys())

    # retrieve the data from the Taxonomy and Kingdom tables
    with Session(bind=connection) as session:
        ec_query = session.query(Genbank, Ec).\
            join(Ec, Genbank.ecs).\
            filter(Genbank.genbank_accession.in_(gbk_accessions)).\
            all()

    if len(ec_query) == 0:
        logger.warning("No EC annotations retrieved for any of the selected GenBank accessions.")
        return query_data

    for record in tqdm(ec_query, desc="Getting EC number annotations"):
        gbk_acc = record[0].genbank_accession

        ec_number = record[1].ec_number

        try:
            query_data[gbk_acc]
            
            try:
                query_data[gbk_acc]['ec_numbers'].add(ec_number)
            except KeyError:
                query_data[gbk_acc]['ec_numbers'] = {ec_number}

        except KeyError:
            query_data[gbk_acc] = {'ec_numbers': {ec_number}}

    return query_data


def get_pdb_accessions(gbk_dict, query_data, connection):
    """Retrieve PDB accessions for the provided Gbks.
    
    :param gbk_dict: dict of selected GenBank accessions {acc: id}
    :param query_data: dict containing all data retrieved from the db
    :param connection: open sqlaclchemy connection for an SQLite db

    Return query_data: dict containing all data retrieved from the db
    """
    logger = logging.getLogger(__name__)
    
    gbk_accessions = list(gbk_dict.keys())

    # retrieve the data from the Taxonomy and Kingdom tables
    with Session(bind=connection) as session:
        pdb_query = session.query(Genbank, Pdb).\
            join(Pdb, Genbank.pdbs).\
            filter(Genbank.genbank_accession.in_(gbk_accessions)).\
            all()

    if len(pdb_query) == 0:
        logger.warning("No PDB accessions retrieved for any of the selected GenBank accessions.")
        return query_data

    for record in tqdm(pdb_query, desc="Getting PDB accessions"):
        gbk_acc = record[0].genbank_accession

        pdb_accession = record[1].pdb_accession

        try:
            query_data[gbk_acc]
            
            try:
                query_data[gbk_acc]['pdb_accessions'].add(pdb_accession)
            except KeyError:
                query_data[gbk_acc]['pdb_accessions'] = {pdb_accession}

        except KeyError:
            query_data[gbk_acc] = {'pdb_accessions': {pdb_accession}}

    return query_data


def get_uniprot_data(gbk_dict, query_data, connection, args):
    """Retrieve UniProt data for the provided Gbks.
    
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
        uniprot_query = session.query(Genbank, Uniprot).\
            join(Uniprot, (Uniprot.genbank_id == Genbank.genbank_id)).\
            filter(Genbank.genbank_accession.in_(gbk_accessions)).\
            all()

    if len(uniprot_query) == 0:
        logger.warning("No UniProt records retrieved for any of the selected GenBank accessions.")
        return query_data

    for record in tqdm(uniprot_query, desc="Getting UniProt data"):
        gbk_acc = record[0].genbank_accession

        if 'uniprot_acc' in args.include:
            uniprot_accession = record[1].uniprot_accession

            try:
                query_data[gbk_acc]
    
                try:
                    query_data[gbk_acc]['uniprot_accession']
                    logger.warning(
                        f"Multiple UniProt records found for GBK acc {gbk_acc}\n"
                        "Retreiving only one."
                    )
                    query_data[gbk_acc]['uniprot_accession'] = uniprot_accession

                except KeyError:
                    query_data[gbk_acc]['uniprot_accession'] = uniprot_accession

            except KeyError:
                query_data[gbk_acc] = {
                    'uniprot_accession': uniprot_accession,
                }

        if 'uniprot_name' in args.include:
            uniprot_name = record[1].uniprot_name

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['uniprot_name']
                    logger.warning(
                        f"Multiple UniProt records found for GBK acc {gbk_acc}\n"
                        "Retreiving only one."
                    )
                    query_data[gbk_acc]['uniprot_name'] = uniprot_name
                    
                except KeyError:
                    query_data[gbk_acc]['uniprot_name'] = uniprot_name

            except KeyError:
                query_data[gbk_acc] = {
                    'uniprot_name': uniprot_name,
                }

        if 'uniprot_seq' in args.include:
            seq = record[1].sequence
            seq_date = record[1].seq_update_date

            try:
                query_data[gbk_acc]
                
                try:
                    query_data[gbk_acc]['uniprot_sequence']
                    logger.warning(
                        f"Multiple UniProt records found for GBK acc {gbk_acc}\n"
                        "Retreiving only one record."
                    )
                    query_data[gbk_acc]['uniprot_sequence'] = seq
                    query_data[gbk_acc]['uniprot_sequence_date'] = seq_date
                except KeyError:
                    query_data[gbk_acc]['uniprot_sequence'] = seq
                    query_data[gbk_acc]['uniprot_sequence_date'] = seq_date

            except KeyError:
                query_data[gbk_acc] = {'sequence': seq, 'sequence_date': seq_date}

    return query_data


def get_gbk_seq(gbk_dict, query_data, connection):
    """Retrieve GenBank protein sequences for the provided Gbks.
    
    :param gbk_dict: dict of selected GenBank accessions {acc: id}
    :param query_data: dict containing all data retrieved from the db
    :param connection: open sqlaclchemy connection for an SQLite db

    Return query_data: dict containing all data retrieved from the db
    """
    logger = logging.getLogger(__name__)
    
    gbk_accessions = list(gbk_dict.keys())

    # retrieve the data from the Taxonomy and Kingdom tables
    with Session(bind=connection) as session:
        gbk_query = session.query(Genbank).\
            filter(Genbank.genbank_accession.in_(gbk_accessions)).\
            all()

    if len(gbk_query) == 0:
        logger.warning("No GenBank records retrieved for any of the selected GenBank accessions.")
        return query_data

    for record in tqdm(gbk_query, desc="Getting GenBank protein sequences"):
        gbk_acc = record[0].genbank_accession
        seq = record[0].sequence
        seq_date = record[0].seq_update_date

        try:
            query_data[gbk_acc]
            
            try:
                logger.warning(
                    f"Multiple GBK records found for GBK acc {gbk_acc}\n"
                    "Retreiving only one gbk sequence."
                )
                query_data[gbk_acc]['gbk_sequence'] = seq
                query_data[gbk_acc]['gbk_sequence_date'] = seq_date
            except KeyError:
                query_data[gbk_acc]['gbk_sequence'] = seq
                query_data[gbk_acc]['gbk_sequence_date'] = seq_date

        except KeyError:
            query_data[gbk_acc] = {'gbk_sequence': seq, 'gbk_sequence_date': seq_date}

    return query_data

