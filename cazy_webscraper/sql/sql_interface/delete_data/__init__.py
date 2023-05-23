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
"""Delete old records in the local CAZyme database."""


import logging
from sqlalchemy import delete, text
from tqdm import tqdm

from cazy_webscraper.sql.sql_orm import Ec, Pdb


def delete_old_relationships(
    uniprot_dict,
    gbk_dict,
    annotation_table,
    relationship_table,
    annotation,
    table_name,
    connection,
    args,
):
    """Delete old relationships that are in the local CAZyme database between a protein and its
    annotations but are not present in the data downloaded today from a remote database, 
    e.g. UniProt.

    :param uniprot_dict: dict containing data retrieved from UniProt
        uniprot_dict[ncbi_acc] = {
            'uniprot_acc': uniprot_acc, - str
            'protein_name': protein_name, -str
            'ec_numbers': ec_numbers, - set
            'sequence': sequence, - str
            'sequence_date': date seq was last updated yyyy-mm-dd
            'pdbs': all_pdbs, - set
        }
    :param gbk_dict: {gbk_acc: gbk_id}, data in the Genbanks table
    :param annotation_table: output from get_pdb_talbe or get_ec_table
        {pdb_accession: pdb_db_id} or {ec_number: ec_id}
    :param relationship_table: dict, data in Genbanks_Pdbs or Genbanks_Ec table
        {gbk_db_id: set(pdb_db_ids) } or {ec_id: set(gbk ids)}
    :param annotation: str, Type of annotation, e.g. ec_numbers or pdbs
        - matches key in uniprot_dict!
    :param table_name: str, name of relationship table in the local db
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing
    """
    # set of tuples (gbk_db_id, annotation_db_id)  # e.g. ec_db_id or pdb_db_id
    relationships_to_del  = set(), set() 

    for ncbi_acc in tqdm(uniprot_dict, desc=f"Identifying {annotation} relationships to delete"):
        gbk_db_id = gbk_dict[ncbi_acc]

        # get all presently linked annotaitons, e.g. ec_numbers or pdb_accessions
        # in the local CAZyme db
        current_annotations = set()  # ec or pdb db ids
        if annotation == 'ec_numbers':
            for ec_id in relationship_table:
                if gbk_db_id in relationship_table[ec_id]:
                    current_annotations.add(ec_id)

        else:
            try:
                current_annotations = relationship_table[gbk_db_id]
            except KeyError:
                continue # gbk/ncbi acc not in the relationship table

        if len(current_annotations) == 0:
            continue # gbk/ncbi acc not in the relationship table

        # identify annotations retrieved today, i.e. the new annotations
        new_annotations = set()
        for anno_acc in uniprot_dict[ncbi_acc][annotation]:
            # get the local db id for the annotation
            anno_db_id = annotation_table[anno_acc]
            new_annotations.add(anno_db_id)

        for anno_db_id in new_annotations:
            if anno_db_id not in current_annotations:
                relationships_to_del.add( (gbk_db_id, anno_db_id) )

    # delete the old relationships and annotations
    if len(relationships_to_del) != 0:
        logger.warning(
            f"Deleteing {len(relationships_to_del)} relationships from the local db {table_name} table"
        )
        if annotation == 'ec_numbers':
            anno_id_name = 'ec_id'
        else:
            anno_id_name = 'pdb_id'
        for record in tqdm(relationships_to_del, desc=f"Deleting old relationships in table {table_name}"):
            # tuple (gbk_db_id, anno_db_id)
            connection.execute(
                text(
                    f"DELETE FROM {table_name} WHERE genbank_id='{record[0]}' AND {anno_id_name}='{record[1]}'"
                )
            )

def delete_old_annotations(annotation_table, relationship_table, table_name, connection, args):
    """Delete old annotations (EC numbers, PDB accessions) that are not linked to any proteins in
    the Genbanks table.

    :param relationship_table: dict, data in Genbanks_Pdbs or Genbanks_Ec table
        {gbk_db_id: set(pdb_db_ids) } or {ec_id: set(gbk ids)}
    :param annotation: str, Type of annotation, e.g. ec_numbers or pdbs
        - matches key in uniprot_dict!
    :param table_name: str, name of relationship table in the local db
    :param connection: open sqlalchemy conenction to an SQLite db engine
    :param args: cmd-line args parser

    Return nothing
    """
    # identify all annotations that are linked to a protein the relationships talbe
    all_linked_annotations = set()
    
    if table_name == 'Ecs':
        all_linked_annotations = set(relationship_table.keys())
    else: 
        for gbk_db_id in relationship_table:
            all_linked_annotations = all_linked_annotations.union(relationship_table[gbk_db_id])

    if len(all_linked_annotations) == 0:
        return
    
    # identify annotations in the annotation table that are not linked to any records in the 
    # relationship table
    unlinked_annotations = set()
    for anno_acc in annotation_table:
        anno_db_id = annotation_table[anno_acc]
        if anno_db_id not in all_linked_annotations:
            unlinked_annotations.add(anno_db_id)
    
    if len(unlinked_annotations) == 0:
        logger.warning(f"No unlinked records found in table {table_name}")
        return

    if table_name == 'Ecs':
        anno_db_id_name = 'ec_id'
    else:
        anno_db_id_name = 'pdb_id'
    
    logger.warning(
        f"Identified {anno_db_id_name}s in the local CAZyme db {table_name} table "
        "that are not linked to any proteins in the Genbanks table\n"
        f"Deleting these {anno_db_id_name}s"
    )
    for anno_db_id in tqdm(unlinked_annotations, desc=f"Deleting unlinked {anno_db_id_name}s"):
            connection.execute(
                text(f"DELETE FROM {table_name} WHERE {anno_db_id_name}='{anno_db_id}'")
            )
