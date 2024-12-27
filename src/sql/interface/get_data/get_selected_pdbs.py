#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
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
#
# Bio.PDB reference:
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
# implemented in Python. Bioinformatics 19: 2308–2310
"""Retrieve PDB accessions matching user criteria from the local CAZyme db"""


from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks, get_table_dicts


def get_pdb_accessions(
    class_filters,
    family_filters,
    taxonomy_filters,
    kingdom_filters,
    ec_filters,
    gbk_table_dict,
    connection,
):
    """Retrieve PDB accessions matching user criteria from the local CAZyme db
    
    :param class_filters: set of CAZy classes to retrieve data for
    :param family_filters: set of CAZy families to retrieve data for
    :param taxonomy_filters: dict of taxonom filters to limit the retrieval of data to
    :param kingdom_filters: set of tax kingdoms to limit the retrieval of data to
    :param ec_filters: set of EC numbers to limit the retrieval of data to
    :param gbk_table_dict: dict, of GenBank accessions and db IDs for GenBank records
    :param connection: open sqlalchemy connection to an SQLite db engine
    
    Return list of PDB accessions
    """

    # retrieve the GenBank accessions of Gbk records matching the user criteria
    selected_gbks = get_selected_gbks.get_genbank_accessions(
        class_filters,
        family_filters,
        taxonomy_filters,
        kingdom_filters,
        ec_filters,
        connection,
    )
    
    # retrieve all PDB accessions for each GenBank accession retrieved for the local CAZyme db

    pdb_table_dict = get_table_dicts.get_pdb_table_dict(connection)  # used to retrieve PDB accs
    # {pdb_accession: pdb_db_id}
    # convert to be keyed by pdb_db_id and valued by pdb_accession
    pdb_id_acc_dict = {}
    for pdb_acc in pdb_table_dict:
        pdb_id_acc_dict[pdb_table_dict[pdb_acc]] = pdb_acc

    gbk_pdb_table_dict = get_table_dicts.get_gbk_pdb_table_dict(connection)
    # {gbk_db_id: {pdb_db_id}}

    pdb_accessions = set()

    for gbk_accession in selected_gbks:
        gbk_id = gbk_table_dict[gbk_accession]['gbk_id']

        try:
            pdb_ids = gbk_pdb_table_dict[gbk_id]

            for pdb_id in pdb_ids:
                # convert pdb_db_id to pdb_accession
                pdb_acc = pdb_id_acc_dict[pdb_id]
                # remove chain annotation if present
                parent_pdb = pdb_acc.split("[")[0]
                pdb_accessions.add(parent_pdb)
        
        except KeyError:
            continue

    return list(pdb_accessions)
