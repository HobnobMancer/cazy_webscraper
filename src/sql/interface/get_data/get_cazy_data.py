#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2025
# (c) University of Strathclyde 2025
# (c) James Hutton Institute 2025
# Author:
# Emma E. M. Hobbs
#
# Contact
# ehobbs@ebi.ac.uk
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


import sqlite3



def get_fams_table_dict(connection: sqlite3.Connection) -> dict:
    """Create dict of objects present in the CazyFamilies table.
    Return dict {family subfamily: db_family_id}
    """
    fam_cur = connection.cursor()
    fam_cur.execute("""SELECT * FROM CazyFamilies""")
    db_fam_dict = {}
    for row in fam_cur:
        # [0] fam_id, [1] fam, [2] subfamily
        subfam = row[2] if not None else '_'
        db_fam_dict[f"{row[1]} {subfam}"] = row[0]
    fam_cur.close()
    return db_fam_dict


def get_prot_fam_table_dict(connection: sqlite3.Connection) -> dict:
    """Build dict representing the records present in the Proteins_CazyFamilies table

    If a GenBank accession is in the db but not has not CazyFamilies instances related to it,
    the GenBank accession is not returned when quering the db.

    Return {protein_id: {family_id}}
    """
    prot_fam_cur = connection.cursor()
    prot_fam_cur.execute("""SELECT * FROM Proteins_CazyFamilies""")

    prot_fam_table_dict = {}  # {protein_id: {family_id}}
    for row in prot_fam_cur:
        # [0] protein_id, [1] family_id
        if row[0] not in prot_fam_table_dict:
            prot_fam_table_dict[row[0]] == set()
        prot_fam_table_dict[row[0]].add(row[1])
    prot_fam_cur.close()

    return prot_fam_table_dict


def get_kingdom_table_dict(connection: sqlite3.Connection) -> dict:
    """Load and parse the Kingdoms table from the db and compile a dict {kgnd: id}

    Return dict {kingdom: kindom_db_id}
    """
    king_cur = connection.cursor()
    king_cur.execute("""SELECT * FROM Kingdoms""")
    kingdom_dict = {}  # {kingdom: kindom_db_id}
    for row in king_cur:
        # row0 = kingdom_id; row1 = kingdom
        kingdom_dict[row[1]] = row[0]
    king_cur.close()
    return kingdom_dict


def get_taxs_table_dict(connection: sqlite3.Connection) -> dict:
    """Create dict of objects present in the Taxs table.

    Return dict {genus species: {'tax_id': db_tax_id, 'kingdom_id': kingdom_id}
    """
    tax_cur = connection.cursor()
    tax_cur.execute("""SELECT * FROM Taxs""")
    db_tax_dict = {}
    for row in tax_cur:
        # [0] tax id, [1] genus, [2] species, [3] kingdom id
        db_tax_dict[f"{row[1]} {row[2]}"] = {
            'tax_id': row[0],
            'kingdom_id': row[3]
        }
    tax_cur.close()    
    return db_tax_dict
