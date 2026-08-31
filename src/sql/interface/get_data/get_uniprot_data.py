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


def get_uniprot_summary(connection: sqlite3.Connection) -> dict:
    """Compile a dict of the data in the UniProts table
    Return dict {uniprot_acc: bool if it has a sequence}
    """
    prot_cur = connection.cursor()
    prot_cur.execute("""SELECT uniprot_accession, sequence FROM UniProts""")
    db_protein_dict = {}  # {uniprot_acc: bool}
    for row in prot_cur:
        db_protein_dict[f"{row[0]}"] = True if row[1] else False
    prot_cur.close()
    return db_protein_dict


def get_uniprot_to_id(conn: sqlite3.Connection, uniprots: set[str]) -> dict[str, int]:
    cur = conn.cursor()
    cur.execute("""SELECT uniprot_id, uniprot_accession FROM UniProts""")
    uniprot2id = {}
    for row in cur:
        if row[1] in uniprots:
            uniprot2id[row[1]] = row[0]
    cur.close()
    return uniprot2id


def get_db_ecs(conn : sqlite3.Connection) -> set[str]:
    cur = conn.cursor()
    cur.execute("""SELECT ec_number FROM Ecs""")
    ec_set = set()
    for row in cur:
        ec_set.add(row[0])
    cur.close()
    return ec_set


def get_ec_to_id(conn: sqlite3.Connection, ecs: set[str]) -> set[str]:
    cur = conn.cursor()
    cur.execute("""SELECT ec_id, ec_number FROM Ecs""")
    ec2id = {}
    for row in cur:
        if row[1] in ecs:
            ec2id[row[1]] = row[0]
    cur.close()
    return ec2id


def get_db_pdbs(conn : sqlite3.Connection) -> dict[str, str]:
    cur = conn.cursor()
    cur.execute("""SELECT pdb_accession, resolution FROM Pdbs""")
    pdbs = {}
    for row in cur:
        pdbs[row[0]] = row[1]
    cur.close()
    return pdbs


def get_pdb_to_id(conn: sqlite3.Connection, pdbs: set[str]) -> set[str]:
    cur = conn.cursor()
    cur.execute("""SELECT pdb_id, pdb_accession FROM Pdbs""")
    pdb2id = {}
    for row in cur:
        if row[1] in pdbs:
            pdb2id[row[1]] = row[0]
    cur.close()
    return pdb2id


def get_db_gos(conn : sqlite3.Connection) -> set[str]:
    cur = conn.cursor()
    cur.execute("""SELECT goterm_id FROM GoTerms""")
    go_set = set()
    for row in cur:
        go_set.add(row[0])
    cur.close()
    return go_set


def get_go_to_id(conn: sqlite3.Connection, gos: set[str]) -> set[str]:
    cur = conn.cursor()
    cur.execute("""SELECT go_id, goterm_id FROM GoTerms""")
    go2id = {}
    for row in cur:
        if row[1] in gos:
            go2id[row[1]] = row[0]
    cur.close()
    return go2id
