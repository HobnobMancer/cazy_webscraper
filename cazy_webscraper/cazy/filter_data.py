#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Filter data in CAZy database dump (txt file)."""


import logging
import re
import sqlite3

from pathlib import Path


logger = logging.getLogger(__name__)


def apply_tax_filters(
    kingdoms: set[str],
    genera: set[str],
    species: set[str],
    strains: set[str],
    db: Path
):
    """
    :param kingdom filter: set of kingdoms
    :param taxonomy_filter_dict: dict
    E.g. {
        'genera': {'Aspergillus', 'AnotherGenus},
        'species': {'Bacteroides cellulosilyticus', 'Genus species'},
        'strains': {'Alternaria alternata SRC1lrK2f v1.0', 'Genus species strain'}
    }
    """
    query = "DELETE FROM TempTable WHERE"
    parameters = []

    if kingdoms:
        query += " kingdom NOT IN ({})".format(', '.join('?' for _ in kingdoms))
        parameters += list(kingdoms)

    if genera:
        if kingdoms:
            query += " AND"
        query += " genus NOT IN ({})".format(','.join('?' for _ in genera))
        parameters += list(genera)

    if species:
        if parameters:
            query += " AND"
        query += f" (species NOT LIKE {' OR '.join(['?'] * len(species))})"
        parameters.extend([f'{sp}%' for sp in species])

    if strains:
        if parameters:
            query += " AND"
        query += " species NOT IN ({})".format(','.join('?' for _ in strains))
        parameters.extend(strains)

    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute(query, parameters)
    cur.close()
    conn.commit()
    conn.close()


def apply_class_and_family_filters(
    class_filter: list[str],
    fam_filter: set[str],
    db: Path,
):
    query = "DELETE FROM TempTable WHERE"

    for i, cazy_class in enumerate(class_filter):
        if i == 0:
            query += " ("
        if i > 0:
            query += " AND "
        query += f"family NOT LIKE '{cazy_class}%'"

    if class_filter and fam_filter:
        query += ") AND"

    for i, fam in enumerate(fam_filter):
        if i == 0:
            query += " ("
        if i > 0:
            query += " AND "

        if re.match(r"^\D{2,3}\d+_\d$", fam):  # subfam, only keep specific subfam
            query += f"family NOT LIKE '{fam}'"
        else:  # keep fam and all it's subfamilies
            query += f"family NOT LIKE '{fam}'"
            query += f" AND family NOT LIKE '{fam}_%'"

    query += ")"

    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute(query)
    conn.commit()
    cur.close()
    conn.close()


def drop_subfamilies(db: Path):
    query = "DELETE FROM TempTable WHERE family LIKE '%\_%' ESCAPE '\'"
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute(query)
    conn.commit()
    cur.close()
    conn.close()
