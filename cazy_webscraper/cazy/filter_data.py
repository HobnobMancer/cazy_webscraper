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
import sqlite3

from pathlib import Path


logger = logging.getLogger(__name__)


def apply_kingdom_filers(kingdom_filter: set[str], db: Path):
    logger.warning("Filtering to kingdoms: %s", kingdom_filter)
    query = "DELETE FROM TempTable WHERE kingdom NOT IN ({})".format(
        ', '.join('?' for _ in kingdom_filter)
    )

    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute(query, list(kingdom_filter))
    conn.commit()
    conn.close()


def apply_tax_filters(genera: set[str], species: set[str], strains: set[str], db: Path):
    """
    :param taxonomy_filter_dict: dict
    E.g. {
        'genera': {'Aspergillus', 'AnotherGenus},
        'species': {'Bacteroides cellulosilyticus', 'Genus species'},
        'strains': {'Alternaria alternata SRC1lrK2f v1.0', 'Genus species strain'}
    }
    """
    query = "DELETE FROM TempTable WHERE "

    if genera:
        query += "genus NOT IN ({})".format(','.join('?' for _ in genera))
        parameters = list(genera)
    else:
        parameters = []

    if species:
        species_clauses = ["species NOT LIKE ?" for _ in species]
        if parameters:
            query += " AND "
        query += " AND ".join(species_clauses)
        species_with_wildcards = [s + '%' for s in species]
        parameters.extend(species_with_wildcards)

    if strains:
        if parameters:
            query += " AND "
        query += "species NOT IN ({})".format(','.join('?' for _ in strains))
        parameters += list(strains)

    print(query)
    print(parameters)

    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute(query, parameters)
    conn.commit()
    conn.close()


def apply_class_and_family_filters(excluded_classes: list[str], fam_filter: set[str], db: Path):
    if excluded_classes and fam_filter:
        class_query = " AND ".join(f"family NOT LIKE '{class_}%'" for class_ in excluded_classes)
        fam_query = ', '.join(f"'{family}'" for family in fam_filter)
        query = f"""
            DELETE FROM TempTable
            WHERE {class_query}
            AND family NOT IN ({fam_query});
        """
    elif excluded_classes and not fam_filter:
        class_query = " AND ".join(f"family NOT LIKE '{class_}%'" for class_ in excluded_classes)
        query = f"""
            DELETE FROM TempTable
            WHERE {class_query};
        """
    elif fam_filter and not excluded_classes:
        fam_query = ', '.join(f"'{family}'" for family in fam_filter)
        query = f"""
            DELETE FROM TempTable
            WHERE family NOT IN ({fam_query});
        """

    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute(query)
    conn.commit()
    conn.close()
