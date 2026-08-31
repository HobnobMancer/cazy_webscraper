#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# Author:
# Emma E. M. Hobbs

# Contact
# ehobbbs@ebi.ac.uk

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
"""Get genome genome data from NCBI Assembly database"""

import gzip
import logging
import re
import urllib.error
import requests
import time
import sqlite3

"""Persist PDB entry metadata from the temp table into the local CAZyme database."""


import logging
import sqlite3

from argparse import Namespace

from src.sql.interface.connect import get_sqlite3_connection
from src.sql.interface.temp_tables import drop_temp_pdb_structure_table


logger = logging.getLogger(__name__)


def persist_pdb_data(args: Namespace) -> dict[str, int]:
    """Merge TEMP_PDB_STRUCTURE into the Pdbs table.

    PDB accessions themselves are added by the UniProt pipeline (they come from UniProt's
    cross-references); this fills in or refreshes the experimental method and resolution from
    RCSB, which is the authoritative source for them. Rows are only updated when --update is
    used or when the existing value is missing, so a plain re-run does not overwrite data.

    Returns a dict of statistics about the operation.
    """
    conn = get_sqlite3_connection(args.database)
    cursor = conn.cursor()

    stats = {"metadata added": 0, "metadata updated": 0, "not in db": 0}

    try:
        # accessions retrieved from PDB that are not in the local db at all: nothing to
        # attach them to, since the Proteins<->Pdbs relationship comes from UniProt
        cursor.execute("""
            SELECT COUNT(*)
            FROM TEMP_PDB_STRUCTURE T
            LEFT JOIN Pdbs P ON T.pdb_accession = P.pdb_accession
            WHERE P.pdb_id IS NULL
        """)
        stats["not in db"] = cursor.fetchone()[0]

        # fill in rows that have no method/resolution yet
        cursor.execute("""
            SELECT COUNT(*)
            FROM TEMP_PDB_STRUCTURE T
            JOIN Pdbs P ON T.pdb_accession = P.pdb_accession
            WHERE (P.method IS NULL AND T.method IS NOT NULL)
               OR (P.resolution IS NULL AND T.resolution IS NOT NULL)
        """)
        stats["metadata added"] = cursor.fetchone()[0]

        cursor.execute("""
            UPDATE Pdbs
            SET method = COALESCE(
                    (SELECT T.method FROM TEMP_PDB_STRUCTURE T
                     WHERE T.pdb_accession = Pdbs.pdb_accession),
                    Pdbs.method
                ),
                resolution = COALESCE(
                    (SELECT T.resolution FROM TEMP_PDB_STRUCTURE T
                     WHERE T.pdb_accession = Pdbs.pdb_accession),
                    Pdbs.resolution
                )
            WHERE (Pdbs.method IS NULL OR Pdbs.resolution IS NULL)
              AND Pdbs.pdb_accession IN (SELECT pdb_accession FROM TEMP_PDB_STRUCTURE)
        """)

        if args.update:
            # refresh rows whose stored value disagrees with what RCSB now reports
            cursor.execute("""
                SELECT COUNT(*)
                FROM TEMP_PDB_STRUCTURE T
                JOIN Pdbs P ON T.pdb_accession = P.pdb_accession
                WHERE (T.method IS NOT NULL AND IFNULL(P.method, '') != T.method)
                   OR (T.resolution IS NOT NULL AND IFNULL(P.resolution, -1) != T.resolution)
            """)
            stats["metadata updated"] = cursor.fetchone()[0]

            cursor.execute("""
                UPDATE Pdbs
                SET method = COALESCE(
                        (SELECT T.method FROM TEMP_PDB_STRUCTURE T
                         WHERE T.pdb_accession = Pdbs.pdb_accession),
                        Pdbs.method
                    ),
                    resolution = COALESCE(
                        (SELECT T.resolution FROM TEMP_PDB_STRUCTURE T
                         WHERE T.pdb_accession = Pdbs.pdb_accession),
                        Pdbs.resolution
                    )
                WHERE Pdbs.pdb_accession IN (SELECT pdb_accession FROM TEMP_PDB_STRUCTURE)
            """)

        conn.commit()

        logger.warning("PDB metadata persistence completed:")
        logger.warning("  - Entries with new method/resolution: %d", stats["metadata added"])
        if args.update:
            logger.warning("  - Entries refreshed: %d", stats["metadata updated"])
        if stats["not in db"]:
            logger.warning(
                "  - Retrieved entries not present in the local db: %d", stats["not in db"]
            )

    except sqlite3.Error as err:
        conn.rollback()
        logger.error("Error persisting PDB metadata: %s", err)

    finally:
        cursor.close()
        conn.close()

    return stats


def rm_temp_pdb_table(db_path) -> None:
    """Drop the PDB metadata temp table once its contents have been merged."""
    conn = get_sqlite3_connection(db_path)
    drop_temp_pdb_structure_table(conn)
    conn.close()
