#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2026
# (c) University of Strathclyde 2026
# (c) James Hutton Institute 2026
#
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
"""Add Pfam domain annotations retrieved from InterPro to a local SQLite database"""


import logging
import sqlite3

from src.sql.interface import insert_data
from src.sql.interface.get_data.get_pfam_data import get_db_pfams
from src.sql.interface.temp_tables import create_temp_pfam_protein_table

logger = logging.getLogger(__name__)


def add_pfam_matches(records: list, conn: sqlite3.Connection) -> tuple[int, int]:
    """Add Pfam family records and their protein matches to the local db.

    :param records: list of PfamMatchRecord objects, each already carrying its target protein_id
    :param conn: open sqlite3 connection

    Return a tuple of (number of new Pfam family rows added, number of new match rows staged)
    """
    db_pfams = get_db_pfams(conn)  # accession: pfam_id
    pfams_to_add = []
    seen_new = set()
    for record in records:
        if record.pfam_accession not in db_pfams and record.pfam_accession not in seen_new:
            pfams_to_add.append((
                record.pfam_accession, record.name, record.annotation_type, record.release
            ))
            seen_new.add(record.pfam_accession)

    if pfams_to_add:
        insert_data(
            conn, "Pfams", ["accession", "name", "annotation_type", "release"], pfams_to_add
        )
        db_pfams = get_db_pfams(conn)

    match_rows = set()
    for record in records:
        pfam_id = db_pfams.get(record.pfam_accession)
        if pfam_id is None:
            logger.warning("Could not resolve db id for Pfam accession %s, skipping match", record.pfam_accession)
            continue
        match_rows.add((
            record.protein_id, pfam_id, record.interpro_accession, record.match_start, record.match_end
        ))

    if match_rows:
        create_temp_pfam_protein_table(conn)
        insert_data(
            conn,
            "TEMP_PFAM_PROTEIN",
            ["protein_id", "pfam_id", "interpro_accession", "match_start", "match_end"],
            list(match_rows),
        )

    conn.commit()

    return len(pfams_to_add), len(match_rows)


def merge_temp_pfam_relationships(db_path: str) -> None:
    """Merge Pfam-protein matches from TEMP_PFAM_PROTEIN into Proteins_Pfams

    Only inserts matches that don't already exist in Proteins_Pfams (same protein, Pfam
    family and match location).

    :param db_path: path to the database
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
        DELETE FROM TEMP_PFAM_PROTEIN
        WHERE rowid NOT IN (
            SELECT MIN(rowid)
            FROM TEMP_PFAM_PROTEIN
            GROUP BY protein_id, pfam_id, match_start, match_end
        );
    """)
    cursor.execute("""
        INSERT OR IGNORE INTO Proteins_Pfams (protein_id, pfam_id, interpro_accession, match_start, match_end)
        SELECT protein_id, pfam_id, interpro_accession, match_start, match_end
        FROM TEMP_PFAM_PROTEIN
    """)

    conn.commit()
    cursor.close()
    conn.close()
