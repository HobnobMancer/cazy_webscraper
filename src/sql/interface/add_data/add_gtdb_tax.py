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
"""Add GTDB Taxonomy data to a local SQLite database"""


import logging
import sqlite3

from argparse import Namespace

from src.sql.interface.connect import get_sqlite3_connection
from src.sql.interface.temp_tables import drop_temp_gtdb_table


logger = logging.getLogger(__name__)


def persist_gtdb_data(args: Namespace) -> dict[str, int]:
    """Merge TEMP_GTDB into the GtdbTaxs table and link the genomes to their lineages.

    GtdbTaxs has a unique constraint over the seven ranks plus the release, so distinct
    lineages are inserted once and shared by every genome that resolves to them.

    Returns a dict of statistics about the operation.
    """
    conn = get_sqlite3_connection(args.database)
    cursor = conn.cursor()

    stats = {"lineages added": 0, "genomes linked": 0, "genomes relinked": 0}

    try:
        cursor.execute("SELECT COUNT(*) FROM GtdbTaxs")
        before = cursor.fetchone()[0]

        # insert each distinct lineage once; the unique constraint makes this idempotent
        cursor.execute("""
            INSERT OR IGNORE INTO GtdbTaxs
                (kingdom, phylum, tax_class, tax_order, family, genus, species, release)
            SELECT DISTINCT kingdom, phylum, tax_class, tax_order, family, genus, species, release
            FROM TEMP_GTDB
        """)

        cursor.execute("SELECT COUNT(*) FROM GtdbTaxs")
        stats["lineages added"] = cursor.fetchone()[0] - before

        # genomes that would have their existing lineage replaced by a different one
        cursor.execute("""
            SELECT COUNT(*)
            FROM Genomes G
            JOIN TEMP_GTDB T ON G.genome_id = T.genome_id
            JOIN GtdbTaxs GT ON
                 IFNULL(GT.kingdom, '')   = IFNULL(T.kingdom, '')
             AND IFNULL(GT.phylum, '')    = IFNULL(T.phylum, '')
             AND IFNULL(GT.tax_class, '') = IFNULL(T.tax_class, '')
             AND IFNULL(GT.tax_order, '') = IFNULL(T.tax_order, '')
             AND IFNULL(GT.family, '')    = IFNULL(T.family, '')
             AND IFNULL(GT.genus, '')     = IFNULL(T.genus, '')
             AND IFNULL(GT.species, '')   = IFNULL(T.species, '')
             AND IFNULL(GT.release, '')   = IFNULL(T.release, '')
            WHERE G.gtdb_tax_id IS NOT NULL AND G.gtdb_tax_id != GT.gtdb_tax_id
        """)
        stats["genomes relinked"] = cursor.fetchone()[0]

        # only overwrite an existing link when --update was given
        where_clause = (
            "" if args.update else " AND G.gtdb_tax_id IS NULL"
        )

        cursor.execute(f"""
            UPDATE Genomes AS G
            SET gtdb_tax_id = (
                SELECT GT.gtdb_tax_id
                FROM TEMP_GTDB T
                JOIN GtdbTaxs GT ON
                     IFNULL(GT.kingdom, '')   = IFNULL(T.kingdom, '')
                 AND IFNULL(GT.phylum, '')    = IFNULL(T.phylum, '')
                 AND IFNULL(GT.tax_class, '') = IFNULL(T.tax_class, '')
                 AND IFNULL(GT.tax_order, '') = IFNULL(T.tax_order, '')
                 AND IFNULL(GT.family, '')    = IFNULL(T.family, '')
                 AND IFNULL(GT.genus, '')     = IFNULL(T.genus, '')
                 AND IFNULL(GT.species, '')   = IFNULL(T.species, '')
                 AND IFNULL(GT.release, '')   = IFNULL(T.release, '')
                WHERE T.genome_id = G.genome_id
            )
            WHERE G.genome_id IN (SELECT genome_id FROM TEMP_GTDB){where_clause}
        """)
        stats["genomes linked"] = cursor.rowcount

        conn.commit()

        logger.warning("GTDB data persistence completed:")
        logger.warning("  - New GTDB lineages added: %d", stats["lineages added"])
        logger.warning("  - Genomes linked to a GTDB lineage: %d", stats["genomes linked"])
        if stats["genomes relinked"] and not args.update:
            logger.warning(
                "  - Genomes with a different lineage already stored (use --update to "
                "overwrite): %d",
                stats["genomes relinked"],
            )

    except sqlite3.Error as err:
        conn.rollback()
        logger.error("Error persisting GTDB data: %s", err)

    finally:
        cursor.close()
        conn.close()

    return stats


def rm_temp_gtdb_table(db_path) -> None:
    """Drop the GTDB temp table once its contents have been merged."""
    conn = get_sqlite3_connection(db_path)
    drop_temp_gtdb_table(conn)
    conn.close()
