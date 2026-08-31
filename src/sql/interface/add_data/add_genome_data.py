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
"""Add genomic genome information to the local CAZyme database."""


import logging
import sqlite3

from argparse import Namespace


logger = logging.getLogger(__name__)


def persist_genome_data(args: Namespace, batch_size: int = 1000) -> dict:
    """Persist genome data from temporary tables to main database tables.

    Args:
        db_path: Path to the local SQLite database
        update: Whether to update existing records
        batch_size: Number of records to process in each batch

    Returns:
        Dictionary with statistics about the operation
    """
    conn = sqlite3.connect(args.database)
    cursor = conn.cursor()

    stats = {
        'genomes_added': 0,
        'genomes_updated': 0,
        'protein_relationships_added': 0
    }

    try:
        # Step 1: Get count of temp genomes for batching
        cursor.execute("SELECT COUNT(*) FROM TEMP_GENOME")
        total_genomes = cursor.fetchone()[0]

        # Step 2: Process genomes in batches
        offset = 0
        while offset < total_genomes:
            batch_stats = _process_genome_batch(conn, offset, batch_size, args.update)

            # Update overall stats
            for key in stats:
                stats[key] += batch_stats[key]

            offset += batch_size

            # Commit after each batch to avoid large transactions
            conn.commit()

            logger.debug("Processed batch %d-%d: added %d, updated %d", 
                        offset - batch_size + 1, min(offset, total_genomes),
                        batch_stats['genomes_added'], batch_stats['genomes_updated'])

        # Step 3: Update TEMP_GENOME2PROTEIN with genome_ids in batches
        _update_temp_protein_genome_ids(conn, batch_size)

        # Step 4: Process protein-genome relationships in batches
        relationship_stats = _process_protein_relationships(conn, batch_size)
        stats['protein_relationships_added'] = relationship_stats['added']

        conn.commit()

        # Log final statistics
        logger.warning("Genome data persistence completed:")
        logger.warning("  - Genomes added: %d", stats['genomes_added'])
        logger.warning("  - Genomes updated: %d", stats['genomes_updated'])
        logger.warning("  - Protein relationships added: %d", stats['protein_relationships_added'])

    except Exception as e:
        conn.rollback()
        logger.error("Error persisting genome data: %s", e)
        raise
    finally:
        conn.close()

    return stats


def _process_genome_batch(conn: sqlite3.Connection, offset: int, batch_size: int, update: bool) -> dict:
    """Process a batch of genomes from TEMP_GENOME table.

    Each temp row is one NCBI assembly carrying both its GenBank (GCA_) and RefSeq (GCF_)
    accession, so a single Genomes row is written per assembly with both accession columns
    populated, rather than one row per accession flavour.

    Args:
        conn: Database connection
        offset: Starting offset for the batch
        batch_size: Size of the batch
        update: Whether to update existing records

    Returns:
        Dictionary with batch statistics
    """
    cursor = conn.cursor()
    batch_stats = {'genomes_added': 0, 'genomes_updated': 0, 'protein_relationships_added': 0}

    cursor.execute("""
        SELECT ncbi_genome_id, gbk_accession, refseq_accession, assembly_name
        FROM TEMP_GENOME
        LIMIT ? OFFSET ?
    """, (batch_size, offset))
    temp_genomes = cursor.fetchall()

    if not temp_genomes:
        return batch_stats

    # look up any assembly we already hold, by either of its accessions
    accessions = [row[1] for row in temp_genomes if row[1]] + [row[2] for row in temp_genomes if row[2]]
    existing_by_accession = {}
    if accessions:
        placeholders = ','.join('?' * len(accessions))
        cursor.execute(f"""
            SELECT genome_id, assembly_name, gbk_version_accession, refseq_version_accession
            FROM Genomes
            WHERE gbk_version_accession IN ({placeholders})
               OR refseq_version_accession IN ({placeholders})
        """, accessions + accessions)

        for genome_id, assembly_name, gbk_acc, refseq_acc in cursor.fetchall():
            record = {
                'genome_id': genome_id,
                'assembly_name': assembly_name,
                'gbk_version_accession': gbk_acc,
                'refseq_version_accession': refseq_acc,
            }
            if gbk_acc:
                existing_by_accession[gbk_acc] = record
            if refseq_acc:
                existing_by_accession[refseq_acc] = record

    genomes_to_add = []
    genomes_to_update = []

    for ncbi_genome_id, gbk_accession, refseq_accession, assembly_name in temp_genomes:
        existing = existing_by_accession.get(gbk_accession) or existing_by_accession.get(refseq_accession)

        if existing is None:
            genomes_to_add.append((
                assembly_name,
                gbk_accession,
                ncbi_genome_id if gbk_accession else None,
                refseq_accession,
                ncbi_genome_id if refseq_accession else None,
            ))
        elif update and (
            existing['assembly_name'] != assembly_name
            or existing['gbk_version_accession'] != gbk_accession
            or existing['refseq_version_accession'] != refseq_accession
        ):
            genomes_to_update.append((
                assembly_name, gbk_accession, refseq_accession, existing['genome_id']
            ))

    if genomes_to_add:
        cursor.executemany("""
            INSERT INTO Genomes (assembly_name, gbk_version_accession, gbk_ncbi_id,
                                 refseq_version_accession, refseq_ncbi_id)
            VALUES (?, ?, ?, ?, ?)
        """, genomes_to_add)
        batch_stats['genomes_added'] = len(genomes_to_add)

    if genomes_to_update:
        cursor.executemany("""
            UPDATE Genomes
            SET assembly_name = ?, gbk_version_accession = ?, refseq_version_accession = ?
            WHERE genome_id = ?
        """, genomes_to_update)
        batch_stats['genomes_updated'] = len(genomes_to_update)

    return batch_stats


def _update_temp_protein_genome_ids(conn: sqlite3.Connection, batch_size: int) -> None:
    """Update TEMP_GENOME2PROTEIN with genome_ids in batches.

    Args:
        conn: Database connection
        batch_size: Size of batches to process
    """
    cursor = conn.cursor()

    # Get count of records to update
    cursor.execute("SELECT COUNT(*) FROM TEMP_GENOME2PROTEIN WHERE genome_id IS NULL")
    total_to_update = cursor.fetchone()[0]

    if not total_to_update:
        logger.info("All TEMP_GENOME2PROTEIN records already have genome_ids")
        return

    logger.info("Updating %d protein-genome mappings with genome_ids in batches of %d", 
               total_to_update, batch_size)

    # Process in batches using ROWID for consistent ordering
    offset = 0
    updated_count = 0

    while offset < total_to_update:
        # Update a batch of records
        cursor.execute("""
            UPDATE TEMP_GENOME2PROTEIN
            SET genome_id = (
                SELECT g.genome_id
                FROM TEMP_GENOME tg
                JOIN Genomes g ON (
                    (tg.gbk_accession IS NOT NULL AND g.gbk_version_accession = tg.gbk_accession)
                    OR (tg.refseq_accession IS NOT NULL AND g.refseq_version_accession = tg.refseq_accession)
                )
                WHERE tg.ncbi_genome_id = TEMP_GENOME2PROTEIN.ncbi_genome_id
            )
            WHERE ROWID IN (
                SELECT ROWID FROM TEMP_GENOME2PROTEIN 
                WHERE genome_id IS NULL 
                LIMIT ? OFFSET ?
            )
        """, (batch_size, offset))

        batch_updated = cursor.rowcount
        updated_count += batch_updated
        offset += batch_size

        # Commit after each batch
        conn.commit()

        logger.debug("Updated batch: %d records (total: %d/%d)", 
                    batch_updated, updated_count, total_to_update)

        if batch_updated == 0:
            break  # No more records to update

    logger.info("Updated %d protein-genome mappings with genome_ids", updated_count)


def _process_protein_relationships(conn: sqlite3.Connection, batch_size: int) -> dict:
    """Process protein-genome relationships in batches.

    Args:
        conn: Database connection
        batch_size: Size of batches to process

    Returns:
        Dictionary with relationship statistics
    """
    cursor = conn.cursor()
    stats = {'added': 0}

    # Get count of temp relationships
    cursor.execute("SELECT COUNT(*) FROM TEMP_GENOME2PROTEIN WHERE genome_id IS NOT NULL")
    total_relationships = cursor.fetchone()[0]

    if not total_relationships:
        logger.info("No protein-genome relationships to process")
        return stats

    logger.info("Processing %d protein-genome relationships in batches of %d", 
               total_relationships, batch_size)

    # Get existing relationships in batches to avoid memory issues
    existing_relationships = set()
    cursor.execute("SELECT COUNT(*) FROM Proteins_Genomes")
    existing_count = cursor.fetchone()[0]

    if existing_count > 0:
        logger.info("Loading %d existing protein-genome relationships", existing_count)
        offset = 0
        while offset < existing_count:
            cursor.execute("SELECT protein_id, genome_id FROM Proteins_Genomes LIMIT ? OFFSET ?", 
                          (batch_size, offset))
            batch_existing = cursor.fetchall()
            existing_relationships.update(batch_existing)
            offset += batch_size

            if len(batch_existing) < batch_size:
                break

    # Process temp relationships in batches
    offset = 0
    while offset < total_relationships:
        cursor.execute("""
            SELECT DISTINCT protein_id, genome_id 
            FROM TEMP_GENOME2PROTEIN 
            WHERE genome_id IS NOT NULL
            LIMIT ? OFFSET ?
        """, (batch_size, offset))

        temp_batch = cursor.fetchall()
        if not temp_batch:
            break

        # Filter out existing relationships
        relationships_to_add = []
        for protein_id, genome_id in temp_batch:
            if (protein_id, genome_id) not in existing_relationships:
                relationships_to_add.append((protein_id, genome_id))
                existing_relationships.add((protein_id, genome_id))  # Avoid duplicates in subsequent batches

        # Add new relationships
        if relationships_to_add:
            cursor.executemany("""
                INSERT INTO Proteins_Genomes (protein_id, genome_id)
                VALUES (?, ?)
            """, relationships_to_add)
            stats['added'] += len(relationships_to_add)

            logger.debug("Added %d new relationships in batch %d-%d", 
                        len(relationships_to_add), offset + 1, offset + len(temp_batch))

        offset += batch_size

        # Commit after each batch
        conn.commit()

    logger.info("Added %d new protein-genome relationships", stats['added'])
    return stats
