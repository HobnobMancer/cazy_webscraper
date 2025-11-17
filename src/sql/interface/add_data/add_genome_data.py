#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Add genomic genome information to the local CAZyme database."""

import logging
from src.sql.interface.get_data.get_genomes import (
    get_genome_table,
    get_gbk_genome_table_data,
)
from src.databases.ncbi.genomes import NcbiGenome

logger = logging.getLogger(__name__)

def update_genome_data(genome_prot_dict: dict, ncbi_genome_dict: dict, gbk_dict: dict, args):
    """Add/update genomic data and link proteins with genomes in the local database.

    Args:
        genome_prot_dict: {genome_name: [protein_accs]}
        ncbi_genome_dict: {genome_name: NcbiGenome}
        gbk_dict: {protein_acc: genbank_id}
        connection: Database connection
        args: Command line arguments
    """
    # Get genome names from the NcbiGenome objects
    genome_names = list(ncbi_genome_dict.keys())
    db_genome_table_dict = get_genome_table(genome_names, connection)
    genomes_to_add, genomes_to_update = [], []

    for genome_name in genome_names:
        (genomes_to_update if genome_name in db_genome_table_dict else genomes_to_add).append(genome_name)

    genomes_of_interest = genomes_to_add + genomes_to_update

    if args.update and genomes_to_update:
        update_genomic_data(genomes_to_update, db_genome_table_dict, ncbi_genome_dict, connection)
        db_genome_table_dict = get_genome_table(genomes_of_interest, connection)

    if genomes_to_add:
        add_genomic_data(ncbi_genome_dict, genomes_to_add, connection)
        db_genome_table_dict = get_genome_table(genomes_of_interest, connection)

    add_prot_genome_relationships(genome_prot_dict, gbk_dict, db_genome_table_dict, connection)

def add_genomic_data(ncbi_genome_dict: dict, genomes_of_interest: list, connection):
    """Add new genomic genome data to the database.

    Args:
        ncbi_genome_dict: Dictionary mapping genome names to NcbiGenome objects
        genomes_of_interest: List of genome names to add
        connection: Database connection
    """
    records = []
    for name in genomes_of_interest:
        if name in ncbi_genome_dict:
            genome = ncbi_genome_dict[name]
            records.append((
                name,
                genome.genbank_acc,
                genome.genome_id,
                genome.refseq_acc,
                genome.genome_id,  # Using genome_id for both since we may not have separate UIDs
            ))
    
    if records:
        try:
            with connection:
                connection.executemany(
                    """INSERT INTO Genomes (
                        genome_name, gbk_version_accession, gbk_ncbi_id,
                        refseq_version_accession, refseq_ncbi_id
                    ) VALUES (?, ?, ?, ?, ?)""",
                    records
                )
                logger.info(f"Added {len(records)} new genome assemblies to database")
        except Exception as e:
            logger.error(f"Failed to insert genomic data: {e}")

def update_genomic_data(genomes_of_interest: list, genome_table_dict: dict, ncbi_genome_dict: dict, connection, unit_test=False):
    """Update existing genomic genome data in the database.
    
    Args:
        genomes_of_interest: List of genome names to update
        genome_table_dict: Dictionary mapping genome names to database IDs
        ncbi_genome_dict: Dictionary mapping genome names to NcbiGenome objects
        connection: Database connection
        unit_test: Whether this is a unit test (for rollback)
    """
    for name in genomes_of_interest:
        if name not in ncbi_genome_dict:
            continue
            
        db_id = genome_table_dict[name]
        genome = ncbi_genome_dict[name]
        
        vals = (
            genome.genbank_acc,
            genome.genome_id,
            genome.refseq_acc,
            genome.genome_id,  # Using genome_id for both since we may not have separate UIDs
            db_id,
        )
        
        try:
            with connection:
                connection.execute(
                    """UPDATE Genomes SET
                        gbk_version_accession = ?, gbk_ncbi_id = ?,
                        refseq_version_accession = ?, refseq_ncbi_id = ?
                        WHERE genome_id = ?""",
                    vals
                )
                if unit_test:
                    connection.rollback()
        except Exception as e:
            logger.error(f"Failed to update genome {name}: {e}")

def add_prot_genome_relationships(genome_prot_dict: dict[str, list], gbk_dict: dict[str, int], db_genome_table_dict: dict[str, int], connection):
    """Link proteins and genomes in the Genbanks_Genomes table.
    
    Args:
        genome_prot_dict: Dictionary mapping genome names to lists of protein accessions
        gbk_dict: Dictionary mapping protein accessions to GenBank IDs
        db_genome_table_dict: Dictionary mapping genome names to genome IDs
        connection: Database connection
    """
    existing = get_gbk_genome_table_data(connection)
    to_add = set()
    
    for genome_name, protein_accs in genome_prot_dict.items():
        if genome_name not in db_genome_table_dict:
            logger.warning(f"Assembly {genome_name} not found in database")
            continue
            
        genome_id = db_genome_table_dict[genome_name]
        
        for protein_acc in protein_accs:
            if protein_acc not in gbk_dict:
                logger.warning(f"Protein {protein_acc} not found in GenBank dictionary")
                continue
                
            genbank_id = gbk_dict[protein_acc]
            relationship = (genbank_id, genome_id)
            
            if relationship not in existing:
                to_add.add(relationship)
    
    if to_add:
        try:
            with connection:
                connection.executemany(
                    "INSERT INTO Genbanks_Genomes (genbank_id, genome_id) VALUES (?, ?)",
                    list(to_add)
                )
                logger.info(f"Added {len(to_add)} protein-genome relationships to database")
        except Exception as e:
            logger.error(f"Failed to insert protein-genome relationships: {e}")
