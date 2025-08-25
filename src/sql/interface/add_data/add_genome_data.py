#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Add genomic assembly information to the local CAZyme database."""

import logging
from tqdm import tqdm
from typing import Dict
from src.sql.interface.get_data.get_assemblies import (
    get_assembly_table,
    get_gbk_genome_table_data,
)
from src.databases.ncbi.genomes import NcbiAssembly

logger = logging.getLogger(__name__)

def update_assembly_data(assembly_prot_dict: Dict[str, list], ncbi_assembly_dict: Dict[str, NcbiAssembly], gbk_dict: Dict[str, int], connection, args):
    """Add/update genomic data and link proteins with genomes in the local database.

    Args:
        assembly_prot_dict: {assembly_name: [protein_accs]}
        ncbi_assembly_dict: {assembly_name: NcbiAssembly}
        gbk_dict: {protein_acc: genbank_id}
        connection: Database connection
        args: Command line arguments
    """
    # Get assembly names from the NcbiAssembly objects
    assembly_names = list(ncbi_assembly_dict.keys())
    db_genome_table_dict = get_assembly_table(assembly_names, connection)
    genomes_to_add, genomes_to_update = [], []

    for assembly_name in assembly_names:
        (genomes_to_update if assembly_name in db_genome_table_dict else genomes_to_add).append(assembly_name)

    genomes_of_interest = genomes_to_add + genomes_to_update

    if args.update and genomes_to_update:
        update_genomic_data(genomes_to_update, db_genome_table_dict, ncbi_assembly_dict, connection)
        db_genome_table_dict = get_assembly_table(genomes_of_interest, connection)

    if genomes_to_add:
        add_genomic_data(ncbi_assembly_dict, genomes_to_add, connection)
        db_genome_table_dict = get_assembly_table(genomes_of_interest, connection)

    add_prot_genome_relationships(assembly_prot_dict, gbk_dict, db_genome_table_dict, connection)

def add_genomic_data(ncbi_assembly_dict: Dict[str, NcbiAssembly], genomes_of_interest: list, connection):
    """Add new genomic assembly data to the database.

    Args:
        ncbi_assembly_dict: Dictionary mapping assembly names to NcbiAssembly objects
        genomes_of_interest: List of assembly names to add
        connection: Database connection
    """
    records = []
    for name in genomes_of_interest:
        if name in ncbi_assembly_dict:
            assembly = ncbi_assembly_dict[name]
            records.append((
                name,
                assembly.genbank_acc,
                assembly.assembly_id,
                assembly.refseq_acc,
                assembly.assembly_id,  # Using assembly_id for both since we may not have separate UIDs
            ))
    
    if records:
        try:
            with connection:
                connection.executemany(
                    """INSERT INTO Genomes (
                        assembly_name, gbk_version_accession, gbk_ncbi_id,
                        refseq_version_accession, refseq_ncbi_id
                    ) VALUES (?, ?, ?, ?, ?)""",
                    records
                )
                logger.info(f"Added {len(records)} new genome assemblies to database")
        except Exception as e:
            logger.error(f"Failed to insert genomic data: {e}")

def update_genomic_data(genomes_of_interest: list, genome_table_dict: Dict[str, int], ncbi_assembly_dict: Dict[str, NcbiAssembly], connection, unit_test=False):
    """Update existing genomic assembly data in the database.
    
    Args:
        genomes_of_interest: List of assembly names to update
        genome_table_dict: Dictionary mapping assembly names to database IDs
        ncbi_assembly_dict: Dictionary mapping assembly names to NcbiAssembly objects
        connection: Database connection
        unit_test: Whether this is a unit test (for rollback)
    """
    for name in genomes_of_interest:
        if name not in ncbi_assembly_dict:
            continue
            
        db_id = genome_table_dict[name]
        assembly = ncbi_assembly_dict[name]
        
        vals = (
            assembly.genbank_acc,
            assembly.assembly_id,
            assembly.refseq_acc,
            assembly.assembly_id,  # Using assembly_id for both since we may not have separate UIDs
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

def add_prot_genome_relationships(assembly_prot_dict: Dict[str, list], gbk_dict: Dict[str, int], db_genome_table_dict: Dict[str, int], connection):
    """Link proteins and genomes in the Genbanks_Genomes table.
    
    Args:
        assembly_prot_dict: Dictionary mapping assembly names to lists of protein accessions
        gbk_dict: Dictionary mapping protein accessions to GenBank IDs
        db_genome_table_dict: Dictionary mapping assembly names to genome IDs
        connection: Database connection
    """
    existing = get_gbk_genome_table_data(connection)
    to_add = set()
    
    for assembly_name, protein_accs in assembly_prot_dict.items():
        if assembly_name not in db_genome_table_dict:
            logger.warning(f"Assembly {assembly_name} not found in database")
            continue
            
        genome_id = db_genome_table_dict[assembly_name]
        
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
