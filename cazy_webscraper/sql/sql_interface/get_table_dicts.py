#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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
# SOFTWARE.
"""Retrieve all objects from a db table and parse the data to build a dict, repr the current table state."""


from tqdm import tqdm

from cazy_webscraper.sql.sql_orm import (
    CazyFamily,
    Ec,
    Genbank,
    Kingdom,
    Pdb,
    Session,
    Taxonomy,
    Uniprot,
)


def get_ec_table_dict(connection):
    """Create dict of objects present in the CazyFamilies table.
    
    :param connection: open sqlalchemy db engine connection
    
    Return dict {ec_number: ec_id}
    """
    with Session(bind=connection) as session:
        db_ec_records = session.query(Ec).all()
    
    ec_table_dict = {}  # {ec_number: ec_id}
    for record in tqdm(db_ec_records, desc="Retrieving existing EC# records"):
        ec_table_dict[record.ec_number] = record.ec_id


def get_ec_gbk_table_dict(connection):
    """Load the Genbanks_Ecs table into memory and compile a dict.

    The the result dict is keyed by EC IDS.
    
    The table contains the current Genbank and EC number relationships in 
    the local CAZyme db.
    
    :param connection: open sqlalchemy connection to an SQLite db
    
    Return dict {ec_id: {gbk ids}}
    """
    with Session(bind=connection) as session:
        all_gbk_ec_records = session.query(Genbank, Ec).\
            join(Ec, Genbank.ecs).\
            all()
        
    ec_gbk_table_dict = {}
    
    for record in all_gbk_ec_records:
        genbank_id = record[0].genbank_id
        ec_id = record[1].ec_id
        
        try:
            ec_gbk_table_dict[ec_id].add(genbank_id)
        except KeyError:
            ec_gbk_table_dict[ec_id] = {genbank_id}
    
    return ec_gbk_table_dict


def get_gbk_ec_table_dict(connection):
    """Load the Genbanks_Ecs table into memory and compile a dict.

    The the result dict is keyed by GENBANK IDS.
    
    The table contains the current Genbank and EC number relationships in 
    the local CAZyme db.
    
    :param connection: open sqlalchemy connection to an SQLite db
    
    Return dict {gbk_id: {ec ids}}
    """
    with Session(bind=connection) as session:
        all_gbk_ec_records = session.query(Genbank, Ec).\
            join(Ec, Genbank.ecs).\
            all()
        
    gbk_ec_table_dict = {}
    
    for record in all_gbk_ec_records:
        genbank_id = record[0].genbank_id
        ec_id = record[1].ec_id
        
        try:
            gbk_ec_table_dict[genbank_id].add(ec_id)
        except KeyError:
            gbk_ec_table_dict[genbank_id] = {ec_id}
    
    return gbk_ec_table_dict


def get_fams_table_dict(connection):
    """Create dict of objects present in the CazyFamilies table.
    
    :param connection: open sqlalchemy db engine connection
    
    Return dict {family subfamily: db_family_id}
    """
    with Session(bind=connection) as session:
        all_families = session.query(CazyFamily).all()
        
    db_fam_dict = {}

    for fam in all_families:
        if fam.subfamily is None:
            subfam = '_'
        else:
            subfam = fam.subfamily
            
        db_fam_dict[f"{fam.family} {subfam}"] = fam.family_id
    
    return db_fam_dict


def get_gbk_table_dict(connection):
    """Compile a dict of the data in the Genbanks table
    
    :param connection: open connection to an SQLite3 database
    
    Return dict {genbank_accession: 'taxa_id': int, 'gbk_id': int}
    """
    with Session(bind=connection) as session:
        all_genbank = session.query(Genbank).all()

    db_gbk_dict = {}  # {genbank_accession: 'taxa_id': str, 'id': int}
    
    for gbk in all_genbank:
        db_gbk_dict[f"{gbk.genbank_accession}"] = {
            'taxa_id': gbk.taxonomy_id,
            'gbk_id': gbk.genbank_id
        }
    
    return db_gbk_dict


def get_gbk_table_seq_dict(connection):
    """Compile a dict of the data in the Genbanks table
    
    :param connection: open connection to an SQLite3 database
    
    Return dict {genbank_accession: 'taxa_id': int, 'gbk_id': int}
    """
    with Session(bind=connection) as session:
        all_genbank = session.query(Genbank).all()

    db_gbk_dict = {}  # {genbank_accession: 'taxa_id': str, 'id': int}
    
    for gbk in all_genbank:
        db_gbk_dict[f"{gbk.genbank_accession}"] = {
            'sequence': gbk.sequence,
            'seq_date': gbk.seq_update_date
        }
    
    return db_gbk_dict


def get_gbk_fam_table_dict(connection):
    """Build dict representing the records present in the Genbanks_CazyFamilies table

    If a GenBank accession is in the db but not has not CazyFamilies instances related to it,
    the GenBank accession is not returned when quering the db.
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return 
    - dict: {gbk_acc: {'families': {'fam subfam': fam_id}}, 'gbk_id': gbk_id }
    - set of tuples: (gbk_id, fam_id), each representing one row in the table
    """
    with Session(bind=connection) as session:
        all_gbk_fam_records = session.query(Genbank, CazyFamily).\
        join(CazyFamily, Genbank.families).\
        all()

    existing_rel_tuples = set()  # set of tuples (gbk_id, fam_id)

    gbk_fam_table_dict = {}
    # {gbk_acc: {'families': {'fam subfam': fam_id}}, 'gbk_id': gbk_id }

    for record in tqdm(all_gbk_fam_records, ' Retreving existing gbk-fam relationships from db'):
        gbk_accession = record[0].genbank_accession
        gbk_id = record[0].genbank_id

        family = record[1].family
        if record[1].subfamily is None:
            subfamily = '_'
        else:
            subfamily = record[1].subfamily
        fam_id = record[1].family_id

        existing_rel_tuples.add( (gbk_id, fam_id) )

        try:
            gbk_fam_table_dict[gbk_accession]

            try:
                gbk_fam_table_dict[gbk_accession][f'{family} {subfamily}']

            except KeyError:
                gbk_fam_table_dict[gbk_accession][f'{family} {subfamily}'] = fam_id

        except KeyError:
            gbk_fam_table_dict[gbk_accession] = {
                'families': {f'{family} {subfamily}': fam_id},
                'gbk_id': gbk_id,
            }

    return gbk_fam_table_dict, existing_rel_tuples


def get_kingdom_table_dict(connection):
    """Load and parse the Kingdoms table from the db and compile a dict {kgnd: id}
    
    :param connection:
    
    Return dict {kingdom: kindom_db_id}
    """
    with Session(bind=connection) as session:
        kingdom_table = session.query(Kingdom).all()
    
    kingdom_dict = {}  # {kingdom: kindom_db_id}
    
    for kingdom_obj in kingdom_table:
        kingdom_dict[kingdom_obj.kingdom] = kingdom_obj.kingdom_id
        
    return kingdom_dict


def get_pdb_table_dict(connection):
    """Create dict of objects present in the Pdbs table.
    
    :param connection: open sqlalchemy db engine connection
    
    Return dict {gbk_id: {pdb accessions}}
    """
    with Session(bind=connection) as session:
        db_pdb_records = session.query(Pdb).all()

    pdb_table_dict = {}  # {gbk_id: {pdb accessions}}

    for record in tqdm(db_pdb_records, desc="Loading existing PDB db records"):
        gbk_id = record.genbank_id
        try:
            pdb_table_dict[gbk_id].add(record.pdb_accession)
            
        except KeyError:
            pdb_table_dict[gbk_id] = {record.pdb_accession}

    pdb_table_dict 


def get_taxs_table_dict(connection):
    """Create dict of objects present in the Taxs table.
    
    :param connection: open sqlalchemy db engine connection
    
    Return dict {genus species: {'tax_id': db_tax_id, 'kingdom_id': kingdom_id}
    """
    with Session(bind=connection) as session:
        all_taxa = session.query(Taxonomy).all()
        
    db_tax_dict = {}
    for taxa in all_taxa:
        if len(taxa.species) == 0:
            db_tax_dict[f"{taxa.genus}"] = {
                'tax_id': taxa.taxonomy_id,
                'kingdom_id': taxa.kingdom_id,
            }
        else:
            db_tax_dict[f"{taxa.genus} {taxa.species}"] = {
                'tax_id': taxa.taxonomy_id,
                'kingdom_id': taxa.kingdom_id,
            }
    
    return db_tax_dict


def get_uniprot_table_dict(connection):
    """Create dict of objects present in the Uniprots table.
    
    :param connection: open sqlalchemy db engine connection
    
    Return dict {acc: {name: str, gbk_id: int, seq: str, seq_date:str } }
    """
    with Session(bind=connection) as session:
        db_uniprot_records = session.query(Uniprot).all()

    uniprot_table_dict = {}  # {acc: {name: str, gbk_id: int, seq: str, seq_date:str } }

    for record in tqdm(db_uniprot_records, desc="Retrieving existing UniProt records from db"):
        uniprot_table_dict[record.uniprot_accession] = {
            "name": record.uniprot_name,
            "genbank_id": record.genbank_id,
            "seq": record.sequence,
            "seq_date": record.seq_update_date,
        }
    
    return uniprot_table_dict


def get_gbk_kingdom_dict(connection):
    """Compile dict of Genbank, Taxonomy and Kingdom records
    
    :param connection: open sqlalchemy db connection
    
    Return dict {kingdom: {genus: {species: {protein_accessions}}}
    """
    with Session(bind=connection) as session:
        query_results = session.query(Genbank, Taxonomy, Kingdom).\
            join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
            join(Genbank, (Genbank.taxonomy_id == Taxonomy.taxonomy_id)).\
            all()

    genbank_kingdom_dict = {}  # kingdom: {genus: {species: {protein_accessions}}}

    for result in tqdm(query_results, desc="Retreving GenBank accessions and taxonomy"):
        genbank_accession = result[0].genbank_accession
        genus = result[1].genus
        species = result[1].species
        kingdom = result[2].kingdom

        try:
            genbank_kingdom_dict[kingdom]

            try:
                genbank_kingdom_dict[kingdom][genus]

                try:
                    genbank_kingdom_dict[kingdom][genus][species].add(genbank_accession)
                
                except KeyError:
                    genbank_kingdom_dict[kingdom][genus][species] = {genbank_accession}

            except KeyError:
                genbank_kingdom_dict[kingdom][genus] = {species: {genbank_accession}}

        except KeyError:
            genbank_kingdom_dict[kingdom] = {
                genus: {
                    species: {genbank_accession},
                },
            }

    return genbank_kingdom_dict
