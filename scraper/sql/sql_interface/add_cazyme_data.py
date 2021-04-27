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
"""Create and add to CAZyme records in a local SQL database"""


import logging

from sqlalchemy.exc import IntegrityError
from scraper.sql.sql_orm import (
    Cazyme,
    Taxonomy,
    CazyFamily,
    Cazymes_Genbanks,
    EC,
    Genbank,
    Kingdom,
    Uniprot,
    Pdb,
)


def add_new_protein_to_db(
    cazyme_name,
    family,
    source_organism,
    tax_kingdom,
    primary_genbank_object,
    session,
    ec_numbers=[],
    gbk_nonprimary=[],
    uni_primary=[],
    uni_nonprimary=[],
    pdb_accessions=[],
):
    """Add a new protein (a CAZyme) record to the database.

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param tax_kingdom: str, taxonomy Kingdom of the source organism
    :param primary_genbank_object: sql_orm.Genbank instance, represents the primary accession
    :param session: open sqlalchemy session to database

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param gbk_nonprimary: list of non-primary GenBank accessions
    :param uni_primary: list, primary accessions of associated records in UniProtKB
    :param uni_nonprimary: list, non-primary accessions of associated records in UniProtKB
    :param pdb_accessions: list, =accessions of associated records in PDB

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    error_message = None

    # add the taxonomy Kingdom
    try:
        new_kingdom = Kingdom(kingdom=f'{tax_kingdom[0].upper()}{tax_kingdom[1:]}')
        session.add(new_kingdom)
        session.commit()
    except IntegrityError:  # raised when Kingdom is already in the Kingdom table
        session.rollback()
        # retrieve the existing Kingdom instance
        tax_kingdom = f'{tax_kingdom[0].upper()}{tax_kingdom[1:]}'
        new_kingdom = session.query(Kingdom).filter(Kingdom.kingdom == tax_kingdom).all()[0]

    # define source organism
    if source_organism == 'unidentified':
        genus = source_organism
        species = source_organism

    else:
        genus_sp_separator = source_organism.find(" ")
        genus = source_organism[:genus_sp_separator]
        species = source_organism[genus_sp_separator:]

    # Add the CAZyme and its taxonomy data
    try:
        new_tax = Taxonomy(genus=genus, species=species)
        new_tax.tax_kingdom = new_kingdom

        new_cazyme = Cazyme(cazyme_name=cazyme_name)
        new_cazyme.taxonomy = new_tax

        session.add(new_cazyme)
        session.commit()

    except IntegrityError:  # raised if Taxonomy record already in the db
        session.rollback()

        try:
            new_cazyme = Cazyme(cazyme_name=cazyme_name)
            tax_query = session.query(Taxonomy).\
                filter(Taxonomy.genus == genus).filter(Taxonomy.species == species).all()
            tax_query[0].cazymes.append(new_cazyme)
            session.commit()

        except IntegrityError:
            logger.warning(
                f"Could not add Taxonomy data for cazyme {cazyme_name} from {genus} {species}"
            )
            error_message = (
                f"Could not add CAZyme {cazyme_name} from {source_organism} to the database"
            )

    # establish relationship between the CAZyme and its primary GenBank accession
    try:
        relationship = Cazymes_Genbanks(
            cazymes=new_cazyme,
            genbanks=primary_genbank_object,
            primary=True,
        )
        session.add(relationship)
        session.commit()
    except IntegrityError:  # raised if CAZyme and GenBank accession are already linked
        session.rollback()

    if family.find("_") != -1:
        add_cazy_subfamily(family, new_cazyme, session)
    else:
        add_cazy_family(family, new_cazyme, session)

    if len(ec_numbers) != 0:
        add_ec_numbers(ec_numbers, new_cazyme, session)

    if len(gbk_nonprimary) != 0:
        add_nonprimary_gbk_accessions(gbk_nonprimary, new_cazyme, session)

    if len(uni_primary) != 0:
        add_uniprot_accessions(uni_primary, new_cazyme, True, session)

    if len(uni_nonprimary) != 0:
        add_uniprot_accessions(uni_nonprimary, new_cazyme, False, session)

    if len(pdb_accessions) != 0:
        add_pdb_accessions(pdb_accessions, new_cazyme, session)

    return error_message


def add_data_to_protein_record(
    cazyme,
    family,
    session,
    ec_numbers=[],
    gbk_nonprimary=[],
    uni_primary=[],
    uni_nonprimary=[],
    pdb_accessions=[],
):
    """Add data to an existing record in the SQL database.

    :param cazyme_name: Cazyme class instance - class imported from scraper.sql.sql_orm
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param session: open sqlalchemy session to database

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param gbk_nonprimary: list of non-primary GenBank accessions
    :param uni_primary: list, primary accessions of associated records in UniProtKB
    :param uni_nonprimary: list, non-primary accessions of associated records in UniProtKB
    :param pdb_accessions: list, accessions of associated records in PDB

    Return nothing.
    """
    if family.find("_") != -1:
        add_cazy_subfamily(family, cazyme, session)
    else:
        add_cazy_family(family, cazyme, session)

    if len(ec_numbers) != 0:
        add_ec_numbers(ec_numbers, cazyme, session)

    if len(gbk_nonprimary) != 0:
        add_nonprimary_gbk_accessions(gbk_nonprimary, cazyme, session)

    if len(uni_primary) != 0:
        add_uniprot_accessions(uni_primary, cazyme, True, session)

    if len(uni_nonprimary) != 0:
        add_uniprot_accessions(uni_nonprimary, cazyme, False, session)

    if len(pdb_accessions) != 0:
        add_pdb_accessions(pdb_accessions, cazyme, session)

    return


def add_cazy_family(family, cazyme, session):
    """Add CAZy family to CAZyme record in the local CAZy database.

    :param family: str, name of a CAZy family/subfamily
    :param cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    try:
        cazyme.families.append(CazyFamily(family=family))
        session.commit()
        return

    except IntegrityError:
        session.rollback()

    # get the record that caused raising the Integrity Error
    family_query = session.query(CazyFamily).filter(CazyFamily.family == family).all()

    try:
        cazyme.families.append(family_query[0])
        session.commit()
    except (IntegrityError, IndexError) as e:
        session.rollback()

    return


def add_cazy_subfamily(subfamily, cazyme, session):
    """Add CAZy family to CAZyme record in the local CAZy database.

    :param family: str, name of a CAZy family/subfamily
    :param cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    family = subfamily[:subfamily.find("_")]

    try:
        cazyme.families.append(CazyFamily(family=family, subfamily=subfamily))
        session.commit()
        return

    except IntegrityError:
        session.rollback()

    # get the record that caused raising the IntegrityError
    subfam_query = session.query(CazyFamily).\
        filter(CazyFamily.family == family).filter(CazyFamily.subfamily == subfamily).all()

    try:
        cazyme.families.append(subfam_query[0])
        session.commit()
    except (IntegrityError, IndexError) as e:
        session.rollback()

    return


# Note: Not all CAZymes have non-primary GenBank accessions, EC numbers, UniProt accessions or
# PDB accessions


def add_nonprimary_gbk_accessions(genbank_accessions, cazyme, session):
    """Add non-primary GenBank protein accessions to the local database.

    :param genbank_accessions: list of non-primary GenBank protein accession numbers (str)
    :param cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    for accession in genbank_accessions:
        try:
            new_genbank = Genbank(genbank_accession=accession)
            session.add(new_genbank)
            session.commit()

        except IntegrityError:
            session.rollback()

            genbank_query = session.query(Genbank, Cazymes_Genbanks).\
                filter(Genbank.genbank_accession == accession).\
                filter(Cazymes_Genbanks.primary == False).all()

            try:
                relationship = Cazymes_Genbanks(
                    cazyme_id=cazyme.cazyme_id,
                    genbank_id=genbank_query[0][0].genbank_id,
                    primary=False,
                )
                session.add(relationship)
                session.commit()
            except (IntegrityError, IndexError) as e:
                session.rollback()
            continue

        try:
            relationship = Cazymes_Genbanks(cazymes=cazyme, genbanks=new_genbank, primary=False)
            session.add(relationship)
            session.commit()
        except IntegrityError:
            session.rollback()

    return


def add_ec_numbers(ec_numbers, cazyme, session):
    """Add EC numbers to CAZyme record in the local CAZy database.

    :param ec_numbers: list of EC numbers (str)
    :param new_cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    for ec in ec_numbers:
        try:
            new_ec = EC(ec_number=ec)
            session.add(new_ec)
            session.commit()

        except IntegrityError:
            session.rollback()

            ec_query = session.query(EC).filter(EC.ec_number == ec).all()
            try:
                cazyme.ecs.append(ec_query[0])
                session.commit()
            except (IntegrityError, IndexError) as e:
                session.rollback()
            continue

        try:
            cazyme.ecs.append(new_ec)
            session.commit()
        except IntegrityError:
            session.rollback()

    return


def add_uniprot_accessions(uniprot_accessions, cazyme, primary, session):
    """Add primary UniProt accessions to CAZyme record in local CAZy database.

    :param uniprot_accessions: list of Uniprot accessions
    :param cazyme: Cazyme class instance
    :param primary: boolean, True = primary accessions, False non-primary accessions
    :param session: open local database session connector

    Return nothing.
    """
    for accession in uniprot_accessions:
        try:
            new_uniprot = Uniprot(uniprot_accession=accession, primary=primary)
            session.add(new_uniprot)
            session.commit()
        except IntegrityError:
            session.rollback()

            uniprot_query = session.query(Uniprot).filter(Uniprot.uniprot_accession == accession).\
                filter(Uniprot.primary == primary).all()

            try:
                cazyme.uniprots.append(uniprot_query[0])
                session.commit()
            except (IntegrityError, IntegrityError) as e:
                session.rollback()
            continue

        try:
            cazyme.uniprots.append(new_uniprot)
            session.commit()
        except IntegrityError:
            session.rollback()

    return


def add_pdb_accessions(pdb_accessions, cazyme, session):
    """Add PDB/3D protein accessions to CAZyme record in the local CAZy database.

    :param pdb_accessions: list of UniProt protein accession numbers (str)
    :param cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    for accession in pdb_accessions:
        try:
            new_pdb = Pdb(pdb_accession=accession)
            session.add(new_pdb)
            session.commit()
        except IntegrityError:
            session.rollback()

            pdb_query = session.query(Pdb).filter(Pdb.pdb_accession == accession).all()

            try:
                cazyme.pdbs.append(pdb_query[0])
                session.commit()
            except (IntegrityError, IntegrityError) as e:
                session.rollback()
            continue

        try:
            cazyme.pdbs.append(new_pdb)
            session.commit()
        except IntegrityError:
            session.rollback()

    return
