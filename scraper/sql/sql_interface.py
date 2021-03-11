#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License
"""Submodule for interacting with a local SQL database"""


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


class SqlInterfaceException(Exception):
    """General exception for SQL interface"""

    def __init__(self, message):
        self.message = message


def add_protein_to_db(
    cazyme_name,
    family,
    source_organism,
    kingdom,
    primary_genbank,
    session,
    ec_numbers=[],
    gbk_nonprimary=[],
    uni_primary=[],
    uni_nonprimary=[],
    pdb_accessions=[],
):
    """Coordinate adding protein (CAZyme) data to the SQL database (db).

    Every CAZyme has a cazyme name (called 'protein name' in CAZy), CAZy (sub)family, source
    organism, and primary GenBank accession. The 'protein_name' assigned by CAZy is not unique for
    each protein. The primary GenBank accession is the only hyperlinked GenBank accession for the
    protein, and believed to be used by CAZy to indicate the source GenBank protein record. It may
    be possible that a GenBank accession is the primary accession for one CAZyme and a non-primary
    accession for another, becuase CAZy does not appear to ID unique proteins  by the GenBank
    accession because duplicate CAZyme entries can be found within CAZy.

    EC numbers, UniProt accessions and PDB accessions may not be given. The UniProt and PDB
    accessions are recorded as the primary accessions, and all other listed accessions as
    non-primary accessions.

    There are a minority of CAZymes in CAZy with NO GenBank accessions. These CAZymes are annotated
    with the GenBank accession (str) 'NA' in the local db, and are logged.

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param kingdom: str, taxonomy Kingdom of the source organism
    :param primary_genbank: str, the hyperlinked GenBank accession from CAZy
    :param session: open sqlalchemy session to database

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param gbk_nonprimary: list of non-primary GenBank accessions
    :param uni_primary: list, primary accessions of associated records in UniProtKB
    :param uni_nonprimary: list, non-primary accessions of associated records in UniProtKB
    :param pdb_accessions: list, accessions of associated records in PDB

    Return nothing.
    """
    error_message = None
    # Each unique protein is identified by a unique primary GenBank accession
    try:
        primary_genbank_object = Genbank(genbank_accession=primary_genbank)
        session.add(primary_genbank_object)
        session.commit()

    except IntegrityError:
        # raised when Genbank accession already in the database
        session.rollback()  # enable continued interation with the database

        parse_unique_genbank_conflict(
            cazyme_name,
            family,
            source_organism,
            kingdom,
            primary_genbank,
            session,
            ec_numbers,
            gbk_nonprimary,
            uni_primary,
            uni_nonprimary,
            pdb_accessions,
        )

        return

    error_message = add_new_protein_to_db(
        cazyme_name,
        family,
        source_organism,
        kingdom,
        primary_genbank_object,
        session,
        ec_numbers,
        gbk_nonprimary,
        uni_primary,
        uni_nonprimary,
        pdb_accessions,
    )

    if error_message is not None:
        raise SqlInterfaceException(error_message)

    return


def parse_unique_genbank_conflict(
    cazyme_name,
    family,
    source_organism,
    kingdom,
    primary_genbank,
    session,
    ec_numbers=[],
    gbk_nonprimary=[],
    uni_primary=[],
    uni_nonprimary=[],
    pdb_accessions=[],
):
    """Called when primary GenBank accession is already in the local database.

    Check if the accession in the database is a primary accession. If it is then add data to the
    existing protein record. If not then add the protein data as a new protein to the database.

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param kingdom: str, taxonomy Kingdom of the source_organism
    :param primary_genbank: str, the hyperlinked GenBank accession from CAZy
    :param session: open sqlalchemy session to database

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param gbk_nonprimary: list of non-primary GenBank accessions
    :param uni_primary: list, primary accessions of associated records in UniProtKB
    :param uni_nonprimary: list, non-primary accessions of associated records in UniProtKB
    :param pdb_acccessions: list, accessions of associated records in PDB

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    error_message = None

    # check if the local GenBank object is for a primary GenBank accession
    primary_genbank_query = session.query(Genbank).\
        filter(Genbank.genbank_accession == primary_genbank).\
        filter(Cazymes_Genbanks.primary == True).all()

    if len(primary_genbank_query) == 0:
        # Accession was not found as a primary accession, inferring the protein is not in the db
        error_message = add_new_protein_to_db(
            cazyme_name,
            family,
            source_organism,
            kingdom,
            primary_genbank_query[0],
            session,
            ec_numbers,
            gbk_nonprimary,
            uni_primary,
            uni_nonprimary,
            pdb_accessions,
        )

    else:
        # Primary GenBank accession found, add additional data to the respective CAZyme
        cazyme_query = session.query(Cazyme, Genbank, Cazymes_Genbanks).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            join(Cazyme, (Cazyme.cazyme_id == Cazymes_Genbanks.cazyme_id)).\
            filter(Genbank.genbank_accession == primary_genbank).\
            filter(Cazymes_Genbanks.primary == True).all()

        if len(cazyme_query) == 0:
            logger.warning(
                f"GenBank accession {primary_genbank}, id={primary_genbank_query[0].genbank_id}\n"
                "was previously added to the local database, but not associated with a CAZyme\n"
                "Adding protein data as a new CAZyme in the local database."
            )
            error_message = add_new_protein_to_db(
                cazyme_name,
                family,
                source_organism,
                kingdom,
                primary_genbank_query[0],
                session,
                ec_numbers,
                gbk_nonprimary,
                uni_primary,
                uni_nonprimary,
                pdb_accessions,
            )

        elif len(cazyme_query) == 1:
            add_data_to_protein_record(
                cazyme_query[0][0],
                family,
                session,
                ec_numbers,
                gbk_nonprimary,
                uni_primary,
                uni_nonprimary,
                pdb_accessions,
            )

        else:
            # multiple CAZymes have the same primary GenBank accession, inferring duplicates
            logger.warning(
                "Potential duplicate CAZymes found in the local database, because they have the "
                f"same primary GenBank accession {primary_genbank}\nPotential duplicate CAZymes:"
            )
            for cazyme in cazyme_query:
                logger.warning(f"Name={cazyme[0].cazyme_name}, id={cazyme[0].cazyme_id}")
            logger.warning(
                f"Protein data added to cazyme {cazyme_query[0][0].cazyme_name}, "
                f"id={cazyme_query[0][0].cazyme_id}"
            )
            add_data_to_protein_record(
                cazyme_query[0][0],
                family,
                session,
                ec_numbers,
                gbk_nonprimary,
                uni_primary,
                uni_nonprimary,
                pdb_accessions,
            )

    return error_message


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
