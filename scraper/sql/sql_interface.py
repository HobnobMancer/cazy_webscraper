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
"""Submodule to interacting with an SQL database"""


import os
import sys

from sqlalchemy import (
    create_engine, Boolean, Column, ForeignKey, Integer, PrimaryKeyConstraint, String, Table
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker


from scraper.sql.sql_orm import (
    Cazyme,
    Taxonomy,
    CazyFamily,
    CazyFamily,
    Cazymes_Genbanks,
    EC,
    Genbank,
    Uniprot,
    Pdb,
    cazymes_families,
    cazymes_ecs,
    cazymes_uniprots,
    cazymes_pdbs,
)


# Use the declarative system
# Database structured in NF1
Base = declarative_base()
Session = sessionmaker()


def build_db(time_stamp, args, logger):
    """Build an empty SQL database and open a session.

    :param time_stamp: str, date and time stamp of when scrape was initated
    :param args: cmd args parser
    :param logger: logger object

    Return an open database session.
    """
    logger.info("Building empty db to store data")

    # build database engine
    if args.database is None:
        if args.output is sys.stdout:
            # write to cwd, this is deleted in scrape is successful
            cwd = os.getcwd()
            db_path = cwd / f"cazy_scrape_temp_{time_stamp}.db"

        else:
            # write to specified output directory
            db_path = args.output / f"cazy_scrape_{time_stamp}.db"

    else:
        # user specificed an existing local CAZy SQL database
        db_path = args.database

    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=False)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)

    return Session()


def add_protein_to_db(
    cazyme_name,
    family,
    source_organism,
    primary_genbank,
    logger,
    session,
    ec_numbers=[],
    genbank_accessions=[],
    uniprot_accessions=[],
    pdb_accessions=[],
):
    """Coordinate adding protein (CAZyme) data to the SQL database (db).

    Every CAZyme will have a cazyme name (called 'protein name' in CAZy), CAZy (sub)family, source
    organism, primary GenBank accession. The primary GenBank accession is the only hyperlinked
    GenBank accession for the protein, and believed to be used by CAZy to indicate the source
    GenBank protein record for the record in CAZy. It can not be guareenteed that a GenBank
    accession will only be recorded as a primary OR a non-primary accession. It may be possible
    that a GenBank accession is the primary accession for one CAZyme and a non-primary accession
    for another. This is believed to be possible becuase CAZy does not appear to ID unique proteins
    by the GenBank accession because duplicate entries for CAZyme can be found within CAZy.

    EC numbers, accessions of associated UniProt and PDB records, non-primary GenBank accessions may
    not be given. For each, if no accessions are collected from CAZy an empty list will be passed.
    For the UniProt and PDB, the first accession listed in each list is recorded as the primary
    accession, and all other listed accessions as non-primary accessions.

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param primary_genbank: str, the hyperlinked GenBank accession from CAZy
    :param logger: logger object
    :param session: open sqlalchemy session to database

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param genbank_accessions: list of non-primary GenBank accessions
    :param uniprot_accessions: list, accessions of associated records in UniProtKB
    :param pdb_accessions: list, accessions of associated records in PDB

    Return nothing or error message.
    """
    # Each unique protein is identified by a unique primary GenBank accession, not its CAZyme name
    # query the local database to see if the current working protein is in the database
    # It is checked that the 'primary_genbank' is logged as a primary accession because if is
    # possible that the accession could be logged as a non-primary accession for another CAZyme
    primary_genbank_query = session.query(Genbank).\
        filter(Genbank.genbank_accession==primary_genbank).\
        filter(Cazymes_Genbanks.primary==True).\
        all()

    if len(primary_genbank_query) == 0:
        # GenBank accession was not found as a primary accession in the database
        # Inferring the protein is not in the local database
        add_new_protein_to_db(
            cazyme_name,
            family,
            source_organism,
            primary_genbank,
            logger,
            session,
            ec_numbers,
            uniprot_accessions,
            pdb_accessions,
        )

    elif len(genbank_query) == 1:
        # GenBank accession found as the primary accession in the database for one protien
        # Identfied that the protein is already catalogued in the local database

        # Get the genbank_id for the GenBank record, then get it's associated Cazyme objects.
        # Above the search filtered for primary==True, so this query only retrieves CAZymes which
        # have the GenBank accession set as their primary accession
        genbank_id = genbank_query[0].genbank_id
        cazyme_query = session.query(Cazyme).join(Cazyme.genbanks).\
            filter(Genbank.genbank_id==genbank_id).all()

        if len(cazyme_query) == 0:
            logger.warning(
                f"The GenBank accession {external_links['GenBank'][0]} with database genbank_id "
                f"{genbank_id}\n"
                "was previously added to the local database, but was not associated with a CAZyme"
            )
            # No protein is associated with the GenBank accession so treat as if it is a new protein
            add_new_protein_to_db(
                cazyme_name,
                family,
                source_organism,
                ec_numbers,
                external_links,
                logger,
                session,
            )

        elif len(cazyme_query) == 1:
            add_data_to_protein_record(
                cazyme_query[0],
                family,
                ec_numbers,
                external_links,
                logger,
                session,
            )

        else:
            # multiple CAZymes are listed with the GenBank accession as their primary GenBank
            # accession, inferring their duplicate records for the same CAZyme
            logger.warning(
                "The following CAZymes found with the same primary GenBank accession "
                f"{external_links['GenBank'][0]}, "
                f"genbank_id={genbank_id},\n"
                "inferring they are the same protein. Therefore, they are potentially duplicate\n"
                "CAZymes in the local database.\n"
                f"Protein data being added to cazyme {cazyme_query[0].cazyme_name}, "
                f"id={cazyme_query[0].cazyme_id}"
            )
            for cazyme in cazyme_query:
                logger.warning(
                    f"Duplicate CAZyme: {cazyme.cazyme_name}, id={cazyme.cazyme_id}"
                )
            add_data_to_protein_record(
                cazyme_query[0],
                family,
                ec_numbers,
                external_links,
                logger,
                session,
            )

    else:
        # The GenBank accession was found multiple times, labeled as the primary accession, in the
        # local database.
        logger.warning(
            "Multiple copies of the GenBank accession "
            f"{external_links['GenBank'][0].genbank_accession} found "
            "in the local database as a primary accession.\n"
            "Checking for potentially duplicate CAZymes in the local database"
        )
        # Check how many CAZymes the GenBank accession is associated with
        duplicate_cazyme_entries = []  # Cazyme objects with the same primary GenBank accession

        for genbank_object in genbank_query:
            # retrieve CAZymes that are linked to the primary GenBank accession
            genbank_id = genbank_object.genbank_id
            cazyme_query = session.query(Cazyme).join(Cazyme.genbanks).\
                filter(Genbank.genbank_id==genbank_id).all()

            if len(cazyme_query) == 0:
                logger.warning(
                    f"The GenBank accession {genbank_object.genbank_accession} with the id "
                    f"{genbank_object.genbank_id}\n"
                    "was previously add to the local database as a primary GenBank accession\n"
                    "but NOT associated with a CAZyme."
                )

            elif len(cazyme_query) == 1:
                duplicate_cazyme_entries.append(cazyme_query[0])

            else:
                # multiple CAZymes have the GenBank accession as their primary GenBank accession
                logger.warning(
                    f"The following CAZymes have the GenBank accession {genbank_object.accession}, "
                    f"genbank_id={genbank_object.genbank_id}\n,"
                    "as their primary GenBank accession,"
                    "infering they are the same CAZyme and are duplicates of each other"
                )
                for cazyme in cazyme_query:
                    logger.warning(
                        f"Duplicate CAZymes: {cazyme.cazyme_name}, cazyme_id={cazyme.cazyme_id}"
                    )
                    duplicate_cazyme_entries.append(cazyme)

        # check how many CAZyme duplicates there are that have the same primary GenBank accession
        if len(duplicate_cazyme_entries) == 0:
            # The GenBank accession has been entered multiple times into the database, but not once
            # associated with a CAZyme.
            logger.warning(
                f"The GenBank accession {external_links['GenBank'][0]} has been entered multiple\n"
                "times as a primary GenBank accession into the local CAZy database with different "
                "'genbank_id's,\n"
                "but it was not associated with a CAZyme, therefore, adding it to the database as "
                "a new CAZyme."
            )
            add_new_protein_to_db(
                cazyme_name,
                family,
                source_organism,
                ec_numbers,
                external_links,
                logger,
                session,
            )

        elif len(duplicate_cazyme_entries) == 1:
            # only one CAZyme found associated with the primary GenBank accession
            # add the protein data to this CAZyme
            logger.warning(
                f"The GenBank accession {external_links['GenBank'][0]} has been entered multiple\n"
                "times as a primary GenBank accession into the local CAZy database with different "
                "'genbank_id's,\n"
                f"but associated with only one CAZyme: {duplicate_cazyme_entries[0].cazyme_name} "
                f"id={duplicate_cazyme_entries[0].cazyme_id}"
            )
            add_data_to_protein_record(
                duplicate_cazyme_entries[0],
                family,
                ec_numbers,
                external_links,
                logger,
                session,
            )

        else:
            logger.warning(
                "The following CAZymes appear to be duplicates in the local database becuase\n"
                "they have the same primary GenBank accession number, "
                f"{external_links['GenBank'][0]}.\n"
                "Although the GenBank accessions have different genbank_ids, therefore, the\n"
                "accession has been added to the database as a primary accession multiple times.\n"
                f"Protein data being added to {duplicate_cazyme_entries[0].cazyme_name} "
                f"id={duplicate_cazyme_entries[0].cazyme_id}"
            )
            for cazyme in duplicate_cazyme_entries:
                logger.warning(
                    f"Duplicate CAZyme: {cazyme.cazyme_name}, id={cazyme.cazyme_id}"
                )
            # add the protein data scraped from CAZy to the first duplicate CAZyme
            add_data_to_protein_record(
                duplicate_cazyme_entries[0],
                family,
                ec_numbers,
                external_links,
                logger,
                session,
            )

    return


def add_new_protein_to_db(
    cazyme_name,
    family,
    source_organism,
    primary_genbank,
    logger,
    session,
    ec_numbers=[],
    genbank_accessions=[],
    uniprot_accessions=[],
    pdb_accessions=[],
):
    """Add a new protein (a CAZyme) record to the database.

    EC numbers, accessions of associated UniProt and PDB records, non-primary GenBank accessions may
    not be given. For each, if no accessions are collected from CAZy an empty list will be passed.
    For the UniProt and PDB, the first accession listed in each list is recorded as the primary
    accession, and all other listed accessions as non-primary accessions.

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param primary_genbank: str, the hyperlinked GenBank accession from CAZy
    :param logger: logger object
    :param session: open sqlalchemy session to database

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param genbank_accessions: list of non-primary GenBank accessions
    :param uniprot_accessions: list, accessions of associated records in UniProtKB
    :param pdb_accessions: list, accessions of associated records in PDB

    Return nothing.
    """
    # define the CAZyme
    new_cazyme = Cazyme(cazyme_name=cazyme_name)

    # define source organism
    genus_sp_separator = source_organism.find(" ")
    genus = source_organism[:genus_sp_separator]
    species = source_organism[genus_sp_separator:]

    # check if the source  organism is already catalogued
    tax_query = session.query(Taxonomy).filter_by(genus=genus, species=species).all()

    if len(tax_query) == 0:
        # create Taxonomy model object
        new_cazyme.taxonomy = Taxonomy(genus=genus, species=species)

    elif len(tax_query) == 1:
        # add exiting Taxonomy object to the new cazyme
        new_cazyme.taxonomy = tax_query[0]

    else:
        logger.warning(
            f"The species {genus} {species} has been loaded into the database "
            f"{len(tax_query)} times.\n"
        )
        new_cazyme.taxonomy = tax_query[0]

    # add the new cazyme to the database
    session.add(new_cazyme)
    session.commit()

    # add Family/Subfamily classifications
    if family.find("_") != -1:
        add_cazy_subfamily(family, new_cazyme, session, logger)
    else:
        add_cazy_family(family, new_cazyme, session, logger)

    # define ec_number
    if len(ec_numbers) != 0:
        add_ec_numbers(ec_numbers, new_cazyme, session, logger)

    # add non-primary GenBank accessions
    if len(genbank_accessions) != 0:
        add_genbank_accessions(genbank_accessions, new_cazyme, session, logger)

    # add UniProt accessions
    if len(genbank_accessions) != 0:
        add_uniprot_accessions(uniprot_accessions, new_cazyme, session, logger)

    # add PDB/3D accessions
    if len(pdb_accessions) != 0:
        add_pdb_accessions(pdb_accessions, new_cazyme, session, logger)

    # final commit to ensure all changes are commited
    session.commit()

    return


def add_data_to_protein_record(
    cazyme_name,
    family,
    source_organism,
    primary_genbank,
    logger,
    session,
    ec_numbers=[],
    genbank_accessions=[],
    uniprot_accessions=[],
    pdb_accessions=[],
):
    """Add data to an existing record in the SQL database.

    EC numbers, accessions of associated UniProt and PDB records, non-primary GenBank accessions may
    not be given. For each, if no accessions are collected from CAZy an empty list will be passed.
    For the UniProt and PDB, the first accession listed in each list is recorded as the primary
    accession, and all other listed accessions as non-primary accessions.

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param primary_genbank: str, the hyperlinked GenBank accession from CAZy
    :param logger: logger object
    :param session: open sqlalchemy session to database

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param genbank_accessions: list of non-primary GenBank accessions
    :param uniprot_accessions: list, accessions of associated records in UniProtKB
    :param pdb_accessions: list, accessions of associated records in PDB

    Return nothing.
    """
    # add Family/Subfamily classifications
    if family.find("_") != -1:
        add_cazy_subfamily(family, cazyme, session, logger)
    else:
        add_cazy_family(family, cazyme, session, logger)

    # define ec_number
    if len(ec_numbers) != 0:
        add_ec_numbers(ec_numbers, new_cazyme, session, logger)

    # add non-primary GenBank accessions
    if len(genbank_accessions) != 0:
        add_genbank_accessions(genbank_accessions, new_cazyme, session, logger)

    # add UniProt accessions
    if len(genbank_accessions) != 0:
        add_uniprot_accessions(uniprot_accessions, new_cazyme, session, logger)

    # add PDB/3D accessions
    if len(pdb_accessions) != 0:
        add_pdb_accessions(pdb_accessions, new_cazyme, session, logger)

    # final commit to ensure all changes are commited
    session.commit()

    return


def add_cazy_family(family, cazyme, session, logger):
    """Add CAZy family to CAZyme record in the local CAZy database.

    :param family: str, name of a CAZy family/subfamily
    :param cazyme: Cazymes class object
    :param session: open local database session connector
    :param logger: logger object

    Return nothing.
    """
    # Check if Family is already in the local database
    family_query = session.query(CazyFamily).filter_by(family=family).all()

    if len(family_query) == 0:
        # add new CAZy to the database
        cazy_family = CazyFamily(family=family)
        session.add(cazy_family)
        session.commit()

        if not(cazy_family in cazyme.families):
            cazyme.families.append(cazy_family)
        session.commit()

    elif len(family_query) == 1:
        # check if it is associated with a CAZy subfamily

        if family_query[0].subfamily is not None:
            # The current family record is associated with a subfamily
            # the current working family is not therefore
            # add new CAZy family without a subfamily association to the database
            cazy_family = CazyFamily(family=family)
            session.add(cazy_family)
            session.commit()

            if not(family_query[0] in cazyme.families):
                cazyme.families.append(cazy_family)
            session.commit()

        else:
            # add existing Family record to current working CAZyme
            if not(family_query[0] in cazyme.families):
                cazyme.families.append(family_query[0])
            session.commit()

    else:
        # check if any of the multiple Family records are associated with subfamilies or not
        nonsubfamily_asociated_family_entries = []

        for entry in family_query:
            # check if retrieved family record is associated with a subfamily
            # retrieve only entries NOT associated with a subfamily
            if entry.subfamily is None:
                nonsubfamily_asociated_family_entries.append(entry)

        if len(nonsubfamily_asociated_family_entries) == 0:
            # add new CAZy family without a subfamily association to the database
            cazy_family = CazyFamily(family=family)
            session.add(cazy_family)
            session.commit()
            if not(cazy_family in cazyme.families):
                cazyme.families.append(cazy_family)
            session.commit()

        elif len(nonsubfamily_asociated_family_entries) == 1:
            # Add existing non-subfamily associated family instance to CAZyme
            if not(nonsubfamily_asociated_family_entries[0] in cazyme.families):
                cazyme.families.append(nonsubfamily_asociated_family_entries[0])
            session.commit()

        else:
            logger.warning(
                f"Duplicate family entries found for the CAZy family {family} "
                "without subfamily association.\n"
                "Linking the family with the id "
                f"{nonsubfamily_asociated_family_entries[0].family_id} to the cazyme "
                f"{cazyme.cazyme_name} id={cazyme.cazyme_id}"
            )
            if not(nonsubfamily_asociated_family_entries[0] in cazyme.families):
                cazyme.families.append(nonsubfamily_asociated_family_entries[0])
            session.commit()

    session.commit()
    return


def add_cazy_subfamily(family, cazyme, session, logger):
    """Add CAZy family to CAZyme record in the local CAZy database.

    :param family: str, name of a CAZy family/subfamily
    :param cazyme: Cazymes class object
    :param session: open local database session connector
    :param logger: logger object

    Return nothing.
    """
    # check if the subfamily is already in the database
    subfam_query = session.query(CazyFamily).filter_by(subfamily=family).all()

    if len(subfam_query) == 0:
        # add new CAZy subfamily to the database
        cazy_family = CazyFamily(family=family[:family.find("_")], subfamily=family)
        session.add(cazy_family)
        session.commit()

        if not(cazy_family in cazyme.families):
            cazyme.families.append(cazy_family)
        session.commit()

    elif len(subfam_query) == 1:
        # add existing subfamily to the current working CAZyme
        if not(subfam_query[0] in cazyme.families):
            cazyme.families.append(subfam_query[0])
        session.commit()

    else:
        logger.warning(
            f"Duplicate subfamily entries found for the CAZy subfamily {family}.\n"
            f"Add the family with the id {subfam_query[0].family_id} to the cazyme "
            f"{cazyme.cazyme_name} id={cazyme.cazyme_id}"
        )
        if not(subfam_query[0] in cazyme.families):
            cazyme.families.append(subfam_query[0])
        session.commit()

    session.commit()
    return


def add_ec_numbers(ec_numbers, cazyme, session, logger):
    """Add EC numbers to CAZyme record in the local CAZy database.

    :param ec_numbers: list of EC numbers (str)
    :param new_cazyme: Cazymes class object
    :param session: open local database session connector
    :param logger: logger object

    Return nothing.
    """
    for ec in ec_numbers:
        # check if the EC number is already in the database
        ec_query = session.query(EC).filter_by(ec_number=ec).all()

        if len(ec_query) == 0:
            # add new EC number to the database
            new_ec = EC(ec_number=ec)
            session.add(new_ec)
            session.commit()

            if not(new_ec in cazyme.ecs):
                cazyme.ecs.append(new_ec)
            session.commit()

        elif len(ec_query) == 1:
            # add existing EC number record to the CAZyme record
            if not(ec_query[0] in cazyme.ecs):
                cazyme.ecs.append(ec_query[0])
            session.commit()

        else:
            # duplicate records found for the current EC number
            logger.warning(
                f"Duplicate entries found for the EC# {ec} in the local CAZy database."
                f"Add the family with the id {ec_query[0].ec_id} to the cazyme "
                f"{cazyme.cazyme_name} id={cazyme.cazyme_id}"
            )
            if not(ec_query[0] in cazyme.ecs):
                cazyme.ecs.append(ec_query[0])
            session.commit()

    session.commit()
    return


def add_genbank_accessions(genbank_accessions, cazyme, session, logger):
    """Add GenBank protein accession numbers to CAZyme record in the local CAZy database.

    :param genbank_accessions: list of GenBank protein accession numbers (str)
    :param cazyme: Cazymes class object
    :param session: open local database session connector
    :param logger: logger object

    Return nothing.
    """
    if len(genbank_accessions) == 0:
        return

    elif len(genbank_accessions) == 1:
        add_primary_genbank(genbank_accessions[0], cazyme, session, logger)
        # check if accession is already in the database

    else:
        add_primary_genbank(genbank_accessions[0], cazyme, session, logger)

        for accession in genbank_accessions[1:]:
            # check if accession is in the database already
            genbank_query = session.query(Genbank).filter_by(genbank_accession=accession).all()

            if len(genbank_query) == 0:
                # add new genbank accession
                new_accession = Genbank(genbank_accession=accession, primary=False)
                session.add(new_accession)
                session.commit()

                if not(new_accession in cazyme.genbanks):
                    cazyme.genbanks.append(new_accession)
                session.commit()

            elif len(genbank_query) == 1:
                # add GenBank record to current working CAZyme
                if not(genbank_query[0] in cazyme.genbanks):
                    cazyme.genbanks.append(genbank_query[0])
                session.commit()

            else:
                logger.warning(
                    f"Duplicate entries for GenBank accession {accession}"
                    f"Adding the accession entry with the ID {genbank_query[0].genbank_id} to the\n"
                    f"cazyme {cazyme.cazyme_name} id={cazyme.cazyme_id}"
                )
                if not(genbank_query[0] in cazyme.genbanks):
                    cazyme.genbanks.append(genbank_query[0])
                session.commit()

    session.commit()
    return


def add_primary_genbank(accession, cazyme, session, logger):
    """Add primary GenBank accession for CAZyme to the local CAZy database.

    :param accession: str, primary GenBank accession of a CAZyme
    :param cazyme: Cazymes class object
    :param session: open local database session connector
    :param logger: logger object

    Return nothing.
    """
    genbank_query = session.query(Genbank).filter_by(genbank_accession=accession).all()

    if len(genbank_query) == 0:
        # add new genbank accession
        new_accession = Genbank(genbank_accession=accession, primary=True)
        session.add(new_accession)
        session.commit()

        if not(new_accession in cazyme.genbanks):
            cazyme.genbanks.append(new_accession)
        session.commit()

    elif len(genbank_query) == 1:
        # add GenBank record to current working CAZyme
        if not(genbank_query[0] in cazyme.genbanks):
            cazyme.genbanks.append(genbank_query[0])
        session.commit()

    else:
        logger.warning(
            f"Duplicate entries for GenBank accession {accession}"
            f"Adding the accession entry with the ID {genbank_query[0].genbank_id} to the\n"
            f"cazyme {cazyme.cazyme_name} id={cazyme.cazyme_id}"
        )
        if not(genbank_query[0] in cazyme.genbanks):
            cazyme.genbanks.append(genbank_query[0])
        session.commit()

    session.commit()
    return


def add_uniprot_accessions(uniprot_accessions, cazyme, session, logger):
    """Add UniProt protein accessions to CAZyme record in the local CAZy database.

    :param uniprot_accessions: list of UniProt protein accession numbers (str)
    :param cazyme: Cazymes class object
    :param session: open local database session connector
    :param logger: logger object

    Return nothing.
    """
    if len(uniprot_accessions) == 0:
        return

    elif len(uniprot_accessions) == 1:
        add_primary_uniprot(uniprot_accessions[0], cazyme, session, logger)
        # check if accession is already in the database

    else:
        add_primary_uniprot(uniprot_accessions[0], cazyme, session, logger)

        for accession in uniprot_accessions[1:]:
            # check if accession is in the database already
            uniprot_query = session.query(Uniprot).filter_by(uniprot_accession=accession).all()

            if len(uniprot_query) == 0:
                # add new uniprot accession
                new_accession = Uniprot(uniprot_accession=accession, primary=False)
                session.add(new_accession)
                session.commit()

                if not(new_accession in cazyme.uniprots):
                    cazyme.uniprots.append(new_accession)
                session.commit()

            elif len(uniprot_query) == 1:
                # add UniProt record to current working CAZyme
                if not(uniprot_query[0] in cazyme.uniprots):
                    cazyme.uniprots.append(uniprot_query[0])
                session.commit()

            else:
                logger.warning(
                    f"Duplicate entries for UniProt accession {accession}"
                    f"Adding the accession entry with the ID {uniprot_query[0].uniprot_id} to the\n"
                    f"cazyme {cazyme.cazyme_name} id={cazyme.cazyme_id}"
                )
                if not(uniprot_query[0] in cazyme.uniprots):
                    cazyme.uniprots.append(uniprot_query[0])
                session.commit()

    session.commit()
    return


def add_primary_uniprot(accession, cazyme, session, logger):
    """Add primary UniProt accession for CAZyme to the local CAZy database.

    :param accession: str, primary UniProt accession of a CAZyme
    :param cazyme: Cazymes class object
    :param session: open local database session connector
    :param logger: logger object

    Return nothing.
    """
    uniprot_query = session.query(Uniprot).filter_by(uniprot_accession=accession).all()

    if len(uniprot_query) == 0:
        # add new uniprot accession
        new_accession = Uniprot(uniprot_accession=accession, primary=True)
        session.add(new_accession)
        session.commit()

        if not(new_accession in cazyme.uniprots):
            cazyme.uniprots.append(new_accession)
        session.commit()

    elif len(uniprot_query) == 1:
        # add UniProt record to current working CAZyme
        if not(uniprot_query[0] in cazyme.uniprots):
            cazyme.uniprots.append(uniprot_query[0])
        session.commit()

    else:
        logger.warning(
            f"Duplicate entries for UniProt accession {accession}"
            f"Adding the accession entry with the ID {uniprot_query[0].uniprot_id} to the\n"
            f"cazyme {cazyme.cazyme_name} id={cazyme.cazyme_id}"
        )
        if not(uniprot_query[0] in cazyme.uniprots):
            cazyme.uniprots.append(uniprot_query[0])
        session.commit()

    session.commit()
    return


def add_pdb_accessions(pdb_accessions, cazyme, session, logger):
    """Add PDB/3D protein accessions to CAZyme record in the local CAZy database.

    :param pdb_accessions: list of UniProt protein accession numbers (str)
    :param cazyme: Cazymes class object
    :param session: open local database session connector
    :param logger: logger object

    Return nothing.
    """
    if len(pdb_accessions) == 0:
        return

    elif len(pdb_accessions) == 1:
        add_primary_pdb(pdb_accessions[0], cazyme, session, logger)
        # check if accession is already in the database

    else:
        add_primary_pdb(pdb_accessions[0], cazyme, session, logger)

        for accession in pdb_accessions[1:]:
            # check if accession is in the database already
            pdb_query = session.query(Pdb).filter_by(pdb_accession=accession).all()

            if len(pdb_query) == 0:
                # add new pdb accession
                new_accession = Pdb(pdb_accession=accession, primary=False)
                session.add(new_accession)
                session.commit()

            if not(new_accession in cazyme.pdbs):
                cazyme.pdbs.append(new_accession)
                session.commit()

            elif len(pdb_query) == 1:
                # add PDB/3D record to current working CAZyme
                if not(pdb_query[0] in cazyme.pdbs):
                    cazyme.pdbs.append(pdb_query[0])
                session.commit()

            else:
                logger.warning(
                    f"Duplicate entries for PDB/3D accession {accession}"
                    f"Adding the accession entry with the ID {pdb_query[0].pdb_id} to the\n"
                    f"cazyme {cazyme.cazyme_name} id={cazyme.cazyme_id}"
                )
                if not(pdb_query[0] in cazyme.pdbs):
                    cazyme.pdbs.append(pdb_query[0])
                session.commit()

    session.commit()
    return


def add_primary_pdb(accession, cazyme, session, logger):
    """Add primary PDB/3D accession for CAZyme to the local CAZy database.

    :param accession: str, primary UniProt accession of a CAZyme
    :param cazyme: Cazymes class object
    :param session: open local database session connector
    :param logger: logger object

    Return nothing.
    """
    pdb_query = session.query(Pdb).filter_by(pdb_accession=accession).all()

    if len(pdb_query) == 0:
        # add new pdb accession
        new_accession = Pdb(pdb_accession=accession, primary=True)
        session.add(new_accession)
        session.commit()

        if not(new_accession in cazyme.pdbs):
            cazyme.pdbs.append(new_accession)
        session.commit()

    elif len(pdb_query) == 1:
        # add PDB/3D record to current working CAZyme
        if not(pdb_query[0] in cazyme.pdbs):
            cazyme.pdbs.append(pdb_query[0])
        session.commit()

    else:
        logger.warning(
            f"Duplicate entries for PDB/3D accession {accession}"
            f"Adding the accession entry with the ID {pdb_query[0].pdb_id} to the\n"
            f"cazyme {cazyme.cazyme_name} id={cazyme.cazyme_id}"
        )
        if not(pdb_query[0] in cazyme.pdbs):
            cazyme.pdbs.append(pdb_query[0])
        session.commit()

    session.commit()
    return
