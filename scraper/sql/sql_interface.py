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
"""Submodule to interacting with a local SQL database"""


import logging

from scraper.sql.sql_orm import (
    Cazyme,
    Taxonomy,
    CazyFamily,
    Cazymes_Genbanks,
    EC,
    Genbank,
    Uniprot,
    Pdb,
)


def add_protein_to_db(
    cazyme_name,
    family,
    source_organism,
    primary_genbank,
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

    There are a minority of CAZymes in CAZy for which NO GenBank accessions are listed.
    These CAZymes will be annotated with the GenBank accession (str) 'NA', and it will be
    logged which CAZymes this is for.

    These proteins are added as new CAZymes to the database becuase the only true method of
    checking for unique or identical proteins is by the GenBank accession. This is becuase
    CAZy retrieves the protein sequences from GenBank. The 'protein_name' assigned by CAZy is
    not unique for each protein, multiple proteins with different GenBank accessions can be
    assigned the same 'protein_name', and aer different proteins. Therefore, there is no way
    to check if a protein has already been catalogued without a GenBank accession. Consequently,
    when a protein is found with no GenBank accession it is added as a new CAZyme with the 
    GenBank accession 'NA'.

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param primary_genbank: str, the hyperlinked GenBank accession from CAZy
    :param session: open sqlalchemy session to database

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param genbank_accessions: list of non-primary GenBank accessions
    :param uniprot_accessions: list, accessions of associated records in UniProtKB
    :param pdb_accessions: list, accessions of associated records in PDB

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # Each unique protein is identified by a unique primary GenBank accession, not its CAZyme name
    # query the local database to see if the current working protein is in the database
    # It is checked that the 'primary_genbank' is logged as a primary accession because it is
    # possible that the accession could be logged as a non-primary accession for another CAZyme
    primary_genbank_query = session.query(Genbank).\
        filter(Genbank.genbank_accession==primary_genbank).\
        filter(Cazymes_Genbanks.primary==True).all()

    if len(primary_genbank_query) == 0:
        # GenBank accession was not found as a primary accession in the database
        # Inferring the protein is not in the local database
        add_new_protein_to_db(
            cazyme_name,
            family,
            source_organism,
            primary_genbank,
            session,
            ec_numbers,
            genbank_accessions,
            uniprot_accessions,
            pdb_accessions,
        )

    elif len(primary_genbank_query) == 1:
        # GenBank accession found as the primary accession for one CAZyme, thus additional data
        # can be added to the CAZyme record

        # retrieve the CAZyme record
        cazyme_query = session.query(Cazyme, Genbank, Cazymes_Genbanks).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            join(Cazyme, (Cazyme.cazyme_id == Cazymes_Genbanks.cazyme_id)).\
            filter(Genbank.genbank_accession == primary_genbank).\
            filter(Cazymes_Genbanks.primary == True).\
            all()

        if len(cazyme_query) == 0:
            logger.warning(
                f"The GenBank accession {primary_genbank} with database genbank_id "
                f"{primary_genbank_query[0].genbank_id}\n"
                "was previously added to the local database, but was not associated with a CAZyme\n"
                "Adding protein data as a new CAZyme in the local database."
            )
            # No protein is associated with the GenBank accession so treat as if it is a new protein
            add_new_protein_to_db(
                cazyme_name,
                family,
                source_organism,
                primary_genbank,
                session,
                ec_numbers,
                genbank_accessions,
                uniprot_accessions,
                pdb_accessions,
            )

        elif len(cazyme_query) == 1:
            add_data_to_protein_record(
                cazyme_query[0][0],
                family,
                source_organism,
                session,
                ec_numbers,
                genbank_accessions,
                uniprot_accessions,
                pdb_accessions,
            )

        else:
            # multiple CAZymes are listed with the same primary GenBank accession, inferring they
            # are duplicate records for the same CAZyme
            logger.warning(
                "Potential duplicate CAZymes found in the local database.\n"
                "The following CAZymes found with the same primary GenBank accession "
                f"{primary_genbank}\n inferring they are the same protein"
            )
            for cazyme in cazyme_query:
                logger.warning(
                    f"Duplicate CAZyme: {cazyme[0].cazyme_name}, id={cazyme[0].cazyme_id}"
                )
            logger.warning(
                f"Protein data added to cazyme {cazyme_query[0][0].cazyme_name}, "
                f"id={cazyme_query[0][0].cazyme_id}"
            )

            add_data_to_protein_record(
                cazyme_query[0][0],
                family,
                source_organism,
                session,
                ec_numbers,
                genbank_accessions,
                uniprot_accessions,
                pdb_accessions,
            )

    else:
        # The primary GenBank accession was found multiple times
        logger.warning(
            f"Multiple copies of the primary GenBank accession {primary_genbank} found.\n"
            "Checking for potentially duplicate CAZymes in the local database"
        )
        # Check if there are duplicate CAZymes in the local database
        duplicate_cazyme_entries = []  # Cazyme objects with the same primary GenBank accession

        for genbank_object in primary_genbank_query:
            # retrieve CAZymes with the same primary GenBank accession
            cazyme_query = session.query(Cazyme).\
                filter(Genbank.genbank_accession == primary_genbank).\
                filter(Cazymes_Genbanks.primary == True).all()

            if len(cazyme_query) == 0:
                logger.warning(
                    f"The GenBank accession {primary_genbank} id={genbank_object.genbank_id},\n"
                    "was previously add to the local database as a primary GenBank accession\n"
                    "but NOT associated with a CAZyme."
                )

            elif len(cazyme_query) == 1:
                duplicate_cazyme_entries.append(cazyme_query[0])

            else:
                # multiple CAZymes have the same primary GenBank accession
                logger.warning(
                    "The following CAZymes have the GenBank accession "
                    f"{genbank_object.genbank_accession}, "
                    f"genbank_id={genbank_object.genbank_id}\n, as their primary GenBank accession,"
                    "infering they are duplicates of each other"
                )
                for cazyme in cazyme_query:
                    logger.warning(
                        f"Duplicate CAZymes: {cazyme.cazyme_name}, cazyme_id={cazyme.cazyme_id}"
                    )
                    duplicate_cazyme_entries.append(cazyme)

        # check how many CAZyme duplicates there are with the same primary GenBank accession
        if len(duplicate_cazyme_entries) == 0:
            logger.warning(
                f"The priary GenBank accession {primary_genbank} has been entered multiple times "
                "into the local CAZy database with different 'genbank_id's,\n"
                "but NOT associated with any CAZymes.\n Adding as new protein to the database."
            )
            add_new_protein_to_db(
                cazyme_name,
                family,
                source_organism,
                primary_genbank,
                session,
                ec_numbers,
                genbank_accessions,
                uniprot_accessions,
                pdb_accessions,
            )

        elif len(duplicate_cazyme_entries) == 1:
            # only one CAZyme found associated with the primary GenBank accession
            # add the protein data to this CAZyme
            logger.warning(
                f"The priary GenBank accession {primary_genbank} has been entered multiple times "
                "into the local CAZy database with different 'genbank_id's,\n"
                f"BUT associated with only ONE CAZyme: {duplicate_cazyme_entries[0].cazyme_name} "
                f"id={duplicate_cazyme_entries[0].cazyme_id}"
            )
            add_data_to_protein_record(
                duplicate_cazyme_entries[0],
                family,
                source_organism,
                session,
                ec_numbers,
                genbank_accessions,
                uniprot_accessions,
                pdb_accessions,
            )

        else:
            logger.warning(
                "The following CAZymes appear to be duplicates in the local database, identified "
                f"by same primary GenBank accession number\n, {primary_genbank}.\n"
            )
            for cazyme in duplicate_cazyme_entries:
                logger.warning(
                    f"Duplicate CAZyme: {cazyme.cazyme_name}, id={cazyme.cazyme_id}"
                )
            logger.warning(
                f"Protein data added to {duplicate_cazyme_entries[0].cazyme_name} "
                f"id={duplicate_cazyme_entries[0].cazyme_id}"
            )
            # add the protein data scraped from CAZy to the first duplicate CAZyme
            add_data_to_protein_record(
                duplicate_cazyme_entries[0],
                family,
                source_organism,
                session,
                ec_numbers,
                genbank_accessions,
                uniprot_accessions,
                pdb_accessions,
            )

    return


def add_new_protein_to_db(
    cazyme_name,
    family,
    source_organism,
    primary_genbank,
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
    :param session: open sqlalchemy session to database

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param genbank_accessions: list of non-primary GenBank accessions
    :param uniprot_accessions: list, accessions of associated records in UniProtKB
    :param pdb_accessions: list, accessions of associated records in PDB

    Return nothing.
    """
    logger = logging.getLogger(__name__)

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
            f"{len(tax_query)} duplicate records of species {genus} {species} found in the local "
            f"database.\n Adding CAZyme to the record tax_id={tax_query[0].taxonomy_id}"
        )
        new_cazyme.taxonomy = tax_query[0]

    # add the new cazyme to the database
    session.add(new_cazyme)
    session.commit()

    # add primary GenBank accession (new CAZyme, therefore new primary GenBank accession)
    new_primary_genbank = Genbank(genbank_accession=primary_genbank)
    session.add(new_primary_genbank)
    session.commit()

    # establish relationship between the CAZyme and its primary GenBank accession
    relationship = Cazymes_Genbanks(cazymes=new_cazyme, genbanks=new_primary_genbank, primary=True)
    session.add(relationship)
    session.commit()

    # add Family/Subfamily classifications
    if family.find("_") != -1:
        add_cazy_subfamily(family, new_cazyme, session)
    else:
        add_cazy_family(family, new_cazyme, session)

    # define ec_number
    if len(ec_numbers) != 0:
        add_ec_numbers(ec_numbers, new_cazyme, session)

    # add non-primary GenBank accessions
    if len(genbank_accessions) != 0:
        add_genbank_accessions(genbank_accessions, new_cazyme, session)

    # add UniProt accessions
    if len(uniprot_accessions) != 0:
        add_uniprot_accessions(uniprot_accessions, new_cazyme, session)

    # add PDB/3D accessions
    if len(pdb_accessions) != 0:
        add_pdb_accessions(pdb_accessions, new_cazyme, session)

    # final commit to ensure all changes are commited
    session.commit()

    return


def add_data_to_protein_record(
    cazyme,
    family,
    source_organism,
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

    :param cazyme_name: Cazyme class instance - class imported from scraper.sql.sql_orm
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
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
        add_cazy_subfamily(family, cazyme, session)
    else:
        add_cazy_family(family, cazyme, session)

    # define ec_number
    if len(ec_numbers) != 0:
        add_ec_numbers(ec_numbers, cazyme, session)

    # add non-primary GenBank accessions
    if len(genbank_accessions) != 0:
        add_genbank_accessions(genbank_accessions, cazyme, session)

    # add UniProt accessions
    if len(uniprot_accessions) != 0:
        add_uniprot_accessions(uniprot_accessions, cazyme, session)

    # add PDB/3D accessions
    if len(pdb_accessions) != 0:
        add_pdb_accessions(pdb_accessions, cazyme, session)

    # final commit to ensure all changes are commited
    session.commit()

    return


def add_cazy_family(family, cazyme, session):
    """Add CAZy family to CAZyme record in the local CAZy database.

    :param family: str, name of a CAZy family/subfamily
    :param cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    logger = logging.getLogger(__name__)

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
            # The current family record is associated with a subfamily, the current working family
            # is not therefore, add new CAZy family without a subfamily association to the database
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
                f"Adding CAZyme {cazyme.cazyme_name} id={cazyme.cazyme_id} to family id="
                f"{nonsubfamily_asociated_family_entries[0].family_id}"
            )
            if not(nonsubfamily_asociated_family_entries[0] in cazyme.families):
                cazyme.families.append(nonsubfamily_asociated_family_entries[0])
            session.commit()

    # final commit to ensure all changes are added to the db
    session.commit()

    return


def add_cazy_subfamily(family, cazyme, session):
    """Add CAZy family to CAZyme record in the local CAZy database.

    :param family: str, name of a CAZy family/subfamily
    :param cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    logger = logging.getLogger(__name__)

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
            f"Adding CAZyme {cazyme.cazyme_name} id={cazyme.cazyme_id} to family "
            f"id={subfam_query[0].family_id}"
        )
        if not(subfam_query[0] in cazyme.families):
            cazyme.families.append(subfam_query[0])
        session.commit()

    session.commit()
    return


# Note: Not all CAZymes have non-primary GenBank accessions, EC numbers, UniProt accessions or
# PDB accessions


def add_genbank_accessions(genbank_accessions, cazyme, session):
    """Add non-primary GenBank protein accessions to the local database.

    :param genbank_accessions: list of non-primary GenBank protein accession numbers (str)
    :param cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    cazyme_id = cazyme.cazyme_id

    for accession in genbank_accessions:
        # check if accession is in the database already
        genbank_query = session.query(Genbank).filter_by(genbank_accession=accession).all()

        if len(genbank_query) == 0:
            # add new genbank accession
            new_accession = Genbank(genbank_accession=accession)
            session.add(new_accession)
            session.commit()

            # establish relationship between the CAZyme and its primary GenBank accession
            genbank_id = new_accession.genbank_id
            caz_gen_query = session.(Cazymes_Genbanks).\
                filter(Cazymes_Genbanks.cazyme_id == cazyme_id).\
                filter(Cazymes_Genbanks.genbank_id == genbank_id).\
                all()

            if len(caz_gen_query) == 0:
                relationship = Cazymes_Genbanks(
                    cazymes=cazyme,
                    genbanks=new_accession,
                    primary=False,
                )
                session.add(relationship)
                session.commit()

        elif len(genbank_query) == 1:
            # add GenBank record to current working CAZyme
            caz_gen_query = session.(Cazymes_Genbanks).\
                filter(Cazymes_Genbanks.cazyme_id == cazyme_id).\
                filter(Cazymes_Genbanks.genbank_id == genbank_query[0].genbank_id).\
                all()

            if len(caz_gen_query) == 0:
                relationship = Cazymes_Genbanks(
                    cazymes=cazyme,
                    genbanks=genbank_query[0],
                    primary=False,
                )
                session.add(relationship)
                session.commit()

        else:
            logger.warning(
                f"Duplicate entries for GenBank accession {accession}"
                f"Adding the accession entry with the ID {genbank_query[0].genbank_id} to the\n"
                f"cazyme {cazyme.cazyme_name} id={cazyme.cazyme_id}"
            )
            caz_gen_query = session.(Cazymes_Genbanks).\
                filter(Cazymes_Genbanks.cazyme_id == cazyme_id).\
                filter(Cazymes_Genbanks.genbank_id == genbank_query[0].genbank_id).\
                all()

            if len(caz_gen_query) == 0:
                relationship = Cazymes_Genbanks(
                    cazymes=cazyme,
                    genbanks=genbank_query[0],
                    primary=False,
                )
                session.add(relationship)
                session.commit()

    session.commit()
    return


def add_ec_numbers(ec_numbers, cazyme, session):
    """Add EC numbers to CAZyme record in the local CAZy database.

    :param ec_numbers: list of EC numbers (str)
    :param new_cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    logger = logging.getLogger(__name__)

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
                f"Duplicate entries found for EC{ec} in the local CAZy database."
                f"Adding EC{ec} id={ec_query[0].ec_id} to the cazyme {cazyme.cazyme_name} "
                f"id={cazyme.cazyme_id}"
            )
            if not(ec_query[0] in cazyme.ecs):
                cazyme.ecs.append(ec_query[0])
            session.commit()

    session.commit()
    return


def add_uniprot_accessions(uniprot_accessions, cazyme, session):
    """Add UniProt protein accessions to CAZyme record in the local CAZy database.

    :param uniprot_accessions: list of UniProt protein accession numbers (str)
    :param cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    if len(uniprot_accessions) == 1:
        add_primary_uniprot(uniprot_accessions[0], cazyme, session)
        # check if accession is already in the database

    else:
        add_primary_uniprot(uniprot_accessions[0], cazyme, session)

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
                    f"Adding the UniProt accession id={uniprot_query[0].uniprot_id} to the\n"
                    f"cazyme {cazyme.cazyme_name} id={cazyme.cazyme_id}"
                )
                if not(uniprot_query[0] in cazyme.uniprots):
                    cazyme.uniprots.append(uniprot_query[0])
                session.commit()

    session.commit()
    return


def add_primary_uniprot(accession, cazyme, session):
    """Add primary UniProt accession for CAZyme to the local CAZy database.

    :param accession: str, primary UniProt accession of a CAZyme
    :param cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    logger = logging.getLogger(__name__)

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
            f"Adding the UniProt accession id={uniprot_query[0].uniprot_id} to the\n"
            f"cazyme {cazyme.cazyme_name} id={cazyme.cazyme_id}"
        )
        if not(uniprot_query[0] in cazyme.uniprots):
            cazyme.uniprots.append(uniprot_query[0])
        session.commit()

    session.commit()
    return


def add_pdb_accessions(pdb_accessions, cazyme, session):
    """Add PDB/3D protein accessions to CAZyme record in the local CAZy database.

    :param pdb_accessions: list of UniProt protein accession numbers (str)
    :param cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    if len(pdb_accessions) == 1:
        add_primary_pdb(pdb_accessions[0], cazyme, session)
        # check if accession is already in the database

    else:
        add_primary_pdb(pdb_accessions[0], cazyme, session)

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
                    f"Adding PDB id={pdb_query[0].pdb_id} to the\n"
                    f"cazyme {cazyme.cazyme_name} id={cazyme.cazyme_id}"
                )
                if not(pdb_query[0] in cazyme.pdbs):
                    cazyme.pdbs.append(pdb_query[0])
                session.commit()

    session.commit()
    return


def add_primary_pdb(accession, cazyme, session):
    """Add primary PDB/3D accession for CAZyme to the local CAZy database.

    :param accession: str, primary UniProt accession of a CAZyme
    :param cazyme: Cazymes class object
    :param session: open local database session connector

    Return nothing.
    """
    logger = logging.getLogger(__name__)

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
            f"Duplicate entries for PDB/3D accession {accession}.\n"
            f"Adding PDB id={pdb_query[0].pdb_id} to cazyme {cazyme.cazyme_name} "
            f"id={cazyme.cazyme_id}"
        )
        if not(pdb_query[0] in cazyme.pdbs):
            cazyme.pdbs.append(pdb_query[0])
        session.commit()

    session.commit()
    return
