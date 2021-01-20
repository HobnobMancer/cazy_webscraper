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
"""Module for building and querying an SQL database."""


import os
import sys

from sqlalchemy import (
    create_engine, Boolean, Column, ForeignKey, Integer, PrimaryKeyConstraint, String, Table
)
from sqlalchemy.exc import InvalidRequestError
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker


# Use the declarative system
# Database structured in NF1
Base = declarative_base()
Session = sessionmaker()


# define association/relationship tables


# linker table between cazymes and CAZy family and subfamilies
cazymes_families = Table(
    "cazymes_families",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("family_id", Integer, ForeignKey("families.family_id")),
    PrimaryKeyConstraint("cazyme_id", "family_id"),
)


# linker table between cazymes and ec numbers
cazymes_ecs = Table(
    "cazymes_ecs",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("ec_id", Integer, ForeignKey("ecs.ec_id")),
    PrimaryKeyConstraint("cazyme_id", "ec_id"),
)


# linker table between cazymes and GenBank accession of source protein sequence
cazymes_genbanks = Table(
    "cazymes_genbanks",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("genbank_id", Integer, ForeignKey("genbanks.genbank_id")),
    PrimaryKeyConstraint("cazyme_id", "genbank_id"),
)


# linker table between cazymes and UniProt accessions of CAZymes
cazymes_uniprots = Table(
    "cazymes_uniprots",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("uniprot_id", Integer, ForeignKey("uniprots.uniprot_id")),
    PrimaryKeyConstraint("cazyme_id", "uniprot_id"),
)


# linker table between CAZymes and PDB structures
cazymes_pdbs = Table(
    "cazymes_pdbs",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("pdb_id", Integer, ForeignKey("pdbs.pdb_id")),
    PrimaryKeyConstraint("cazyme_id", "pdb_id"),
)


# define models


class Cazyme(Base):
    """Describes a CAZyme, a protein single entry in CAZy."""
    __tablename__ = "cazymes"

    cazyme_id = Column(Integer, primary_key=True)
    cazyme_name = Column(String)
    taxonomy_id = Column(Integer, ForeignKey("taxs.taxonomy_id"))

    taxonomy = relationship("Taxonomy", back_populates="cazymes")

    families = relationship(
        "CazyFamily",
        secondary=cazymes_families,
        back_populates="cazymes",
        lazy="dynamic",
    )
    ecs = relationship(
        "EC",
        secondary=cazymes_ecs,
        back_populates="cazymes",
        lazy="dynamic",
    )
    genbanks = relationship(
        "Genbank",
        secondary=cazymes_genbanks,
        back_populates="cazymes",
        lazy="dynamic",
    )
    uniprots = relationship(
        "Uniprot",
        secondary=cazymes_uniprots,
        back_populates="cazymes",
        lazy="dynamic",
    )
    pdbs = relationship(
        "Pdb",
        secondary=cazymes_pdbs,
        back_populates="cazymes",
        lazy="dynamic",
    )

    def __repr__(self):
        """Return string representation of Cazyme table object."""
        return f"<CAZyme name={self.cazyme_name}, id={self.cazyme_id}>"


class Taxonomy(Base):
    """Describes the source organism of CAZymes."""
    __tablename__ = "taxs"

    taxonomy_id = Column(Integer, primary_key=True)
    genus = Column(String)
    species = Column(String)

    cazymes = relationship("Cazyme", back_populates="taxonomy")

    def __repr__(self):
        """Return string representation of source organism."""
        return f"<Source organism, taxonomy, Genus={self.genus}, Species={self.species}>"


class CazyFamily(Base):
    """Describes the source organism of CAZymes."""
    __tablename__ = "families"
    family_id = Column(Integer, primary_key=True)
    family = Column(String)
    subfamily = Column(String)

    cazymes = relationship(
        "Cazyme",
        secondary=cazymes_families,
        back_populates="families",
        lazy="dynamic",
    )

    def __repr__(self):
        """Return string representation of source organism."""
        if self.subfamily is None:
            return f"<CAZy family id={self.family_id} name={self.family}>"
        else:
            return f"<CAZy subfamily id={self.family_id} name={self.subfamily}>"


class EC(Base):
    """Describe EC numbers."""
    __tablename__ = "ecs"

    ec_id = Column(Integer, primary_key=True)
    ec_number = Column(String)

    cazymes = relationship("Cazyme", secondary=cazymes_ecs, back_populates="ecs", lazy="dynamic")

    def __repr__(self):
        """Return string representation of EC number (EC) object."""
        return f"<EC#={self.ec_number}>"


class Genbank(Base):
    """Describe GenBank accession numbers of protein sequences."""
    __tablename__ = "genbanks"

    genbank_id = Column(Integer, primary_key=True)
    genbank_accession = Column(String)
    primary = Column(Boolean)

    cazymes = relationship(
        "Cazyme",
        secondary=cazymes_genbanks,
        back_populates="genbanks",
        lazy="dynamic",
    )

    def __repr__(self):
        """Return representation of GenBank accession number."""
        return f"<GenBank acc={self.genbank_accession}, primary={self.primary}>"


class Uniprot(Base):
    """Describe UniProt accession number."""
    __tablename__ = "uniprots"

    uniprot_id = Column(Integer, primary_key=True)
    uniprot_accession = Column(String)
    primary = Column(Boolean)

    cazymes = relationship(
        "Cazyme",
        secondary=cazymes_uniprots,
        back_populates="uniprots",
        lazy="dynamic",
    )

    def __repr__(self):
        """Return representation of UniProt accession number."""
        return f"<UniProt={self.uniprot_accession}, primary={self.primary}>"


class Pdb(Base):
    """Describe PDB accession number of protein structure."""
    __tablename__ = "pdbs"

    pdb_id = Column(Integer, primary_key=True)
    pdb_accession = Column(String)
    primary = Column(Boolean)

    cazymes = relationship("Cazyme", secondary=cazymes_pdbs, back_populates="pdbs", lazy="dynamic")

    def __repr__(self):
        """Return representation of accesison of PDB structure record."""
        return f"<PDB={self.pdb_accession}, primary={self.primary}>"


def build_db(time_stamp, args, logger):
    """Build an empty SQL database and open a session.

    :param time_stamp: str, date and time stamp of when scrape was initated
    :param args: cmd args parser
    :param logger: logger object

    Return an open database session.
    """
    logger.info("Building empty db to store data")

    # build database engine

    if args.output is sys.stdout:
        # write to cwd, this is deleted in scrape is successful
        cwd = os.getcwd()
        db_path = cwd / f"cazy_scrape_temp_{time_stamp}.db"

    else:
        # write to specified output directory
        db_path = args.output / f"cazy_scrape_{time_stamp}.db"

    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=False)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)

    return Session()


def add_protein_to_db(
    cazyme_name,
    family,
    source_organism,
    ec_numbers,
    external_links,
    logger,
    session,
):
    """Coordinate adding protein (CAZyme) data to the SQL database (db).

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param external_links: dict, links to external databases
    :param session: open sqlalchemy session to database

    Return nothing or duplicate error message.
    """
    # Each unique protein is identified by a unique primary GenBank accession, not its CAZyme name
    # query the local database to see if the current working protein is in the database
    genbank_query = session.query(Genbank).\
        filter_by(genbank_accession=external_links['GenBank"'][0]).\
        filter(Genbank.primary==True).\
        all()

    if len(genbank_query) == 0:
        # GenBank accession was not found as a primary accession in the database
        # Inferring the protein is not in the local database
        add_new_protein_to_db(
            cazyme_name,
            family,
            source_organism,
            ec_numbers,
            external_links,
            logger,
            session,
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
                f"{external_links['GenBank"'][0]}, "
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
            f"{external_links['GenBank"'][0].genbank_accession} found "
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
                    "was previously add to the local database as a primary GenBank accession\n
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
                f"The GenBank accession {external_links['GenBank"'][0]} has been entered multiple\n"
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
                f"The GenBank accession {external_links['GenBank"'][0]} has been entered multiple\n"
                "times as a primary GenBank accession into the local CAZy database with different "
                "'genbank_id's,\n"
                f"but associated with only one CAZyme: {duplicate_cazyme_entries[0].cazyme_name} "
                f"id={duplicate_cazyme_entries[0].cazyme_id}"
            )
            add_data_to_protein_record(
                duplicate_cazyme_entries[0],
                family,
                external_links,
                logger,
                session,
            )

        else:
            logger.warning(
                "The following CAZymes appear to be duplicates in the local database becuase\n"
                "they have the same primary GenBank accession number, "
                f"{external_links['GenBank"'][0]}.\n"
                "Although the GenBank accessions have different genbank_ids, therefore, the\n"
                "accession has been added to the database as a primary accession multiple times.\n"
                f"Protein data being added to {duplicate_cazyme_entries[0].cazyme_name} "
                f"id={duplicate_cazyme_entries[0].cazyme_id}""
            )
            for cazyme in duplicate_cazyme_entries:
                logger.warning(
                    f"Duplicate CAZyme: {cazyme.cazyme_name}, id={cazyme.cazyme_id}"
                )
            # add the protein data scraped from CAZy to the first duplicate CAZyme
            add_data_to_protein_record(
                duplicate_cazyme_entries[0],
                family,
                external_links,
                logger,
                session,
            )

    return


def add_new_protein_to_db(
    cazyme_name,
    family,
    source_organism,
    ec_numbers,
    external_links,
    logger,
    session,
):
    """Add a new protein (a CAZyme) record to the database.

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param external_links: dict, links to external databases
    :param session: open sqlalchemy session to database

    Return nothing.
    """
    # define the CAZyme
    new_cazyme = Cazyme(cazyme_name=cazyme_name)

    # define source organism
    genus_sp_separator = source_organism.find(" ")
    genus = source_organism[:genus_sp_separator]
    species = source_organism[genus_sp_separator:]

    # check if the source  organism is already catalogued
    query = session.query(Taxonomy).filter_by(genus=genus, species=species).all()

    if len(query) == 0:
        # create Taxonomy model object
        new_cazyme.taxonomy = Taxonomy(genus=genus, species=species)

    elif len(query) == 1:
        # add exiting Taxonomy object to the new cazyme
        new_cazyme.taxonomy = query[0]

    else:
        logger.warning(
            f"The species {genus} {species} has been loaded into the database "
            f"{len(query)} times.\n"
        )
        new_cazyme.taxonomy = query[0]

    # add the new cazyme to the database
    session.add(new_cazyme)
    session.commit()

    # add Family/Subfamily classifications
    if family.find("_") != -1:
        add_cazy_subfamily(family, new_cazyme, session, logger)
    else:
        add_cazy_family(family, new_cazyme, session, logger)

    # define ec_number
    if ec_numbers is not None:
        add_ec_numbers(ec_numbers, new_cazyme, session, logger)

    # add GenBank accessions
    try:
        genbank_accessions = external_links["GenBank"]
        add_genbank_accessions(genbank_accessions, new_cazyme, session, logger)
    except KeyError:
        pass

    # add UniProt accessions
    try:
        uniprot_accessions = external_links["UniProt"]
        add_uniprot_accessions(uniprot_accessions, new_cazyme, session, logger)
    except KeyError:
        pass

    # add PDB/3D accessions
    try:
        pdb_accessions = external_links["PDB/3D"]
        add_pdb_accessions(pdb_accessions, new_cazyme, session, logger)
    except KeyError:
        pass

    # final commit to ensure all changes are commited
    session.commit()

    return


def add_data_to_protein_record(
    cazyme,
    family,
    ec_numbers,
    external_links,
    logger,
    session,
):
    """Add data to an existing record in the SQL database.

    :param cazyme: CAZyme class object
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param external_links: dict, links to external databases
    :param session: open sqlalchemy session to database

    Return nothing.
    """
    # add Family/Subfamily classifications
    if family.find("_") != -1:
        add_cazy_subfamily(family, cazyme, session, logger)
    else:
        add_cazy_family(family, cazyme, session, logger)

    # define ec_number
    if ec_numbers is not None:
        add_ec_numbers(ec_numbers, cazyme, session, logger)

    # add GenBank accessions
    try:
        genbank_accessions = external_links["GenBank"]
        add_genbank_accessions(genbank_accessions, cazyme, session, logger)
    except KeyError:
        pass

    # add UniProt accessions
    try:
        uniprot_accessions = external_links["GenBank"]
        add_uniprot_accessions(uniprot_accessions, cazyme, session, logger)
    except KeyError:
        pass

    # add PDB/3D accessions
    try:
        pdb_accessions = external_links["GenBank"]
        add_pdb_accessions(pdb_accessions, cazyme, session, logger)
    except KeyError:
        pass

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
    query = session.query(CazyFamily).filter_by(family=family).all()

    if len(query) == 0:
        # add new CAZy to the database
        cazy_family = CazyFamily(family=family)
        session.add(cazy_family)
        session.commit()

        if not(cazy_family in cazyme.families):
            cazyme.families.append(cazy_family)
        session.commit()

    elif len(query) == 1:
        # check if it is associated with a CAZy subfamily

        if query[0].subfamily is not None:
            # The current family record is associated with a subfamily
            # the current working family is not therefore
            # add new CAZy family without a subfamily association to the database
            cazy_family = CazyFamily(family=family)
            session.add(cazy_family)
            session.commit()

            if not(query[0] in cazyme.families):
                cazyme.families.append(cazy_family)
            session.commit()

        else:
            # add existing Family record to current working CAZyme
            if not(query[0] in cazyme.families):
                cazyme.families.append(query[0])
            session.commit()

    else:
        # check if any of the multiple Family records are associated with subfamilies or not
        nonsubfamily_asociated_family_entries = []

        for entry in query:
            # check if retrieved family record is associated with a subfamily
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
                "without subfamily association."""
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
    query = session.query(CazyFamily).filter_by(subfamily=family).all()

    if len(query) == 0:
        # add new CAZy subfamily to the database
        cazy_family = CazyFamily(family=family[:family.find("_")], subfamily=family)
        session.add(cazy_family)
        session.commit()

        if not(cazy_family in cazyme.families):
            cazyme.families.append(cazy_family)
        session.commit()

    elif len(query) == 1:
        # add existing subfamily to the current working CAZyme
        if not(query[0] in cazyme.families):
            cazyme.families.append(query[0])
        session.commit()

    else:
        logger.warning(
            f"Duplicate subfamily entries found for the CAZy subfamily {family}."""
        )
        if not(query[0] in cazyme.families):
            cazyme.families.append(query[0])
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
        query = session.query(EC).filter_by(ec_number=ec).all()

        if len(query) == 0:
            # add new EC number to the database
            new_ec = EC(ec_number=ec)
            session.add(new_ec)
            session.commit()

            if not(new_ec in cazyme.ecs):
                cazyme.ecs.append(new_ec)
            session.commit()

        elif len(query) == 1:
            # add existing EC number record to the CAZyme record
            if not(query[0] in cazyme.ecs):
                cazyme.ecs.append(query[0])
            session.commit()

        else:
            # duplicate records found for the current EC number
            logger.warning(
                f"Duplicate entries found for the EC# {ec} in the local CAZy database."
            )
            if not(query[0] in cazyme.ecs):
                cazyme.ecs.append(query[0])
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
            query = session.query(Genbank).filter_by(genbank_accession=accession).all()

            if len(query) == 0:
                # add new genbank accession
                new_accession = Genbank(genbank_accession=accession, primary=False)
                session.add(new_accession)
                session.commit()

                if not(new_accession in cazyme.genbanks):
                    cazyme.genbanks.append(new_accession)
                session.commit()

            elif len(query) == 1:
                # add GenBank record to current working CAZyme
                if not(query[0] in cazyme.genbanks):
                    cazyme.genbanks.append(query[0])
                session.commit()

            else:
                logger.warning(
                    f"Duplicate entries for GenBank accession {accession}"
                )
                if not(query[0] in cazyme.genbanks):
                    cazyme.genbanks.append(query[0])
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
    query = session.query(Genbank).filter_by(genbank_accession=accession).all()

    if len(query) == 0:
        # add new genbank accession
        new_accession = Genbank(genbank_accession=accession, primary=True)
        session.add(new_accession)
        session.commit()

        if not(new_accession in cazyme.genbanks):
            cazyme.genbanks.append(new_accession)
        session.commit()

    elif len(query) == 1:
        # add GenBank record to current working CAZyme
        if not(query[0] in cazyme.genbanks):
            cazyme.genbanks.append(query[0])
        session.commit()

    else:
        logger.warning(
            f"Duplicate entries for GenBank accession {accession}"
        )
        if not(query[0] in cazyme.genbanks):
            cazyme.genbanks.append(query[0])
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
            query = session.query(Uniprot).filter_by(uniprot_accession=accession).all()

            if len(query) == 0:
                # add new uniprot accession
                new_accession = Uniprot(uniprot_accession=accession, primary=False)
                session.add(new_accession)
                session.commit()

                if not(new_accession in cazyme.uniprots):
                    cazyme.uniprots.append(new_accession)
                session.commit()

            elif len(query) == 1:
                # add UniProt record to current working CAZyme
                if not(query[0] in cazyme.uniprots):
                    cazyme.uniprots.append(query[0])
                session.commit()

            else:
                logger.warning(
                    f"Duplicate entries for UniProt accession {accession}"
                )
                if not(query[0] in cazyme.uniprots):
                    cazyme.uniprots.append(query[0])
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
    query = session.query(Uniprot).filter_by(uniprot_accession=accession).all()

    if len(query) == 0:
        # add new uniprot accession
        new_accession = Uniprot(uniprot_accession=accession, primary=True)
        session.add(new_accession)
        session.commit()

        if not(new_accession in cazyme.uniprots):
            cazyme.uniprots.append(new_accession)
        session.commit()

    elif len(query) == 1:
        # add UniProt record to current working CAZyme
        if not(query[0] in cazyme.uniprots):
            cazyme.uniprots.append(query[0])
        session.commit()

    else:
        logger.warning(
            f"Duplicate entries for UniProt accession {accession}"
        )
        if not(query[0] in cazyme.uniprots):
            cazyme.uniprots.append(query[0])
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
            query = session.query(Pdb).filter_by(pdb_accession=accession).all()

            if len(query) == 0:
                # add new pdb accession
                new_accession = Pdb(pdb_accession=accession, primary=False)
                session.add(new_accession)
                session.commit()

            if not(new_accession in cazyme.pdbs):
                cazyme.pdbs.append(new_accession)
                session.commit()

            elif len(query) == 1:
                # add PDB/3D record to current working CAZyme
                if not(query[0] in cazyme.pdbs):
                    cazyme.pdbs.append(query[0])
                session.commit()

            else:
                logger.warning(
                    f"Duplicate entries for PDB/3D accession {accession}"
                )
                if not(query[0] in cazyme.pdbs):
                    cazyme.pdbs.append(query[0])
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
    query = session.query(Pdb).filter_by(pdb_accession=accession).all()

    if len(query) == 0:
        # add new pdb accession
        new_accession = Pdb(pdb_accession=accession, primary=True)
        session.add(new_accession)
        session.commit()

        if not(new_accession in cazyme.pdbs):
            cazyme.pdbs.append(new_accession)
        session.commit()

    elif len(query) == 1:
        # add PDB/3D record to current working CAZyme
        if not(query[0] in cazyme.pdbs):
            cazyme.pdbs.append(query[0])
        session.commit()

    else:
        logger.warning(
            f"Duplicate entries for PDB/3D accession {accession}"
        )
        if not(query[0] in cazyme.pdbs):
            cazyme.pdbs.append(query[0])
        session.commit()

    session.commit()
    return
