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

from sqlalchemy import create_engine, Boolean, Column, ForeignKey, Integer, String, Table
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.ext.declarative import declarative_base


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
)

# linker table between cazymes and ec numbers
cazymes_ecs = Table(
    "cazymes_ecs",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("ec_id", Integer, ForeignKey("ecs.ec_id")),
)

# linker table between cazymes and GenBank accession of source protein sequence
cazymes_genbanks = Table(
    "cazymes_genbanks",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("genbank_id", Integer, ForeignKey("genbanks.genbank_id")),
)

# linker table between cazymes and UniProt accessions of CAZymes
cazymes_uniprots = Table(
    "cazymes_uniprots",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("uniprot_id", Integer, ForeignKey("uniprots.uniprot_id")),
)

# linker table between CAZymes and PDB structures
cazymes_pdbs = Table(
    "cazymes_pdbs",
    Base.metadata,
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("pdb_id", Integer, ForeignKey("pdbs.pdb_id")),
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
    # query database to see if CAZyme is already present
    query_result = session.query(Cazyme).filter_by(cazyme_name=cazyme_name).all()

    if len(query_result) == 0:
        add_new_protein_to_db(
            cazyme_name,
            family,
            source_organism,
            ec_numbers,
            external_links,
            logger,
            session,
        )
        return

    for record in query_result:
        if isinstance(record, Cazyme):
            # check if the primary GenBank accession is the same
            duplicate_records = record.genbanks.\
                filter(Cazyme.genbanks.any(genbank_accession=external_links["GenBank"][0])).\
                filter(Genbank.primary==True).\
                all()

            # Different primary GenBank accession identify these CAZymes are
            # derived from different GenBank protein records
            if len(duplicate_records) == 0:
                add_new_protein_to_db(
                    cazyme_name,
                    family,
                    source_organism,
                    ec_numbers,
                    external_links,
                    logger,
                    session,
                )

            elif len(duplicate_records) > 1:
                logger.warning(
                    "Duplicate records found in SQL database,\n"
                    f"under the CAZyme name {cazyme_name} and "
                    f"the primary GenBank accession {external_links['GenBank'][0]}.\n"
                    "The new CAZyme will not be added to the SQL database but\n"
                    "written out to the sql errors log file for manual inspection."
                )

                return(
                    "Duplicate records found in SQL database,\n"
                    f"under the CAZyme name {cazyme_name} and "
                    f"the primary GenBank accession {external_links['GenBank'][0]}.\n"
                    "The new CAZyme was not be added to the SQL database.\n"
                    f"family={family}, source={source_organism}, ec#s={ec_numbers}\n"
                    f"external links={external_links}"
                )

            else:
                add_data_to_protein_record(
                    query_result[0],
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
        uniprot_accessions = external_links["GenBank"]
        add_uniprot_accessions(uniprot_accessions, new_cazyme, session, logger)
    except KeyError:
        pass

    # add PDB/3D accessions
    try:
        pdb_accessions = external_links["GenBank"]
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
    """Add EC numbers to CAZyme record in the local CAZy database.

    :param family: str, name of a CAZy family/subfamily
    :param cazyme: Cazymes class object
    :param session: open local database session connector
    :param logger: logger object

    Return nothing.
    """
    # Check if a subfamily was given
    if family.find("_") != -1:  # cazyme is classifed under a subfamily
        # check if the subfamily is already in the database
        query = session.query(CazyFamily).filter_by(subfamily=family)

        if len(query) == 0:
            # add new CAZy subfamily to the database
            cazy_family = CazyFamily(family=family[:family.find("_")], subfamily=family)
            session.add(cazy_family)
            session.commit()

            cazyme.families.append(cazy_family)
            session.commit()

        elif len(query) == 1:
            # add existing subfamily to the current working CAZyme
            cazyme.families.append(query[0])
            session.commit()

        else:
            logger.warning(
                f"Duplicate subfamily entries found for the CAZy subfamily {family}."""
            )
            cazyme.families.append(query[0])
            session.commit()

    else:
        # check if the family is already in the database
        query = session.query(CazyFamily).filter_by(family=family)

        if len(query) == 0:
            # add new CAZy to the database
            cazy_family = CazyFamily(family=family)
            session.add(cazy_family)
            session.commit()

            cazyme.families.append(cazy_family)
            session.commit()

        elif len(query) == 1:
            # check if it is associated with a CAZy subfamily
            if query.subfamily is not None:
                # add new CAZy family without a subfamily association to the database
                cazy_family = CazyFamily(family=family)
                session.add(cazy_family)
                session.commit()

                cazyme.families.append(cazy_family)
                session.commit()

        else:
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

                cazyme.families.append(cazy_family)
                session.commit()

            else:
                logger.warning(
                    f"Duplicate family entries found for the CAZy family {family} "
                    "without subfamily association."""
                )
                cazyme.families.append(entry[0])
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

            cazyme.ecs.append(new_ec)
            session.commit()

        elif len(query) == 1:
            # add existing EC number record to the CAZyme record
            cazyme.ecs.append(query[0])
            session.commit()

        else:
            # duplicate records found for the current EC number
            logger.warning(
                f"Duplicate entries found for the EC# {ec} in the local CAZy database."
            )
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
        add_primary_genbank(genbank_accessions, cazyme, session, logger)
        # check if accession is already in the database

    else:
        add_primary_genbank(genbank_accessions[0], cazyme, session, logger)

        for accession in genbank_accessions[1:]:
            # check if accession is in the database already
            query = session.query(Genbank).filter_by(genbank_accession=accession).all()

            if len(query) == 0:
                # add new genbank accession
                new_accession = Genbank(genbank_accession=accession, primary=True)
                session.add(new_accession)
                session.commit()

                cazyme.genbanks.append(new_accession)
                session.commit()

            elif len(query) == 1:
                # add GenBank record to current working CAZyme
                cazyme.genbanks.append(query[0])
                session.commit()

            else:
                logger.warning(
                    f"Duplicate entries for GenBank accession {accession}"
                )
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

        cazyme.genbanks.append(new_accession)
        session.commit()

    elif len(query) == 1:
        # add GenBank record to current working CAZyme
        cazyme.genbanks.append(query[0])
        session.commit()

    else:
        logger.warning(
            f"Duplicate entries for GenBank accession {accession}"
        )
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
        add_primary_uniprot(uniprot_accessions, cazyme, session, logger)
        # check if accession is already in the database

    else:
        add_primary_uniprot(uniprot_accessions[0], cazyme, session, logger)

        for accession in uniprot_accessions[1:]:
            # check if accession is in the database already
            query = session.query(Uniprot).filter_by(uniprot_accession=accession).all()

            if len(query) == 0:
                # add new uniprot accession
                new_accession = Uniprot(uniprot_accession=accession, primary=True)
                session.add(new_accession)
                session.commit()

                cazyme.uniprots.append(new_accession)
                session.commit()

            elif len(query) == 1:
                # add UniProt record to current working CAZyme
                cazyme.uniprots.append(query[0])
                session.commit()

            else:
                logger.warning(
                    f"Duplicate entries for UniProt accession {accession}"
                )
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

        cazyme.uniprots.append(new_accession)
        session.commit()

    elif len(query) == 1:
        # add UniProt record to current working CAZyme
        cazyme.uniprots.append(query[0])
        session.commit()

    else:
        logger.warning(
            f"Duplicate entries for UniProt accession {accession}"
        )
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
        add_primary_pdb(pdb_accessions, cazyme, session, logger)
        # check if accession is already in the database

    else:
        add_primary_pdb(pdb_accessions[0], cazyme, session, logger)

        for accession in pdb_accessions[1:]:
            # check if accession is in the database already
            query = session.query(Pdb).filter_by(pdb_accession=accession).all()

            if len(query) == 0:
                # add new pdb accession
                new_accession = Pdb(pdb_accession=accession, primary=True)
                session.add(new_accession)
                session.commit()

                cazyme.pdbs.append(new_accession)
                session.commit()

            elif len(query) == 1:
                # add PDB/3D record to current working CAZyme
                cazyme.pdbs.append(query[0])
                session.commit()

            else:
                logger.warning(
                    f"Duplicate entries for PDB/3D accession {accession}"
                )
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

        cazyme.pdbs.append(new_accession)
        session.commit()

    elif len(query) == 1:
        # add PDB/3D record to current working CAZyme
        cazyme.pdbs.append(query[0])
        session.commit()

    else:
        logger.warning(
            f"Duplicate entries for PDB/3D accession {accession}"
        )
        cazyme.pdbs.append(query[0])
        session.commit()

    session.commit()
    return
