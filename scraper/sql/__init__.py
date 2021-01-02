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
from sqlalchemy.orm import backref, relationship, sessionmaker
from sqlalchemy.ext.declarative import declarative_base


def build_db(time_stamp, args, logger):
    """Build an empty SQL database.

    :param time_stamp: str, date and time stamp of when scrape was initated
    :param args: cmd args parser
    :param logger: logger object

    Return open database session.
    """
    logger.info("Building empty db to store data")

    # build database engine

    if args.output is sys.stdout:
        # write to cwd, this is deleted in scrape is successful
        cwd = os.getcwd()
        db_path = cwd / f"cazy_scrape_temp_{time_stamp}.db"
        engine = create_engine(f"sqlite+pysqlite:///{cwd}", echo=False)

    else:
        # write to specified output directory
        db_path = args.output / f"cazy_scrape_{time_stamp}.db"
        engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=False)

    Base = declarative_base()

    # define all relation tables

    cazyme_tax = Table(
        "cazyme_tax",
        Base.metadata,
        Column("ID", Integer, primary_key=True, autoincrement=True),
        Column("cazymeID", Integer, ForeignKey("cazyme_tb.ID")),
        Column("taxonomyID", Integer, ForeignKey("taxonomy_tb.ID")),
    )

    cazyme_ec = Table(
        "cazyme_ec",
        Base.metadata,
        Column("ID", Integer, primary_key=True, autoincrement=True),
        Column("cazymeID", Integer, ForeignKey("cazyme.ID")),
        Column("ecID", Integer, ForeignKey("ec.ID")),
    )

    cazyme_genbank = Table(
        "cazyme_genbank",
        Base.metadata,
        Column("ID", Integer, primary_key=True, autoincrement=True),
        Column("cazymeID", Integer, ForeignKey("cazyme.ID")),
        Column("genbankID", Integer, ForeignKey("genbank.ID")),
    )

    cazyme_uniprot = Table(
        "cazyme_uniprot",
        Base.metadata,
        Column("ID", Integer, primary_key=True, autoincrement=True),
        Column("cazymeID", Integer, ForeignKey("cazyme.ID")),
        Column("uniprotID", Integer, ForeignKey("uniprot.ID")),
    )

    cazyme_pdb = Table(
        "cazyme_pdb",
        Base.metadata,
        Column("ID", Integer, primary_key=True, autoincrement=True),
        Column("cazymeID", Integer, ForeignKey("cazyme.ID")),
        Column("pdbID", Integer, ForeignKey("pdb.ID")),
    )

    # define all other tables
    class Cazyme(Base):
        __tablename__ = "cazyme"
        ID = Column(Integer, primary_key=True, autoincrement=True)
        name = Column(String)
        ecs = relationship("EC", backref=backref("cazyme"))
        genbanks = relationship("Genbank", secondary=cazyme_genbank, back_populates="cazyme")
        uniprots = relationship("Uniprot", secondary=cazyme_uniprot, back_populates="cazyme")
        pdbs = relationship("Pdb", backref=backref("cazyme"))

        def __repr__(self):
            return f"<CAZyme({self.name})>"

    class Taxonomy(Base):
        __tablename__ = "taxonomy"
        ID = Column(Integer, primary_key=True, autoincrement=True)
        genus = Column(String)
        species = Column(String)
        cazymes = relationship("Cazyme", backref=backref("taxonomy"))

    class EC(Base):
        __tablename__ = "ec"
        ID = Column(Integer, primary_key=True, autoincrement=True)
        ec_number = Column(String)

    class Genbank(Base):
        __tablename__ = "genbank"
        ID = Column(Integer, primary_key=True, autoincrement=True)
        accession = Column(String)
        primary_accession = Column(Boolean)
        cazymes = relationship("Cazyme", secondary=cazyme_genbank, back_populates="genbank")

    class Uniprot(Base):
        __tablename__ = "uniprot"
        ID = Column(Integer, primary_key=True, autoincrement=True)
        accession = Column(String)
        primary_accession = Column(Boolean)
        cazymes = relationship("Cazyme", secondary=cazyme_uniprot, back_populates="genbank")

    class Pdb(Base):
        __tablename__ = "pdb"
        ID = Column(Integer, primary_key=True, autoincrement=True)
        accession = Column(String)

    # create a session
    Session = sessionmaker(bind=engine)
    session = Session()

    return session


def data_to_db(protein, engine):
    """Add protein data to SQL database (db).

    :param protein: Protein class object

    Return nothing.
    """
    # create a connection to the db
    conn = engine.connect()

    # insert data into db
    conn.execute(organism_tb.insert(), "DATA HERE")

    return
