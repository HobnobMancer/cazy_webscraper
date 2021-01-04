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

    Base = declarative_base()

    # define association/relationship tables

    # linker table between cazymes and source organisms
    cazymes_taxs = Table(
        "cazymes_taxs",
        Base.metadata,
        Column("ID", Integer),
        Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
        Column("taxonomy_id", Integer, ForeignKey("taxs.taxonomy_id"))
    )

    # linker table between cazymes and ec numbers
    cazymes_ecs = Table(
        "cazymes_ecs",
        Base.metadata,
        Column("ID", Integer),
        Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
        Column("ec_id", Integer, ForeignKey("ecs.ec_id")),
    )

    # linker table between cazymes and GenBank accession of source protein sequence
    cazymes_genbanks = Table(
        "cazymes_genbanks",
        Base.metadata,
        Column("ID", Integer),
        Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
        Column("genbank_id", Integer, ForeignKey("genbanks.genbank_id")),
    )

    # linker table between cazymes and UniProt accessions of CAZymes
    cazymes_uniprots = Table(
        "cazymes_uniprots",
        Base.metadata,
        Column("ID", Integer),
        Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
        Column("uniprot_id", Integer, ForeignKey("uniprots.uniprot_id")),
    )

    # linker table between CAZymes and PDB structures
    cazymes_pdbs = Table(
        "cazymes_pdbs",
        Base.metadata,
        Column("ID", Integer),
        Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
        Column("pdb_id", Integer, ForeignKey("pdbs.pdb_id")),
    )

    # define models

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

    # Add a new book
    add_new_book(
        session,
        author_name="Stephen King",
        book_title="The Stand",
        publisher_name="Random House",
    )

    return
