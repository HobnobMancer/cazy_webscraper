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
        """Describes a CAZyme, a protein single entry in CAZy."""
        __tablename__ = "cazymes"

        cazyme_id = Column(Integer, primary_key=True)
        name = Column(String)

        taxs = relationship("Taxonomy", secondary=cazymes_taxs, back_populates="cazymes", lazy="dynamic")
        ecs = relationship("EC", secondary=cazymes_ecs, back_populates="cazymes", lazy="dynamic")
        genbanks = relationship("Genbank", secondary=cazymes_genbanks, back_populates="cazymes", lazy="dynamic")
        uniprots = relationship("Uniprot", secondary=cazymes_uniprots, back_populates="cazymes", lazy="dynamic")
        pdbs = relationship("Pdb", secondary=cazymes_pdbs, back_populates="cazymes", lazy="dynamic")

        def __repr__(self):
            """Return string representation of Cazyme table object."""
            return f"<CAZyme name={self.name}, id={self.cazyme_id}>"

    class Taxonomy(Base):
        """Describes the source organism of CAZymes."""
        __tablename__ = "taxs"
        taxonomy_id = Column(Integer, primary_key=True)
        genus = Column(String)
        species = Column(String)

        cazymes = relationship("Cazyme", secondary=cazymes_taxs, back_populates="taxs", lazy="dynamic")

        def __repr__(self):
            """Return string representation of source organism."""
            return f"<Source organism, taxonomy, Genus={self.genus}, Species={self.species}>"

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

        cazymes = relationship("Cazyme", secondary=cazymes_genbanks, back_populates="genbanks", lazy="dynamic")

        def __repr__(self):
            """Return representation of GenBank accession number."""
            return f"<GenBank acc={self.genbank_accession}, primary={self.primary}>"

    class Uniprot(Base):
        """Describe UniProt accession number."""
        __tablename__ = "uniprots"
        uniprot_id = Column(Integer, primary_key=True)
        uniprot_accession = Column(String)
        primary = Column(Boolean)

        cazymes = relationship("Cazyme", secondary=cazymes_uniprots, back_populates="uniprots", lazy="dynamic")

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

    # create engine
    Base.metadata.create_all(engine)

    # create database session
    Session = sessionmaker(bind=engine)
    session = Session()

    return session


def protein_to_db(cazyme_name, source_organism, ec_numbers, primary_genbank, genbanks, primary_uniprot, uniprots, primary_pdb, pdbs, session):
    """Add protein (CAZyme) data to SQL database (db).

    :param protein: Protein class object

    Return nothing.
    """
    # define the CAZyme
    new_cazyme = Cazyme(name=cazyme_name)
    session.add(new_cazyme)
    session.commit()

    # define source organism
    genus_sp_separator = source_organism.find(" ")
    genus = source_organism[:genus_sp_separator]
    species = source_organism[genus_sp_separator:]

    organism = Taxonomy(genus=genus, species=species)
    session.add(organism)
    session.commit()
    new_cazyme.taxs.append(organism)

    # define ec_number
    for ec in ec_numbers:
        ec_num = EC(ec_number=ec)
        session.add(ec_num)
        session.commit()
        new_cazyme.ecs.append(ec_num)

    # define genbank accession
    new_genbank = Genbank(genbank_accession=primary_genbank, primary=True)
    session.add(new_genbank)
    session.commit()
    new_cazyme.genbanks.append(new_genbank)

    if genbanks is not None:
        for accession in genbanks:
            new_genbank = Genbank(genbank_accession=accession, primary=False)
            session.add(new_genbank)
            session.commit()
            new_cazyme.genbanks.append(new_genbank)

    # define uniprot accessions
    new_uniprot = Uniprot(uniprot_accession=primary_uniprot, primary=True)
    session.add(new_uniprot)
    session.commit()
    new_cazyme.uniprots.append(new_uniprot)

    if uniprots is not None:
        for accession in uniprots:
            new_uniprot = Uniprot(uniprot_accession=accession, primary=False)
            session.add(new_uniprot)
            session.commit()
            new_cazyme.uniprots.append(new_uniprot)

    # define accessions of pdb structure records
    new_pdb = Pdb(pdb_accession=primary_pdb, primary=True)
    session.add(new_pdb)
    session.commit()
    new_cazyme.pdbs.append(new_pdb)

    if pdbs is not None:
        for accession in pdbs:
            new_pdb = Pdb(pdb_accession=accession, primary=False)
            session.add(new_pdb)
            session.commit()
            new_cazyme.pdbs.append(new_pdb)

    # final commit to ensure all changes are commited
    session.commit()

    return
