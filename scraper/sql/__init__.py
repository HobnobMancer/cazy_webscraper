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
    Column("ID", Integer),
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("family_id", Integer, ForeignKey("families.family_id")),
    Column("subfamily_id", Integer, ForeignKey("subfamilies.subfamily_id")),
)

# linker table between cazymes and source organisms
cazymes_taxs = Table(
    "cazymes_taxs",
    Base.metadata,
    Column("ID", Integer),
    Column("cazyme_id", Integer, ForeignKey("cazymes.cazyme_id")),
    Column("taxonomy_id", Integer, ForeignKey("taxs.taxonomy_id")),
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

    families = relationship(
        "CazyFamily",
        secondary=cazymes_families,
        back_populates="cazymes",
        lazy="dynamic",
    )
    subfamilies = relationship(
        "CazySubFamily",
        secondary=cazymes_families,
        back_populates="cazymes",
        lazy="dynamic",
    )
    taxs = relationship(
        "Taxonomy",
        secondary=cazymes_taxs,
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
        return f"<CAZyme name={self.name}, id={self.cazyme_id}>"


class CazyFamily(Base):
    """Describes the source organism of CAZymes."""
    __tablename__ = "families"
    family_id = Column(Integer, primary_key=True)
    name = Column(String)

    cazymes = relationship(
        "Cazyme",
        secondary=cazymes_families,
        back_populates="families",
        lazy="dynamic",
    )
    subfamilies = relationship("CazySubFamily")

    def __repr__(self):
        """Return string representation of source organism."""
        return f"<CAZy family id={self.family_id} name={self.family_name}>"


class CazySubFamily(Base):
    """Describes the source organism of CAZymes."""
    __tablename__ = "subfamilies"
    subfamily_id = Column(Integer, primary_key=True)
    name = Column(String)
    family_id = Column(Integer, ForeignKey("families.family_id"))

    cazymes = relationship(
        "Cazyme",
        secondary=cazymes_families,
        back_populates="subfamilies",
        lazy="dynamic",,
    )

    def __repr__(self):
        return f"<Subfamily {self.subfamily_name} under Family {self.family_id}>"


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


def protein_to_db(
    cazyme_name,
    family,
    source_organism,
    ec_numbers,
    external_links,
    session,
):
    """Add protein (CAZyme) data to SQL database (db).

    :param cazyme_name: str, name of the protein/CAZyme
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param external_links: dict, links to external databases
    :param session: open sqlalchemy session to database

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
    if ec_numbers is not None:
        for ec in ec_numbers:
            ec_num = EC(ec_number=ec)
            session.add(ec_num)
            session.commit()
            new_cazyme.ecs.append(ec_num)

    # add CAZy family
    if family.find("_") != -1:  # subfamily
        cazy_family = CazyFamily(name=family[:family.find("_")])
        cazy_subfamily = CazySubFamily(name=family)

        session.add(cazy_subfamily)
        cazy_family.subfamilies = [cazy_subfamily]
        session.add(cazy_family)
        session.commit()

        new_cazyme.families.append(cazy_family)
        new_cazyme.subfamilies.append(cazy_subfamily)

    else:
        cazy_family = CazyFamily(name=family)
        session.add(cazy_family)
        session.commit()

        new_cazyme.families.append(cazy_family)

    # add accession numbers of records linked to the protein/CAZyme in external databases

    # define GenBank accessions
    try:
        genbank_accessions = external_links["GenBank"]

        if len(genbank_accessions) != 0:
            if len(genbank_accessions) > 1:
                new_genbank = Genbank(genbank_accession=genbank_accessions[0], primary=True)
                session.add(new_genbank)
                session.commit()
                new_cazyme.genbanks.append(new_genbank)

                for accession in genbank_accessions[1:]:
                    new_genbank = Genbank(genbank_accession=accession, primary=False)
                    session.add(new_genbank)
                    session.commit()
                    new_cazyme.genbanks.append(new_genbank)

            else:
                new_genbank = Genbank(genbank_accession=genbank_accessions[0], primary=True)
                session.add(new_genbank)
                session.commit()
                new_cazyme.genbanks.append(new_genbank)

    except KeyError:
        pass

    # define UniProt accessions
    try:
        uniprot_accessions = external_links["UniProt"]

        if len(uniprot_accessions) != 0:
            if len(uniprot_accessions) > 1:
                new_uniprot = Genbank(genbank_accession=uniprot_accessions[0], primary=True)
                session.add(new_uniprot)
                session.commit()
                new_cazyme.genbanks.append(new_uniprot)

                for accession in uniprot_accessions[1:]:
                    new_uniprot = Genbank(genbank_accession=accession, primary=False)
                    session.add(new_uniprot)
                    session.commit()
                    new_cazyme.genbanks.append(new_uniprot)

            else:
                new_uniprot = Genbank(genbank_accession=uniprot_accessions[0], primary=True)
                session.add(new_uniprot)
                session.commit()
                new_cazyme.genbanks.append(new_uniprot)

    except KeyError:
        pass

    # define PDB accessions
    try:
        pdb_accessions = external_links["PDB/3D"]

        if len(pdb_accessions) != 0:
            if len(pdb_accessions) > 1:
                new_pdb = Genbank(genbank_accession=pdb_accessions[0], primary=True)
                session.add(new_pdb)
                session.commit()
                new_cazyme.genbanks.append(new_pdb)

                for accession in pdb_accessions[1:]:
                    new_pdb = Genbank(genbank_accession=accession, primary=False)
                    session.add(new_pdb)
                    session.commit()
                    new_cazyme.genbanks.append(new_pdb)

            else:
                new_pdb = Genbank(genbank_accession=pdb_accessions[0], primary=True)
                session.add(new_pdb)
                session.commit()
                new_cazyme.genbanks.append(new_pdb)

    except KeyError:
        pass

    # final commit to ensure all changes are commited
    session.commit()

    return
