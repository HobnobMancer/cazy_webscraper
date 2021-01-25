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
"""Submodule to build a local SQL database"""


from sqlalchemy import (
    Boolean, Column, ForeignKey, Integer, PrimaryKeyConstraint, String, Table
)
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
