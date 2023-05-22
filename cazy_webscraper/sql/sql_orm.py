#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Submodule to build a local SQL database"""

from distutils.command.config import LANG_EXT
import logging
import re

import sqlite3

from sqlalchemy import (
    Column,
    ForeignKey,
    Index,
    Integer,
    PrimaryKeyConstraint,
    String,
    Table,
    UniqueConstraint,
    MetaData,
    create_engine,
    event,
    exc,
)
from sqlalchemy.engine import Engine
from sqlalchemy.ext.compiler import compiles
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.sql.expression import BinaryExpression, func, literal
from sqlalchemy.sql.operators import custom_op


# Use the declarative system
# Database structured in NF1
metadata_obj = MetaData()
Base = declarative_base()
Session = sessionmaker()


# Enable regular expression searching of the database
class ReString(String):
    """Enchanced version of standard SQLAlchemy's :class:`String`.

    Supports additional operators that can be used while constructing filter expressions.
    """
    class comparator_factory(String.comparator_factory):
        """Contains implementation of :class:`String` operators related to regular expressions."""
        def regexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('~'))

        def iregexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('~*'))

        def not_regexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('!~'))

        def not_iregexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('!~*'))


class RegexMatchExpression(BinaryExpression):
    """Represents matching of a column againsts a regular expression."""


@compiles(RegexMatchExpression, 'sqlite')
def sqlite_regex_match(element, compiler, **kw):
    """Compile the SQL expression representing a regular expression match for the SQLite engine."""
    # determine the name of a custom SQLite function to use for the operator
    operator = element.operator.opstring
    try:
        func_name, _ = SQLITE_REGEX_FUNCTIONS[operator]
    except (KeyError, ValueError) as e:
        would_be_sql_string = ' '.join((compiler.process(element.left),
                                        operator,
                                        compiler.process(element.right)))
        raise exc.StatementError(
            f"unknown regular expression match operator: {operator} {would_be_sql_string} {e}"
        )

    # compile the expression as an invocation of the custom function
    regex_func = getattr(func, func_name)
    regex_func_call = regex_func(element.left, element.right)
    return compiler.process(regex_func_call)


@event.listens_for(Engine, 'connect')
def sqlite_engine_connect(dbapi_connection, connection_record):
    """Listener for the event of establishing connection to a SQLite database.

    Creates the functions handling regular expression operators
    within SQLite engine, pointing them to their Python implementations above.
    """
    if not isinstance(dbapi_connection, sqlite3.Connection):
        return

    for name, function in SQLITE_REGEX_FUNCTIONS.values():
        dbapi_connection.create_function(name, 2, function)


# Mapping from the regular expression matching operators
# to named Python functions that implement them for SQLite.
SQLITE_REGEX_FUNCTIONS = {
    '~': ('REGEXP',
          lambda value, regex: bool(re.match(regex, value))),
    '~*': ('IREGEXP',
           lambda value, regex: bool(re.match(regex, value, re.IGNORECASE))),
    '!~': ('NOT_REGEXP',
           lambda value, regex: not re.match(regex, value)),
    '!~*': ('NOT_IREGEXP',
            lambda value, regex: not re.match(regex, value, re.IGNORECASE)),
}


# define linker/relationship tables


genbanks_genomes = Table(
    'Genbanks_Genomes',
    Base.metadata,
    Column("genbank_id", Integer, ForeignKey("Genbanks.genbank_id")),
    Column("genome_id", Integer, ForeignKey("Genomes.genome_id")),
    PrimaryKeyConstraint("genbank_id", "genome_id"),
)


genbanks_families = Table(
    'Genbanks_CazyFamilies',
    Base.metadata,
    Column("genbank_id", Integer, ForeignKey("Genbanks.genbank_id")),
    Column("family_id", Integer, ForeignKey("CazyFamilies.family_id")),
    PrimaryKeyConstraint("genbank_id", "family_id"),
)

genbanks_ecs = Table(
    "Genbanks_Ecs",
    Base.metadata,
    Column("genbank_id", Integer, ForeignKey("Genbanks.genbank_id")),
    Column("ec_id", Integer, ForeignKey("Ecs.ec_id")),
    PrimaryKeyConstraint("genbank_id", "ec_id"),
)


genbanks_pdbs = Table(
    "Genbanks_Pdbs",
    Base.metadata,
    Column("genbank_id", Integer, ForeignKey("Genbanks.genbank_id")),
    Column("pdb_id", Integer, ForeignKey("Pdbs.pdb_id")),
    PrimaryKeyConstraint("genbank_id", "pdb_id"),
)


# Define class tables
class Genbank(Base):
    """Represents a protein GenBank accession number and protein seq.

    The GenBank accession is used to identify unique proteins in the database.
    """
    __tablename__ = 'Genbanks'

    __table_args__ = (
        UniqueConstraint("genbank_accession"),
    )

    genbank_id = Column(Integer, primary_key=True)
    genbank_accession = Column(String, index=True)
    sequence = Column(ReString)
    seq_update_date = Column(ReString)
    taxonomy_id = Column(Integer, ForeignKey("Taxs.taxonomy_id"))
    ncbi_tax_id = Column(Integer, ForeignKey("NcbiTaxs.ncbi_tax_id"))
    uniprot_id = Column(Integer, ForeignKey("Uniprots.uniprot_id"))

    uniprot = relationship(
        "Uniprot",
        back_populates="genbank",
    )

    ncbi_taxs = relationship(
        "NcbiTax",
        back_populates="genbanks",
    )

    genomes = relationship(
        "Genome",
        secondary=genbanks_genomes,
        back_populates="genbanks",
        lazy="dynamic",
    )

    organism = relationship(
        "Taxonomy",
        back_populates="genbanks",
    )

    families = relationship(
        "CazyFamily",
        secondary=genbanks_families,
        back_populates="genbanks",
        lazy="dynamic",
    )

    ecs = relationship(
        "Ec",
        secondary=genbanks_ecs,
        back_populates="genbanks",
        lazy="dynamic",
    )

    pdbs = relationship(
        "Pdb",
        secondary=genbanks_pdbs,
        back_populates="genbank",
        lazy="dynamic",
    )

    def __str__(self):
        return f"-Genbank accession={self.genbank_accession}-"

    def __repr__(self):
        return f"<Class GenBank acc={self.genbank_accession}>"


class Taxonomy(Base):
    """Represent the taxonomy of an organism."""
    __tablename__ = "Taxs"

    __table_args__ = (
        UniqueConstraint("genus", "species"),
        Index("organism_option", "taxonomy_id", "genus", "species")
    )

    taxonomy_id = Column(Integer, primary_key=True)
    genus = Column(String)
    species = Column(String)
    kingdom_id = Column(Integer, ForeignKey("Kingdoms.kingdom_id"))

    genbanks = relationship("Genbank", back_populates="organism")
    tax_kingdom = relationship("Kingdom", back_populates="taxonomy")

    def __str__(self):
        return f"-Source organism, Genus={self.genus}, Species={self.species}-"

    def __repr__(self):
        return (
            f"<Class Taxonomy: genus={self.genus}, species={self.species}, "
            f"id={self.taxonomy_id}, kndgm={self.kingdom_id}>"
        )


class Kingdom(Base):
    """Describes a taxonomy Kingdom. Data retrieved from NCBI"""
    __tablename__ = "Kingdoms"

    __table_args__ = (
        UniqueConstraint("kingdom"),
    )

    kingdom_id = Column(Integer, primary_key=True)
    kingdom = Column(String)

    taxonomy = relationship("Taxonomy", back_populates="tax_kingdom")

    def __str__(self):
        return f"-Kingdom, kingdom={self.kingdom}-"

    def __repr__(self):
        return f"<Class Kingdom, kingdom={self.kingdom}, kingdom_id={self.kingdom_id}>"


class Genome(Base):
    """Represent the genomic assembly from which a protein is sourced."""
    __tablename__ = "Genomes"

    __table_args__ = (
        UniqueConstraint("assembly_name", "gbk_version_accession", "refseq_version_accession"),
        Index(
            "genome_options", "assembly_name", "gbk_version_accession", "refseq_version_accession"
        )
    )

    genome_id = Column(Integer, primary_key=True)
    assembly_name = Column(String)
    gbk_version_accession = Column(String)
    gbk_ncbi_id = Column(Integer)
    refseq_version_accession = Column(String)
    refseq_ncbi_id = Column(Integer)
    gtdb_tax_id = Column(Integer, ForeignKey("GtdbTaxs.gtdb_tax_id"))
    
    genbanks = relationship(
        "Genbank",
        secondary=genbanks_genomes,
        back_populates="genomes",
        lazy="dynamic",
    )

    genome_gtdb_tax = relationship("GtdbTax", back_populates="gtdb_genomes")

    def __str__(self):
        return f"-Genome, Gbk={self.gkb_version_accession}, RefSeq={self.refseq_version_accession}-"

    def __repr__(self):
        return (
            f"<Class Genome: Gbk={self.gkb_version_accession}, RefSeq={self.refseq_version_accession}>"
        )


class GtdbTax(Base):
    """Represent the GTDB taxonomic lineage."""
    __tablename__ = "GtdbTaxs"

    __table_args__ = (
        UniqueConstraint(
            "kingdom",
            "phylum",
            "tax_class",
            "tax_order",
            "family",
            "genus",
            "species",
            "release",
        ),
    )

    gtdb_tax_id = Column(Integer, primary_key=True)
    kingdom = Column(ReString)
    phylum = Column(ReString)
    tax_class = Column(ReString)
    tax_order = Column(ReString)
    family = Column(ReString)
    genus = Column(ReString)
    species = Column(ReString)
    release = Column(ReString)  # GTDB release

    gtdb_genomes = relationship(
        "Genome",
        back_populates="genome_gtdb_tax",
    )

    def __str__(self):
        return f"-GtdbTax, id={self.gtdb_tax_id}, kingdom={self.kingdom}-"

    def __repr__(self):
        return (
            f"<Class GtdbTax: id={self.gtdb_tax_id}, kingdom={self.kingdom}, genus={self.genus}>"
        )


class CazyFamily(Base):
    """Describes a CAZy family, and subfamily if applicable.

    Every unique CAZy family-subfamily pair is represented as a unique instance
    in the database.
    """
    __tablename__ = "CazyFamilies"

    # define columns before table_args so subfam column can be called
    family_id = Column(Integer, primary_key=True)
    family = Column(ReString, nullable=False)  # make this an ReString later
    subfamily = Column(String, nullable=True)

    __table_args__ = (
        UniqueConstraint("family", "subfamily"),
        Index("fam_index", "family", "subfamily"),
    )

    genbanks = relationship(
        "Genbank",
        secondary=genbanks_families,
        back_populates="families",
        lazy="dynamic",
    )

    def __str__(self):
        return (
            f"-CAZy Family, Family={self.family}, Subfamily={self.subfamily}, id={self.family_id}-"
        )

    def __repr__(self):
        """Return string representation of source organism."""
        return(
            f"<Class Family, family={self.family}, subfamily={self.subfamily}, id={self.family_id}>"
        )


class NcbiTax(Base):
    """Describes a NCBI Taxonomy lineage."""
    __tablename__ = "NcbiTaxs"

    # define columns before table_args so subfam column can be called
    ncbi_tax_id = Column(Integer, primary_key=True)
    kingdom = Column(ReString)
    phylum = Column(ReString)
    tax_class = Column(ReString)
    tax_order = Column(ReString)
    family = Column(ReString)
    genus = Column(ReString)
    species = Column(ReString)
    strain = Column(ReString)

    __table_args__ = (
        UniqueConstraint("ncbi_tax_id"),
        Index("ncbi_index", "ncbi_tax_id", "genus", "species"),
    )

    genbanks = relationship(
        "Genbank",
        back_populates="ncbi_taxs",
        lazy="dynamic",
    )

    def __str__(self):
        return f"-Ncbi Tax, Kingdom={self.kingdom}, genus={self.genus}, species={self.species}-"

    def __repr__(self):
        """Return string representation of NCBI Tax record."""
        return(
            f"<Class NcbiTax, Kingdom={self.kingdom}, genus={self.genus}, species={self.species}>"
        )


class Uniprot(Base):
    """Table containing UniProt accessions and protein sequences retrieved from UniProtKB"""
    __tablename__ = "Uniprots"

    __table_args__ = (
        UniqueConstraint("uniprot_accession",),
        Index("uniprot_option", "uniprot_id", "uniprot_accession")
    )

    uniprot_id = Column(Integer, primary_key=True)
    uniprot_accession = Column(String)
    uniprot_name = Column(ReString)
    sequence = Column(ReString)
    seq_update_date = Column(ReString)
    uniprot_tax_id = Column(Integer, ForeignKey("UniprotTaxs.uniprot_tax_id"))

    genbank = relationship("Genbank", back_populates="uniprot")
    taxs = relationship("UniprotTax", back_populates="uniprots")

    def __str__(self):
        return (
            f"-Uniprot, accession={self.uniprot_accession}, name={self.uniprot_name}, "
            f"id={self.uniprot_id}-"
        )

    def __repr__(self):
        """Return string representation of source organism."""
        return(
            f"<Uniprot, accession={self.uniprot_accession}, "
            f"name={self.uniprot_name}, id={self.uniprot_id}>"
        )


class UniprotTax(Base):
    """Describes a NCBI Taxonomy lineage."""
    __tablename__ = "UniprotTaxs"

    uniprot_tax_id = Column(Integer, primary_key=True)
    genus = Column(ReString)
    species = Column(ReString)

    __table_args__ = (
        UniqueConstraint("genus", "species"),
        Index("uniprot_tax_index", "uniprot_tax_id", "genus", "species"),
    )

    uniprots = relationship(
        "Uniprot",
        back_populates="taxs",
        lazy="dynamic",
    )

    def __str__(self):
        return f"-UniProt Tax,  genus={self.genus}, species={self.species}-"

    def __repr__(self):
        """Return string representation of UniProt Tax record."""
        return(
            f"<Class UniProtTax, genus={self.genus}, species={self.species}>"
        )


class Ec(Base):
    """Describe EC numbers."""
    __tablename__ = "Ecs"
    __table_args__ = (
        UniqueConstraint("ec_number"),
    )

    ec_id = Column(Integer, primary_key=True)
    ec_number = Column(String, index=True)

    genbanks = relationship(
        "Genbank",
        secondary=genbanks_ecs,
        back_populates="ecs",
        lazy="dynamic",
    )

    def __str__(self):
        return f"-EC{self.ec_number}-ec_id={self.ec_number}-"

    def __repr__(self):
        return f"<Class EC, EC{self.ec_number}, ec_id={self.ec_number}>"


class Pdb(Base):
    """Describe a PDB accession number of protein structure."""
    __tablename__ = "Pdbs"
    __table_args__ = (
        UniqueConstraint("pdb_accession"),
    )

    pdb_id = Column(Integer, primary_key=True)
    pdb_accession = Column(String)

    Index('pdb_idx', pdb_accession)

    genbank = relationship(
        "Genbank",
        secondary=genbanks_pdbs,
        back_populates="pdbs",
        lazy="dynamic",
    )

    def __str__(self):
        return f"-PDB accession={self.pdb_accession}, id={self.pdb_id}-"

    def __repr__(self):
        return f"<Class Pdb accession={self.pdb_accession}, id={self.pdb_id}>"


class Log(Base):
    """Record what data was added to the database and when."""
    __tablename__ = "Logs"

    log_id = Column(Integer, primary_key=True)
    date = Column(String)  # date CAZy scrape was initiated
    time = Column(String)  # time scrape was initated
    database = Column(String)
    retrieved_annotations = Column(String)
    classes = Column(String)  # CAZy classes scraped
    families = Column(String)  # CAZy families scraped
    kingdoms = Column(String)  # Taxonomy Kingdoms to retrieve CAZymes from
    genera_filter = Column(String)
    species_filter = Column(String)
    strains_filter = Column(String)
    ec_filter = Column(String)
    cmd_line = Column(String)  # command line arguments

    def __str__(self):
        return(
            f"Log: date={self.date}, scraped classes={self.classes}, "
            f"scraped families={self.families}, cmd line commands={self.cmd_line}"
        )

    def __repr__(self):
        return(
            f"Class Log: date={self.date}, scraped classes={self.classes}, "
            f"scraped families={self.families}, cmd line commands={self.cmd_line}>"
        )


def get_db_connection(db_path, sql_echo, new=False):
    """Create open connection to local CAZy SQL database.

    :param db_path: cmd args parser
    :param sql_echo: bool, value to set to sql_echo to
    :param new: bool, whether it is a new or an existing database being connected to

    Return an open database session.
    """
    logger = logging.getLogger(__name__)

    if new:
        logger.info("Building a new empty database")
    else:
        logger.info("Opening session to an existing local database")

    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=sql_echo, future=True)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)  # allows for calls to session later on when required
    connection = engine.connect()

    return connection
