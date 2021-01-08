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
    session,
):
    """Coordinate adding protein (CAZyme) data to the SQL database (db).

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param external_links: dict, links to external databases
    :param session: open sqlalchemy session to database

    Return nothing.
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
            session,
        )

    else:
        for record in query_result:
            if isinstance(record, Cazyme):
                # check if the primary GenBank accession is the same
                duplicate_records = record.genbanks.\
                    filter(Cazyme.genbanks.any(genbank_accession=external_links["GenBank"][])).\
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
                    return "RETURN MESSAGE SO PICKED UP ABOVE AND PROTEIN ADDED TO LOG ERRORS"

                else:
                    add_data_to_protein_record(
                        family,
                        source_organism,
                        ec_numbers,
                        external_links,
                        session,
                    )

    return


def add_data_to_protein_record(
    family,
    source_organism,
    ec_numbers,
    external_links,
    session,
):
    """Add data to an existing record in the SQL database.

    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param external_links: dict, links to external databases
    :param session: open sqlalchemy session to database

    Return nothing.
    """

    # final commit to ensure all changes are commited
    session.commit()

    return


def add_new_protein_to_db(
    cazyme_name,
    family,
    source_organism,
    ec_numbers,
    external_links,
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

    # add taxonomy data to the new cazyme
    new_cazyme.taxonomy = Taxonomy(genus=genus, species=species)

    session.add(new_cazyme)
    session.commit()

    # define ec_number
    if ec_numbers is not None:
        for ec in ec_numbers:
            ec_num = EC(ec_number=ec)
            session.add(ec_num)
            session.commit()
            new_cazyme.ecs.append(ec_num)
        session.commit()

    # add CAZy family
    if family.find("_") != -1:  # cazyme is classifed under a subfamily
        cazy_family = CazyFamily(family=family[:family.find("_")], subfamily=family)
        session.add(cazy_family)
        session.commit()

        new_cazyme.families.append(cazy_family)
        session.commit()

    else:
        cazy_family = CazyFamily(family=family)
        session.add(cazy_family)
        session.commit()

        new_cazyme.families.append(cazy_family)
        session.commit()

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
