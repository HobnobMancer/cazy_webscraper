#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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
"""Submodule to interact with local SQLite database, and adding data other than CAZyme records."""


import logging

from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import aliased
from scraper.sql import sql_orm
from scraper.sql.sql_orm import (
    Cazyme,
    CazyFamily,
    Cazymes_Genbanks,
    Genbank,
)
from scraper.sql.sql_interface.add_cazyme_data import (
    add_new_protein_to_db,
    add_data_to_protein_record,
    add_nonprimary_gbk_accessions,
    add_cazy_family,
    add_ec_numbers,
    add_uniprot_accessions,
    add_pdb_accessions,
)


class SqlInterfaceException(Exception):
    """General exception for SQL interface"""

    def __init__(self, message):
        self.message = message


def log_scrape_in_db(
    time_stamp,
    config_dict,
    taxonomy_filters,
    kingdoms,
    ec_filters,
    session,
    args,
):
    """Add a log of scraping CAZy to the local database.

    :param time_stamp: str, date and time cazy_webscraper was invoked
    :param config_dict: dict of CAZy classes and families to be scraped
    :param taxonomy_filters: dict of genera, species and strains to restrict the scrape to
    :param kingdoms: list of taxonomy Kingdoms to restrict scrape to
    :param ec_filters: set of EC numbers to limit the scrape to
    :param session: open SQL database session
    :param args: cmd arguments

    Return nothing."""
    date = time_stamp[:time_stamp.find("--")]
    time = time_stamp[((time_stamp.find("--")) + 2):].replace("-", ":")

    new_log = sql_orm.Log(date=date, time=time)

    if config_dict is not None:
        # get classes that user named to be scraped
        try:
            classes = config_dict["classes"]
            if classes is not None:
                classes = str(classes).replace("[", "").replace("]", "").replace("'", "")
                new_log.classes = classes
        except KeyError:
            pass

        # create a list of families instructed to be scraped
        families = []
        for key in config_dict:
            if key == "classes":
                continue
            if config_dict[key] is not None:
                families.append(config_dict[key])

        if len(families) != 0:
            families = str(families).replace("[", "").replace("]", "").replace("'", "")
            new_log.families = families

    # get taxonomy filters defined by user, and separate into genera, species and strains
    try:
        if len(taxonomy_filters["genera"]) != 0:
            genera = str(taxonomy_filters["genera"]).replace("[", "").\
                replace("]", "").replace("'", "")
            new_log.genera = genera
    except TypeError:
        pass

    try:
        if len(taxonomy_filters["species"]) != 0:
            species = str(taxonomy_filters["species"])
            species = species.replace("[", "").replace("]", "").replace("'", "")
            new_log.species = species
    except TypeError:
        pass

    try:
        if len(taxonomy_filters["strains"]) != 0:
            strains = str(taxonomy_filters["strains"])
            strains = strains.replace("[", "").replace("]", "").replace("'", "")
            new_log.strains = strains
    except TypeError:
        pass

    # get Taxonomy Kingdoms defined by user to be scraped
    if kingdoms is not None:
        new_log.kingdoms = str(kingdoms).replace("[", "").replace("]", "").replace("'", "")
    else:
        new_log.kingdoms = "ALL (Archaea, Bacteria, Eukaryota, Viruses, Unclassified"

    # get EC numbers defined by user to be scraped
    if len(ec_filters) != 0:
        new_log.ec_numbers = str(ec_filters).replace("[", "").replace("]", "").replace("'", "")

    # retrieve commands from the command line
    cmd_line = ""
    for cmd in [
        [args.classes, " --classes '"],
        [args.families, " --families '"],
        [args.kingdoms, " --kingdoms"],
        [args.genera, " --genera '"],
        [args.species, " --species '"],
        [args.strains, " --strains '"],
        [args.ec, " --ec '"],
        [args.streamline, "--streamline '"],
    ]:
        try:
            cmd_line = cmd_line + cmd[1] + cmd[0] + "'"
        except TypeError:
            pass

    if len(cmd_line) != 0:
        cmd_line = cmd_line.strip()
        new_log.cmd_line = cmd_line

    session.add(new_log)
    session.commit()

    return


def add_protein_to_db(
    cazyme_name,
    family,
    source_organism,
    kingdom,
    primary_genbank,
    session,
    args,
    ec_numbers=[],
    gbk_nonprimary=[],
    uni_primary=[],
    uni_nonprimary=[],
    pdb_accessions=[],
):
    """Coordinate adding protein (CAZyme) data to the SQL database (db).

    Every CAZyme has a cazyme name (called 'protein name' in CAZy), CAZy (sub)family, source
    organism, and primary GenBank accession. The 'protein_name' assigned by CAZy is not unique for
    each protein. The primary GenBank accession is the only hyperlinked GenBank accession for the
    protein, and believed to be used by CAZy to indicate the source GenBank protein record. It may
    be possible that a GenBank accession is the primary accession for one CAZyme and a non-primary
    accession for another, becuase CAZy does not appear to ID unique proteins  by the GenBank
    accession because duplicate CAZyme entries can be found within CAZy.

    EC numbers, UniProt accessions and PDB accessions may not be given. The UniProt and PDB
    accessions are recorded as the primary accessions, and all other listed accessions as
    non-primary accessions.

    There are a minority of CAZymes in CAZy with NO GenBank accessions. These CAZymes are annotated
    with the GenBank accession (str) 'NA' in the local db, and are logged.

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param kingdom: str, taxonomy Kingdom of the source organism
    :param primary_genbank: str, the hyperlinked GenBank accession from CAZy
    :param session: open sqlalchemy session to database
    :param args: cmd-line args parser

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param gbk_nonprimary: list of non-primary GenBank accessions
    :param uni_primary: list, primary accessions of associated records in UniProtKB
    :param uni_nonprimary: list, non-primary accessions of associated records in UniProtKB
    :param pdb_accessions: list, accessions of associated records in PDB

    Return nothing.
    """
    error_message = None
    # Each unique protein is identified by a unique primary GenBank accession
    try:
        primary_genbank_object = Genbank(genbank_accession=primary_genbank)
        session.add(primary_genbank_object)
        session.commit()

    except IntegrityError:  # raised when Genbank accession already in the database
        session.rollback()  # allows continued interation with the database

        # check if assuming for each family a protein appears, it's data is identical
        if args.streamline is not None:
            streamline_addition(
                cazyme_name,
                family,
                source_organism,
                kingdom,
                primary_genbank,
                session,
                args,
                ec_numbers,
                gbk_nonprimary,
                uni_primary,
                uni_nonprimary,
                pdb_accessions,
            )

        error_message = parse_unique_genbank_conflict(
            cazyme_name,
            family,
            source_organism,
            kingdom,
            primary_genbank,
            session,
            ec_numbers,
            gbk_nonprimary,
            uni_primary,
            uni_nonprimary,
            pdb_accessions,
        )

        return

    error_message = add_new_protein_to_db(
        cazyme_name,
        family,
        source_organism,
        kingdom,
        primary_genbank_object,
        session,
        ec_numbers,
        gbk_nonprimary,
        uni_primary,
        uni_nonprimary,
        pdb_accessions,
    )

    if error_message is not None:
        raise SqlInterfaceException(error_message)

    return


def streamline_addition(
    cazyme_name,
    family,
    source_organism,
    kingdom,
    primary_genbank,
    session,
    args,
    ec_numbers=[],
    gbk_nonprimary=[],
    uni_primary=[],
    uni_nonprimary=[],
    pdb_accessions=[],
):
    """Apply assumption that for each family a protein appears in, it's data is identcal.

    Found protein already exists in the local database, check that the new family is associated,
    and ignore all associated data user has specified to ignored in these cases.

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param kingdom: str, taxonomy Kingdom of the source organism
    :param primary_genbank: str, the hyperlinked GenBank accession from CAZy
    :param session: open sqlalchemy session to database
    :param args: cmd-line args parser

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param gbk_nonprimary: list of non-primary GenBank accessions
    :param uni_primary: list, primary accessions of associated records in UniProtKB
    :param uni_nonprimary: list, non-primary accessions of associated records in UniProtKB
    :param pdb_accessions: list, accessions of associated records in PDB

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    error_message = None

    # check all GenBank accessions have been recorded in the db previously
    G = aliased(Genbank)
    CG = aliased(Cazymes_Genbanks)
    C = aliased(Cazyme)
    C1 = aliased(Cazyme)
    CG1 = aliased(Cazymes_Genbanks)
    G1 = aliased(Genbank)
    F = aliased(CazyFamily)

    # retrieve all GenBank accessions for CAZyme identfied by its primary GenBank accession and name
    query = session.query(G.genbank_accession, C.cazyme_id, G1.genbank_accession, F.family).\
        join(CG, (CG.genbank_id == G.genbank_id)).\
        join(C, (C.cazyme_id == CG.cazyme_id)).\
        join(F, C.families).\
        join(C1, (C1.cazyme_id == CG1.cazyme_id)).\
        join(CG1, (CG1.genbank_id == G1.genbank_id)).\
        filter(G.genbank_accession == primary_genbank).\
        filter(C.cazyme_name == cazyme_name).\
        filter(C1.cazyme_id == C.cazyme_id).\
        all()

    if len(query) == 0:  # couldn't find unique CAZyme so parse as per usual
        logger.warning(
            f"Primary GenBank {primary_genbank} added before as primary but for a different CAZyme."
        )
        error_message = parse_unique_genbank_conflict(
            cazyme_name,
            family,
            source_organism,
            kingdom,
            primary_genbank,
            session,
            ec_numbers,
            gbk_nonprimary,
            uni_primary,
            uni_nonprimary,
            pdb_accessions,
        )
        return error_message

    record_accessions = set([_[2] for _ in query])
    record_accessions.remove(primary_genbank)

    for acc in gbk_nonprimary:
        if acc in record_accessions:
            record_accessions.remove(acc)

    cazyme = query[0][1]

    if len(record_accessions) != 0:  # Add new GenBank accession to CAZyme
        add_nonprimary_gbk_accessions(gbk_nonprimary, cazyme, session)

    # check cazyme is associated with the current working family
    current_families = set([_[-1] for _ in query])
    if family not in current_families:
        add_cazy_family(family, cazyme, session)

    # check if remaining associated data was set to be assumed to be identical each time
    # the CAZyme is presented in a HTML table

    if ((args.streamline).find("ec") != -1) and (len(ec_numbers) != 0):
        add_ec_numbers(ec_numbers, cazyme, session)

    if (args.streamline).find("uniprot") != -1:
        if len(uni_primary) != 0:
            add_uniprot_accessions(uni_primary, cazyme, True, session)
        if len(uni_nonprimary) != 0:
            add_uniprot_accessions(uni_nonprimary, cazyme, True, session)

    if ((args.streamline).find("pdb") != -1) and (len(pdb_accessions) != 0):
        add_pdb_accessions(pdb_accessions, cazyme, session)

    return error_message


def parse_unique_genbank_conflict(
    cazyme_name,
    family,
    source_organism,
    kingdom,
    primary_genbank,
    session,
    ec_numbers=[],
    gbk_nonprimary=[],
    uni_primary=[],
    uni_nonprimary=[],
    pdb_accessions=[],
):
    """Called when primary GenBank accession is already in the local database.

    Check if the accession in the database is a primary accession. If it is then add data to the
    existing protein record. If not then add the protein data as a new protein to the database.

    :param cazyme_name: str, name of the protein/CAZyme
    :param family: str, CAZy family or subfamily the protein is catalogued under in CAZy
    :param source_organism: str, the scientific name of the organism from which the CAZy is derived
    :param kingdom: str, taxonomy Kingdom of the source_organism
    :param primary_genbank: str, the hyperlinked GenBank accession from CAZy
    :param session: open sqlalchemy session to database

    ::optional parameters::
    :param ec_numbers: list of EC numbers which the CAZyme is annotated with
    :param gbk_nonprimary: list of non-primary GenBank accessions
    :param uni_primary: list, primary accessions of associated records in UniProtKB
    :param uni_nonprimary: list, non-primary accessions of associated records in UniProtKB
    :param pdb_acccessions: list, accessions of associated records in PDB

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    error_message = None

    # check if the local GenBank object is for a primary GenBank accession
    primary_genbank_query = session.query(Genbank).\
        filter(Genbank.genbank_accession == primary_genbank).\
        filter(Cazymes_Genbanks.primary == True).all()

    if len(primary_genbank_query) == 0:
        # Accession was not found as a primary accession, inferring the protein is not in the db
        error_message = add_new_protein_to_db(
            cazyme_name,
            family,
            source_organism,
            kingdom,
            primary_genbank_query[0],
            session,
            ec_numbers,
            gbk_nonprimary,
            uni_primary,
            uni_nonprimary,
            pdb_accessions,
        )

    else:
        # Primary GenBank accession found, add additional data to the respective CAZyme
        cazyme_query = session.query(Cazyme, Genbank, Cazymes_Genbanks).\
            join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
            join(Cazyme, (Cazyme.cazyme_id == Cazymes_Genbanks.cazyme_id)).\
            filter(Genbank.genbank_accession == primary_genbank).\
            filter(Cazymes_Genbanks.primary == True).all()

        if len(cazyme_query) == 0:
            logger.warning(
                f"GenBank accession {primary_genbank}, id={primary_genbank_query[0].genbank_id}\n"
                "was previously added to the local database, but not associated with a CAZyme\n"
                "Adding protein data as a new CAZyme in the local database."
            )
            error_message = add_new_protein_to_db(
                cazyme_name,
                family,
                source_organism,
                kingdom,
                primary_genbank_query[0],
                session,
                ec_numbers,
                gbk_nonprimary,
                uni_primary,
                uni_nonprimary,
                pdb_accessions,
            )

        elif len(cazyme_query) == 1:
            add_data_to_protein_record(
                cazyme_query[0][0],
                family,
                session,
                ec_numbers,
                gbk_nonprimary,
                uni_primary,
                uni_nonprimary,
                pdb_accessions,
            )

        else:
            # multiple CAZymes have the same primary GenBank accession, inferring duplicates
            logger.warning(
                "Potential duplicate CAZymes found in the local database, because they have the "
                f"same primary GenBank accession {primary_genbank}\nPotential duplicate CAZymes:"
            )
            for cazyme in cazyme_query:
                logger.warning(f"Name={cazyme[0].cazyme_name}, id={cazyme[0].cazyme_id}")
            logger.warning(
                f"Protein data added to cazyme {cazyme_query[0][0].cazyme_name}, "
                f"id={cazyme_query[0][0].cazyme_id}"
            )
            add_data_to_protein_record(
                cazyme_query[0][0],
                family,
                session,
                ec_numbers,
                gbk_nonprimary,
                uni_primary,
                uni_nonprimary,
                pdb_accessions,
            )

    return error_message


def add_deleted_cazy_family(family, session):
    """Add CAZy family to CAZyme record in the local CAZy database.

    :param family: str, name of a CAZy family/subfamily
    :param session: open local database session connector

    Return nothing.
    """
    if family.find("_") != -1:
        subfamily = family[:family.find("_")]
        deleted_fam = CazyFamily(family=family[:family.find("_")], subfamily=subfamily)
    else:
        deleted_fam = CazyFamily(family=family)

    try:
        session.add(deleted_fam)
        session.commit()
        return

    except IntegrityError:
        session.rollback()

    return
