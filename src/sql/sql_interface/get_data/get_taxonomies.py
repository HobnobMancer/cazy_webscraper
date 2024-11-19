#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
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
#
# Bio.PDB reference:
# Hamelryck, T., Manderick, B. (2003) PDB parser and structure class
# implemented in Python. Bioinformatics 19: 2308–2310
"""Retrieve taxonomy data for proteins matching user criteria"""


import logging
import sys

from tqdm import tqdm

from cazy_webscraper.sql.sql_orm import (
    Ec,
    Genbank,
    Session,
    Taxonomy,
    Uniprot,
)
from cazy_webscraper.sql.sql_interface.get_data.get_selected_gbks import (
    get_class_fam_genbank_accessions
)


def get_taxonomies(
    class_filters,
    family_filters,
    taxonomy_filter_dict,
    kingdom_filters,
    ec_filters,
    connection,
    args,
):
    """Retrieve the taxonomy data for a set of GenBank protein accessions]

    :param genbank_accessions: list of GenBank protein accessions
    :param connection: connection to local sql db
    :param args: cmd-line args parser

    Return dict {genus: {species: {strain}}}
    """
    logger = logging.getLogger(__name__)

    organisms = []

    if args.genbank_accessions is not None or args.uniprot_accessions is not None:
        gbk_tax_dict, uniprt_tax_dict = get_uni_gbk_tax_dict(connection)

        if args.genbank_accessions is not None:
            logger.info(
                "Retrieving PDB structures for GenBank accessions "
                f"listed in\n{args.genbank_accessions}"
            )
            organisms += get_taxs_for_user_gbks(gbk_tax_dict, args)

        if args.uniprot_accessions is not None:
            logger.info(
                "Extracting protein sequences for UniProt accessions "
                f"listed in\n{args.uniprot_accessions}"
            )
            organisms += get_taxs_for_uniprots(uniprt_tax_dict, args)

    else:
        organisms = get_filtered_taxs(
            class_filters,
            family_filters,
            taxonomy_filter_dict,
            kingdom_filters,
            ec_filters,
            connection,
        )

    return organisms


def get_uni_gbk_tax_dict(connection):
    """Compile a dict of the data in the Genbanks and Taxs tables

    :param connection: open connection to an SQLite3 database

    Return dict {genbank_accession: organism}
    """
    with Session(bind=connection) as session:
        uniprt_gbk_tax_records = session.query(Uniprot, Genbank, Taxonomy).\
            join(Genbank, (Taxonomy.taxonomy_id == Genbank.taxonomy_id)).\
            join(Uniprot, (Genbank.uniprot_id == Uniprot.uniprot_id)).\
            all()

    gbk_tax_dict = {}  # {genbank_accession: organism}

    uniprt_tax_dict = {}  # {uniprot acc: organism}

    for record in uniprt_gbk_tax_records:
        uniprt_tax_dict[f"{record[0].genbank_accession}"] = f"{record[-1].genus} {record[-1].species}"
        gbk_tax_dict[f"{record[1].genbank_accession}"] = f"{record[-1].genus} {record[-1].species}"

    return gbk_tax_dict, uniprt_tax_dict


def get_taxs_for_user_gbks(gbk_tax_dict, args):
    """Extract taxonomy accs for GenBank accessions listed in a file

    :param gbk_tax_dict: dict {genbank_accession: organism}
    :param args: cmd-line args parser

    Return list of organisms
    """
    logger = logging.getLogger(__name__)

    try:
        with open(args.genbank_accessions, "r") as fh:
            lines = fh.read().splitlines()
    except FileNotFoundError:
        logger.warning(
            f"Could not find file of GenBank accessions at {args.genbank_accessions}\n"
            "Check the path is correct\n"
            "Terminating program"
        )
        sys.exit(1)

    gbk_accessions = [line.strip() for line in lines]

    organisms = set()

    for gbk_accession in tqdm(gbk_accessions, desc="Getting database IDs for provided GenBank IDs"):
        try:
            organisms.add(gbk_tax_dict[gbk_accession])
        except KeyError:
            logging.warning(
                f"Genbank accession {gbk_accession} retrieved from list in file\n"
                "But accession not in the local CAZyme database\n"
                f"Not extracted protein sequences for {gbk_accession}"
            )

    return list(organisms)


def get_taxs_for_uniprots(uni_tax_dict, args):
    """Extract taxonomy accs for GenBank accessions listed in a file

    :param uni_tax_dict: dict {uniprot_acc: organism}
    :param args: cmd-line args parser

    Return list of organisms
    """
    logger = logging.getLogger(__name__)

    try:
        with open(args.uniprot_accessions, "r") as fh:
            lines = fh.read().splitlines()
    except FileNotFoundError:
        logger.warning(
            f"Could not find file of UniProt accessions at {args.uniprot_accessions}\n"
            "Check the path is correct\n"
            "Terminating program"
        )
        sys.exit(1)

    uniprot_accessions = [line.strip() for line in lines]

    organisms = set()

    for uniprot_accession in tqdm(
        uniprot_accessions,
        desc="Getting organisms for provided UniProt IDs",
    ):
        try:
            organisms.add(uni_tax_dict[uniprot_accession])
        except KeyError:
            logging.warning(
                f"UniProt accession {uniprot_accession} retrieved from list in file\n"
                "But accession not in the local CAZyme database\n"
                f"Not extracted protein sequences for {uniprot_accession}"
            )
            continue

    return list(organisms)


def get_filtered_taxs(
    class_filters,
    family_filters,
    taxonomy_filter_dict,
    kingdom_filters,
    ec_filters,
    connection,
):
    """Retrieve the taxs for proteins matching the user criteria.

    :param class_filters: set of CAZy classes to retrieve data for
    :param family_filters: set of CAZy families to retrieve data for
    :param taxonomy_filters: dict of taxonom filters to limit the retrieval of data to
    :param kingdom_filters: set of tax kingdoms to limit the retrieval of data to
    :param ec_filters: set of EC numbers to limit the retrieval of data to
    :param connection: open sqlaclchemy connection for an SQLite db

    Return dict of organisms
    """
    logger = logging.getLogger(__name__)

    # retrieve GenBank accessions of proteins in user selected CAZy classes and (sub)families
    initially_selected_gbk = get_class_fam_genbank_accessions(
        class_filters,
        family_filters,
        connection,
    )

    if len(initially_selected_gbk) == 0:
        logger.error(
            "Retrieved NO proteins for the user selected CAZy classes and (sub)families\n"
            "Ensure proteins belonging to these classes and (sub)families are catalouged "
            "into the local CAZyme db\n"
            "Terminating program"
        )
        sys.exit(1)

    logger.info(
        f"Retrieved {len(initially_selected_gbk)} from user selected CAZy class and (sub)families"
    )

    # Retrieve the db ID numbers of taxonomy entries matching the users taxonomy/kingdom filters
    filtered_records = apply_tax_filters(
        initially_selected_gbk,
        taxonomy_filter_dict,
        kingdom_filters,
    )

    if len(ec_filters) != 0:
        filtered_records = apply_ec_filters(
            filtered_records,
            ec_filters,
            connection
        )

    logger.info(f"Identified {len(filtered_records)} organisms of interest")

    organisms = {}  # {genus: {species: {strain}}} or {genus: None if virus}

    for record in tqdm(filtered_records, desc="Compiling dict of db taxonomies"):
        genus = record[1].genus
        species = record[1].species

        species_strain = record[1].species
        if species_strain.split(" ")[0] == 'sp.':
            species = species_strain
            strain = ""
        else:
            species = species_strain.split(" ")[0]
            strain = species_strain.replace(species, "").strip()

        # add organism to dict of organisms {genus: {species: {strain}}}
        if record[-1].kingdom == 'Viruses':
            organism = f"{record[1].genus} {record[1].species}"
            organism = organism.split("(")[0].strip()
            organism = " ".join([_.strip() for _ in organism.split(" ") if _.find("/") == -1 and _.find("]") == -1 and _.find("[") == -1]).strip()

            if organism.startswith("Influenza A virus"):
                organism = "Influenza A virus"

            elif organism.startswith("Influenza B virus"):
                organism = "Influenza B virus"

            organisms[organism] = None

        else:
            try:
                organisms[genus]

                try:
                    organisms[genus][species].add(strain)

                except KeyError:
                    organisms[genus][species] = {strain}

            except KeyError:
                organisms[genus] = {species: {strain}}

    return organisms


def apply_tax_filters(initially_selected_records, taxonomy_filters, kingdom_filters):
    """Filter retrieved GenBank accessions by taxonomy filters.

    :param initially_selected_records: list of db objs retrieved from the db
        including a Genbank, Taxonomy and Kingdom record
    :param taxonomy_filters: dict of taxonom filters to limit the retrieval of data to
    :param kingdom_filters: set of tax kingdoms to limit the retrieval of data to
    :param connection: open sqlaclchemy connection for an SQLite db

    Return set of db objs, each containing the gbk obj and tax obj
    """
    logger = logging.getLogger(__name__)

    if len(taxonomy_filters['genera']) == 0 and \
        len(taxonomy_filters['species']) == 0 and \
        len(taxonomy_filters['strains']) == 0 and \
            len(kingdom_filters) == 0:
        logger.warning("Applying no taxonomic filters")
        return list(set(initially_selected_records))

    filtered_records = set()

    if len(kingdom_filters) == 0:
        logger.warning("Not applying kingdom filter(s)")
    for kingdom in tqdm(kingdom_filters, desc="Applying kingdom filter(s)"):
        for obj in initially_selected_records:
            if kingdom == obj[2].kingdom:
                filtered_records.add(obj)

    if len(taxonomy_filters['genera']) == 0:
        logger.warning("Npt applying genera filters")
    for genus in taxonomy_filters['genera']:
        for obj in initially_selected_records:
            if genus == obj[1].genus:
                filtered_records.add(obj)

    if len(taxonomy_filters['species']) == 0:
        logger.warning("Not applying species filters")
    for species in taxonomy_filters['species']:
        for obj in initially_selected_records:
            if species == obj[1].species.split(" "):
                filtered_records.add(obj)

    if len(taxonomy_filters['strains']) == 0:
        logger.warning("Not applying strains filters")
    for strain in taxonomy_filters['strains']:
        for obj in initially_selected_records:
            if strain == obj[1].species:
                filtered_records.add(obj)

    return filtered_records


def apply_ec_filters(
    current_objs,
    ec_filters,
    connection,
):
    """Apply EC number filter to the retrieved Genbank records.

    :param current_gbk_objs: list of db Genbank objs retrieved from the db
    :param ec_filters: set of EC numbers to limit the retrieval of data to
    :param connection: open sqlaclchemy connection for an SQLite db

    Return set of db Genbank objects including associated tax record
    """
    logger = logging.getLogger(__name__)

    filtered_records = set()
    ec_gbk_ids = set()

    # Retrieve all Genbank.genbank_ids for each EC number
    for ec in tqdm(ec_filters, desc="Retrieving gbks for EC# filters"):
        with Session(bind=connection) as session:
            gbk_query = session.query(Genbank.genbank_id).\
                join(Ec, Genbank.ecs).\
                filter(Ec.ec_number == ec).\
                all()

        for gbk_id in gbk_query:
            ec_gbk_ids.add(gbk_id)

    if len(ec_gbk_ids) == 0:
        logger.error(
            "Retrieved NO proteins matching the provided EC numbers\n"
            "Check the local CAZyme db contains the EC numbers provided\n"
            "Terminating program"
        )
        sys.exit(1)

    ec_filtered_gbks = set()

    for record in tqdm(current_objs, desc="Checking gbk records against EC filters"):
        if (record[0].genbank_id,) in ec_gbk_ids:
            filtered_records.add(record)

    return ec_filtered_gbks
