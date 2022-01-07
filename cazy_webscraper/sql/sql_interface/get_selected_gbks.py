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
"""Funcs for retrieving GenBank accessions of interest, matching user filter criteria"""


import logging
import sys

from tqdm import tqdm

from cazy_webscraper.sql.sql_orm import (
    CazyFamily,
    Ec,
    Genbank,
    Kingdom,
    Taxonomy,
    Session,
)


def get_genbank_accessions(
    class_filters,
    family_filters,
    taxonomy_filters,
    kingdom_filters,
    ec_filters,
    connection,
):
    """Retrieve the GenBank accessions of proteins of interest (to retrieve data for).
    
    :param class_filters: set of CAZy classes to retrieve data for
    :param family_filters: set of CAZy families to retrieve data for
    :param taxonomy_filters: dict of taxonom filters to limit the retrieval of data to
    :param kingdom_filters: set of tax kingdoms to limit the retrieval of data to
    :param ec_filters: set of EC numbers to limit the retrieval of data to
    :param connection: open sqlaclchemy connection for an SQLite db
    
    Return dict {gbk_acc: gbk_id}
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
            "Ensure proteins belonging to these classes and (sub)families are catalouged into the local CAZyme db\n"
            "Terminating program"
        )
        sys.exit(1)
    
    logger.info(
        f"Retrieved {len(initially_selected_gbk)} from user selected CAZy class and (sub)families"
    )
    
    # Retrieve the db ID numbers of taxonomy entries matching the users taxonomy/kingdom filters
    filtered_gbk_accessions = apply_tax_filters(
        initially_selected_gbk,
        taxonomy_filters,
        kingdom_filters,
        connection
    )
    
    if len(filtered_gbk_accessions) == 0:
        logger.error(
            "Retrieved NO proteins for the user selected taxonomy and kingdom filters\n"
            "Ensure proteins belonging to these taxa are catalouged into the local CAZyme db\n"
            "Terminating program"
        )
        sys.exit(1)
    
    # Apply EC number filter if provided
    if len(ec_filters) != 0:
        filtered_gbk_accessions = apply_ec_filters(
            ec_filters,
            connection,
        )
    
    # extract the accession numbers from the db Genbank objects and their db genbank_id
    gbk_dict = {}
    for obj in filtered_gbk_accessions:
        gbk_dict[obj.genbank_accession] = obj.genbank_id
    
    if len(list(gbk_dict.keys())) == 0:
        logger.error(
            "No proteins in the local CAZyme db matched the provided critiera.\n"
            "Check the critieria matches data in the local CAZyme db.\n"
            "Terminating program"
        )
        sys.exit(1)
    
    return gbk_dict

    
def get_class_fam_genbank_accessions(
    class_filters,
    family_filters,
    connection,
):
    """Retrieve the GenBank accessions of proteins from user selected CAZy classes and (sub)families

    :param class_filters: set of CAZy classes to retrieve data for
    :param family_filters: set of CAZy families to retrieve data for
    :param connection: open sqlaclchemy connection for an SQLite db
    
    Return list of db objects containing a Genbank obj, Taxonomy obj and Kingdom obj.
    """
    initially_selected_gbk = []

    if len(class_filters) == 0 and len(family_filters) == 0:
        # could retrieve all GenBank accessions
        with Session(bind=connection) as session:
            gbk_query = session.query(Genbank, Taxonomy, Kingdom).\
                join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                join(Genbank, (Genbank.taxonomy_id == Taxonomy.taxonomy_id)).\
                join(CazyFamily, Genbank.families).\
                all()

            initially_selected_gbk = gbk_query

    else:
        for cazy_class in tqdm(class_filters, desc="Retrieving GenBank accessions for selected CAZy classes"):
            with Session(bind=connection) as session:
                class_subquery = session.query(Genbank.genbank_id).\
                join(CazyFamily, Genbank.families).\
                filter(CazyFamily.family.like(f'{cazy_class}%')).\
                subquery()

            with Session(bind=connection) as session:
                gbk_query = session.query(Genbank, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Genbank, (Genbank.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(CazyFamily, Genbank.families).\
                    filter(Genbank.genbank_id.in_(class_subquery)).\
                    all()

            initially_selected_gbk += gbk_query
    
        for cazy_fam in tqdm(family_filters, desc="Retrieving GenBank accessions for selected CAZy families"):
            if cazy_fam.find('_') != -1:  # subfamily
                with Session(bind=connection) as session:
                    fam_subquery = session.query(Genbank.genbank_id).\
                    join(CazyFamily, Genbank.families).\
                    filter(CazyFamily.subfamily == cazy_fam).\
                    subquery()

            else:  # family
                with Session(bind=connection) as session:
                    fam_subquery = session.query(Genbank.genbank_id).\
                    join(CazyFamily, Genbank.families).\
                    filter(CazyFamily.family == cazy_fam).\
                    subquery()

            with Session(bind=connection) as session:
                gbk_query = session.query(Genbank, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Genbank, (Genbank.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(CazyFamily, Genbank.families).\
                    filter(Genbank.genbank_id.in_(fam_subquery)).\
                    all()

            initially_selected_gbk += gbk_query
        
    return initially_selected_gbk


def apply_tax_filters(
    initially_selected_gbk,
    taxonomy_filters,
    kingdom_filters,
    connection,
):
    """Filter retrieved GenBank accessions by taxonomy filters.
    
    :param initally_selected_gbk: list of db Genbank objs retrieved from the db
    :param taxonomy_filters: dict of taxonom filters to limit the retrieval of data to
    :param kingdom_filters: set of tax kingdoms to limit the retrieval of data to
    :param connection: open sqlaclchemy connection for an SQLite db
    
    Return set of db Genbank objs
    """
    logger = logging.getLogger(__name__)
    
    if len(taxonomy_filters) == 0 and len(kingdom_filters) == 0:
        gbks = [obj[0] for obj in initially_selected_gbk]
        return set(gbks)
    
    tax_ids = set()
    
    for kingdom in tqdm(kingdom_filters, desc="Retrieving IDs of species from selected kingdoms"):
        with Session(bind=connection) as session:
            kingdom_query = session.query(Taxonomy.taxonomy_id).\
                join(Kingdom, (Kingdom.kingdom_id == Taxonomy.kingdom_id)).\
                filter(Kingdom.kingdom == kingdom).\
                all()
            for taxa in kingdom_query:
                tax_ids.add(taxa[0])

    genera = taxonomy_filters['genus']
    species = taxonomy_filters['species']
    strains = taxonomy_filters['strains']

    for genus in tqdm(genera, desc="Retrieving IDs of species from selected genera"):
        with Session(bind=connection) as session:
            tax_query = session.query(Taxonomy.taxonomy_id).\
                filter(Taxonomy.genus == genus).\
                all()
            for taxa in tax_query:
                tax_ids.add(taxa[0])

    for species in tqdm(genera, desc="Retrieving IDs of species from selected species"):
        with Session(bind=connection) as session:
            tax_query = session.query(Taxonomy.taxonomy_id).\
                filter(Taxonomy.species.like(f'{species}%')).\
                all()
            for taxa in tax_query:
                tax_ids.add(taxa[0])

    for strain in tqdm(strains, desc="Retrieving IDs of species from selected strains"):
        with Session(bind=connection) as session:
            tax_query = session.query(Taxonomy.taxonomy_id).\
                filter(Taxonomy.species == strain).\
                all()
            for taxa in tax_query:
                tax_ids.add(taxa[0])
                
    if len(tax_ids) == 0:
        logger.error(
            "Retrieve NO taxonomy objects matching the provided kingdom and tax filters\n"
            "Therefore, retrieved NO proteins matching the provided criteria.\n"
            "Check the database contains the selected kingdoms, genera, species and strains\n"
            "Terminating program"
        )
        sys.exit(1)
    
    filtered_gbk = set()
    for gbk in initially_selected_gbk:
        if gbk[1].taxonomy_id in tax_ids:
            filtered_gbk.add(gbk[0])
            
    return filtered_gbk


def apply_ec_filters(
    current_gbk_objs,
    ec_filters,
    connection,
):
    """Apply EC number filter to the retrieved Genbank records.
    
    :param current_gbk_objs: list of db Genbank objs retrieved from the db
    :param ec_filters: set of EC numbers to limit the retrieval of data to
    :param connection: open sqlaclchemy connection for an SQLite db
    
    Return set of db Genbank objects.
    """
    logger = logging.getLogger(__name__)
    
    # retrieve the Genbank db IDs for all GenBank accessions associated with each user selected EC num in the db
    gbk_ids = set()

    for ec in tqdm(ec_filters, desc="Retrieving selected EC numbers from db"):
        with Session(bind=connection) as session:
            gbk_query = session.query(Genbank.genbank_id).\
                join(Ec, Genbank.ecs).\
                filter(Ec.ec_number == ec).\
                all()
        for gbk_id in gbk_query:
            gbk_ids.add(gbk_id)

    if len(gbk_ids) == 0:
        logger.error(
        "Retrieved NO proteins matching the provided EC numbers\n"
        "Therefore, not retrievign UniProt data for any proteins.\n"
        "Check the local CAZyme db contains the EC numbers provided\n"
        "Terminating program"
    )
        sys.exit(1)
    
    filtered_gbks = set()
    for gbk in current_gbk_objs:
        if gbk.genbank_id in gbk_ids:
            filtered_gbks.add(gbk)
    
    return filtered_gbks

