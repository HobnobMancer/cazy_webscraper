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
"""Funcs for retrieving GenBank accessions of interest, matching user filter criteria"""


import logging
import sys

from sqlalchemy import select
from sqlalchemy.orm import aliased
from tqdm import tqdm

from cazy_webscraper.sql.sql_orm import (
    CazyFamily,
    Ec,
    Genbank,
    Kingdom,
    Taxonomy,
    Session,
)


CLASS_ABBREVIATIONS = {
    'Glycoside Hydrolases (GHs)': 'GH',
    'GlycosylTransferases (GTs)': 'GT',
    'Polysaccharide Lyases (PLs)': 'PL',
    'Carbohydrate Esterases (CEs)': 'CE',
    'Auxiliary Activities (AAs)': 'AA',
    'Carbohydrate-Binding Modules (CBMs)': 'CBM',
}


def get_ids(genbank_accessions, connection, cache_dir):
    """Get the local CAZyme database IDs for the list of provided GenBank accessions.
    
    :param genbank_accessions: set of GenBank accessions
    :param connection: open sqlalchemy engine connection
    :param cache_dir: path to cache directory
    
    Return dict, keyed by GenBank accession and valued by database record ID.
    """
    cache_path = cache_dir / "seqs_with_no_db_id"
    logger = logging.getLogger(__name__)
    gbk_dict = {}

    for accession in tqdm(genbank_accessions, desc="Getting local db record IDs"):
        with Session(bind=connection) as session:
            gbk_query = session.query(Genbank).\
                filter(Genbank.genbank_accession == accession).\
                first()
        
        try:
            gbk_dict[accession] = gbk_query.genbank_id
        except AttributeError:
            logger.error(
                f"Could not retrieve record with accessions {accession}\n"
                "from the local CAZyme database.\n"
                "Not adding this protein to the local CAZyme database."
                "Logging accessoin in:\n"
                f"{cache_path}"
            )
            with open(cache_path, "a") as fh:
                fh.write(f"{accession}\n")

    return gbk_dict
    

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
            filtered_gbk_accessions,
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
    logger = logging.getLogger(__name__)

    initially_selected_gbk = []

    if len(class_filters) == 0 and len(family_filters) == 0:
        logger.warning("No class or family filters applied")
        # could retrieve all GenBank accessions
        with Session(bind=connection) as session:
            gbk_query = session.query(Genbank, Taxonomy, Kingdom).\
                join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                join(Genbank, (Genbank.taxonomy_id == Taxonomy.taxonomy_id)).\
                join(CazyFamily, Genbank.families).\
                all()

            initially_selected_gbk = gbk_query
        
        return initially_selected_gbk

    if len(class_filters) != 0:
        logger.warning("Applying CAZy class filter(s)")
    for cazy_class in tqdm(class_filters, desc="Retrieving GenBank accessions for selected CAZy classes"):
        if cazy_class not in list(CLASS_ABBREVIATIONS.values()):
            class_abbrev = CLASS_ABBREVIATIONS[cazy_class]
        else:
            class_abbrev = cazy_class
        
        logger.warning(f"Retrieving CAZymes for CAZy class {cazy_class}")

        # perform a subquery to retrieve all CAZy families in the CAZy class
        inner_stmt = select(CazyFamily.family).where(CazyFamily.family.like(f'{class_abbrev}%'))
        subq = inner_stmt.subquery()
        aliased_families = aliased(CazyFamily, subq)
        stmt = select(aliased_families)

        # perform query to retrieve proteins in the CAZy families
        with Session(bind=connection) as session:
            gbk_query = session.query(Genbank, Taxonomy, Kingdom).\
                join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                join(Genbank, (Genbank.taxonomy_id == Taxonomy.taxonomy_id)).\
                join(CazyFamily, Genbank.families).\
                filter(CazyFamily.family.in_(stmt)).\
                all()

        initially_selected_gbk += gbk_query

    if len(family_filters) != 0:
        logger.warning("Applying CAZy family filter(s)")
    for cazy_family in tqdm(family_filters, desc="Retrieving GenBank accessions for selected CAZy families"):

        logger.warning(f"Retrieving CAZymes for CAZy family {cazy_family}")

        inner_stmt = select(CazyFamily.family).where(CazyFamily.family == cazy_family)
        subq = inner_stmt.subquery()
        aliased_families = aliased(CazyFamily, subq)
        stmt = select(aliased_families)

        if cazy_family.find('_') != -1:  # subfamily
            with Session(bind=connection) as session:
                gbk_query = session.query(Genbank, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Genbank, (Genbank.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(CazyFamily, Genbank.families).\
                    filter(CazyFamily.subfamily.in_(stmt)).\
                    all()
        
        else:
            with Session(bind=connection) as session:
                gbk_query = session.query(Genbank, Taxonomy, Kingdom).\
                    join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                    join(Genbank, (Genbank.taxonomy_id == Taxonomy.taxonomy_id)).\
                    join(CazyFamily, Genbank.families).\
                    filter(CazyFamily.family.in_(stmt)).\
                    all()

        initially_selected_gbk += gbk_query
        
    return list(set(initially_selected_gbk))


def apply_tax_filters(
    initally_selected_records,
    taxonomy_filters,
    kingdom_filters,
):
    """Filter retrieved GenBank accessions by taxonomy filters.
    
    :param initally_selected_records: list of db objs retrieved from the db
        including a Genbank, Taxonomy and Kingdom record
    :param taxonomy_filters: dict of taxonom filters to limit the retrieval of data to
    :param kingdom_filters: set of tax kingdoms to limit the retrieval of data to
    :param connection: open sqlaclchemy connection for an SQLite db
    
    Return set of db Genbank objs
    """
    logger = logging.getLogger(__name__)
    
    if len(taxonomy_filters['genera']) == 0 and \
        len(taxonomy_filters['species']) == 0 and \
        len(taxonomy_filters['strains']) == 0 and \
        len(kingdom_filters) == 0:
        logger.warning("Applying no taxonomic filters")
        gbks = [obj[0] for obj in initally_selected_records]
        return set(gbks)
 
    tax_filtered_gbks = set()  # Set of Genbank records from the local database

    if len(kingdom_filters) == 0:
        logger.warning("Not applying kingdom filter(s)")
    for kingdom in tqdm(kingdom_filters, desc="Applying kingdom filter(s)"):
        for obj in initally_selected_records:
            if kingdom == obj[2].kingdom:
                tax_filtered_gbks.add(obj[0])

    if len(taxonomy_filters['genera']) == 0:
        logger.warning("Npt applying genera filters")
    for genus in taxonomy_filters['genera']:
        for obj in initally_selected_records:
            if genus == obj[1].genus:
                tax_filtered_gbks.add(obj[0])
                
    if len(taxonomy_filters['species']) == 0:
        logger.warning("Not applying species filters")
    for species in taxonomy_filters['species']:
        for obj in initally_selected_records:
            if species == obj[1].species.split(" "):
                tax_filtered_gbks.add(obj[0])

    if len(taxonomy_filters['strains']) == 0:
        logger.warning("Not applying strains filters")
    for strain in taxonomy_filters['strains']:
        for obj in initally_selected_records:
            if strain == obj[1].species:
                tax_filtered_gbks.add(obj[0])

    return tax_filtered_gbks


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

    for gbk_record in tqdm(current_gbk_objs, desc="Checking gbk records against EC filters"):
        if (gbk_record.genbank_id,) in ec_gbk_ids:
            ec_filtered_gbks.add(gbk_record)
        
    return ec_filtered_gbks

