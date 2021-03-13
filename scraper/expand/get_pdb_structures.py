#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
"""Retrieve PDB structures from RSCB PDB and write to disk"""


import logging
import os
import sys
import time

import pandas as pd

from datetime import datetime
from typing import List, Optional

from Bio.PDB import PDBList
from tqdm import tqdm

from scraper.sql.sql_orm import (
    Cazyme,
    CazyFamily,
    Kingdom,
    Pdb,
    Taxonomy,
    get_db_session,
)
from scraper.utilities import config_logger, file_io, parse_configuration
from scraper.utilities.parsers import build_pdb_structures_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up programme and initate run."""
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # parse cmd-line arguments
    if argv is None:
        parser = build_pdb_structures_parser()
        args = parser.parse_args()
    else:
        parser = build_pdb_structures_parser(argv)
        args = parser.parse_args()

    # build logger
    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)

    # create session to local database
    if os.path.isfile(args.database) is False:
        logger.error(
            "Could not find local CAZy database. Check path is correct. Terminating programme."
        )
        sys.exit(1)
    session = get_db_session(args)

    # create output directory
    if args.outdir is None:
        # save structure files to the cwd
        outdir = os.getcwd()
    else:
        outdir = args.outdir
    file_io.make_output_directory(outdir, args.force, args.nodelete)

    # check if any classes or families were specified to retrieve the sequences only for them
    file_io_path = file_io.__file__
    config_dict, taxonomy_filters, kingdoms = parse_configuration.get_configuration(
        file_io_path,
        args,
    )

    # retrieve protein structures from PDB

    # retrieve sequences for all CAZymes
    if config_dict is None:
        get_every_cazymes_structures(outdir, taxonomy_filters, kingdoms, session, args)

    # retrieve sequences for specific CAZy classes and/or families
    else:
        get_structures_for_specific_cazymes(
            outdir,
            config_dict,
            taxonomy_filters,
            kingdoms,
            session,
            args,
        )

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    end_time = pd.to_datetime(start_time)
    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    logger.info(
        "Finished dowloading protein structures from PDB."
        "Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}"
    )

    print(
        "=====================cazy_webscraper-expand-pdb_structures=====================\n"
        "Finished dowloading protein structures from PDB."
        "Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}\n"
    )


def get_every_cazymes_structures(outdir, taxonomy_filters, kingdoms, session, args):
    """Get PDB accessions in the db, and coordinate downloading the structure from pdb.

    :param outdir: path to output directory
    :param taxonomy_filters: set of genera, species, and strains to retrieve structures for
    :param kingdoms: set of taxonomy Kingdoms to retrieve structures for
    :param session: open SQLite db session
    :param args: cmd-line argument parser

    Return nothing.
    """
    # retrieve PDB accessions
    pdb_query = session.query(Pdb, Cazyme, Taxonomy, Kingdom).\
        join(Cazyme.pdbs).\
        join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
        join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
        all()

    if taxonomy_filters is None:
        for query_result in pdb_query:
            pdb_accession = query_result[0].pdb_accession
            pdb_accession = pdb_accession[:pdb_accession.find("[")]
            download_pdb_structures(pdb_accession, outdir, args)

    else:
        for query_result in pdb_query:
            source_organism = query_result[-2].genus + query_result[-2].species
            if any(filter in source_organism for filter in taxonomy_filters):
                pdb_accession = query_result[0].pdb_accession
                pdb_accession = pdb_accession[:pdb_accession.find("[")]
                download_pdb_structures(pdb_accession, outdir, args)
            elif query_result[-1].kingdom in kingdoms:
                pdb_accession = query_result[0].pdb_accession
                pdb_accession = pdb_accession[:pdb_accession.find("[")]
                download_pdb_structures(pdb_accession, outdir, args)

    return


def get_structures_for_specific_cazymes(outdir, config_dict, taxonomy_filters, kingdoms, session, args):
    """Retrieve primary PDB accessions for CAZymes meeting criteria in the config_dict.

    :param outdir: path to output directory
    :param config_dict: dict, defines CAZy classes and families to retrieve accessions from
    :param taxonomy_filters: set of genera, species, and strains to retrieve structures for
    :param kingdoms: set of taxonomy Kingdoms to retrieve structures for
    :param session: open SQLite db session
    :param args: cmd-line argument parser

    Return nothing.
    """
    # start with the classes
    if len(config_dict["classes"]) != 0:
        # create a dictionary to convert full class name to abbreviation
        cazy_classes = config_dict["classes"]

        for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):
            # retrieve class name abbreviation
            cazy_class = cazy_class[((cazy_class.find("(")) + 1):((cazy_class.find(")")) - 1)]

            # retrieve the CAZymes from the specified class
            class_subquery = session.query(Cazyme.cazyme_id).\
                join(CazyFamily, Cazyme.families).\
                filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
                subquery()

            # Retrieve PDB accessions for the selected CAZymes
            pdb_query = session.query(Pdb, Cazyme, Taxonomy, Kingdom).\
                join(Cazyme.pdbs).\
                join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                filter(Cazyme.cazyme_id.in_(class_subquery)).all()

            if taxonomy_filters is None:
                for query_result in pdb_query:
                    pdb_accession = query_result[0].pdb_accession
                    pdb_accession = pdb_accession[:pdb_accession.find("[")]
                    download_pdb_structures(pdb_accession, outdir, args)

            else:
                for query_result in pdb_query:
                    source_organism = query_result[-2].genus + query_result[-2].species
                    if any(filter in source_organism for filter in taxonomy_filters):
                        pdb_accession = query_result[0].pdb_accession
                        pdb_accession = pdb_accession[:pdb_accession.find("[")]
                        download_pdb_structures(pdb_accession, outdir, args)
                    elif query_result[-1].kingdom in kingdoms:
                        pdb_accession = query_result[0].pdb_accession
                        pdb_accession = pdb_accession[:pdb_accession.find("[")]
                        download_pdb_structures(pdb_accession, outdir, args)

    # retrieve protein structure for specified families
    for key in config_dict:
        if key == "classes":
            continue
        if config_dict[key] is None:
            continue

        for family in tqdm(config_dict[key], desc=f"Parsing families in {key}"):
            # Select the CAZymes under the specified (sub)family
            if family.find("_") != -1:  # subfamily
                # Retrieve GenBank accessions catalogued under the subfamily
                family_subquery = session.query(Cazyme.cazyme_id).\
                    join(CazyFamily, Cazyme.families).\
                    filter(CazyFamily.subfamily == family).\
                    subquery()

            else:  # family
                # Retrieve GenBank accessions catalogued under the family
                family_subquery = session.query(Cazyme.cazyme_id).\
                    join(CazyFamily, Cazyme.families).\
                    filter(CazyFamily.family == family).\
                    subquery()

            # Retrieve PDB accessions of the selected CAZymes
            pdb_query = session.query(Pdb, Cazyme, Taxonomy, Kingdom).\
                join(Cazyme.pdbs).\
                join(Cazyme, (Cazyme.taxonomy_id == Taxonomy.taxonomy_id)).\
                join(Taxonomy, (Taxonomy.kingdom_id == Kingdom.kingdom_id)).\
                filter(Cazyme.cazyme_id.in_(family_subquery)).all()

            if taxonomy_filters is None:
                for query_result in pdb_query:
                    pdb_accession = query_result[0].pdb_accession
                    pdb_accession = pdb_accession[:pdb_accession.find("[")]
                    download_pdb_structures(pdb_accession, outdir, args)

            else:
                for query_result in pdb_query:
                    source_organism = query_result[-2].genus + query_result[-2].species
                    if any(filter in source_organism for filter in taxonomy_filters):
                        pdb_accession = query_result[0].pdb_accession
                        pdb_accession = pdb_accession[:pdb_accession.find("[")]
                        download_pdb_structures(pdb_accession, outdir, args)
                    elif query_result[-1].kingdom in kingdoms:
                        pdb_accession = query_result[0].pdb_accession
                        pdb_accession = pdb_accession[:pdb_accession.find("[")]
                        download_pdb_structures(pdb_accession, outdir, args)

    return


def download_pdb_structures(pdb_accession, outdir, args):
    """Download protein structure from the RSCB PDB database

    :param pdb_accession: str, accession of record in the PDB database
    :param outdir: path to output directory
    :param args: cmd-line args parser

    Return nothing.
    """
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(f"{pdb_accession}", file_format=args.pdb, pdir=args.outdir)
    time.sleep(2)  # to prevent bombarding the system

    return


if __name__ == "__main__":
    main()
