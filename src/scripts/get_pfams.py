#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2026
# (c) University of Strathclyde 2026
# (c) James Hutton Institute 2026
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# ehobbs@ebi.ac.uk
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
"""Coordinate retrieving Pfam domain annotations from InterPro."""


import logging

from argparse import Namespace

from saintBioutils.utilities.file_io import make_output_directory

from src.databases.interpro import get_pfam_data
from src.sql import sql_orm
from src.sql.interface.add_data.scrape_log import log_scrape_in_db
from src.sql.interface.connect import connect_existing_db
from src.sql.interface.get_data.get_proteins import get_ncbi_acc_for_uniprot_acc, get_uniprot_accessions
from src.utilities.parse_configuration import get_expansion_configuration
from src.utilities.parse_configuration.accession_files import get_acc_from_file


logger = logging.getLogger(__name__)


def main(args: Namespace, time_stamp: str, start_time):
    connection, logger_name, cache_dir = connect_existing_db(args, time_stamp, start_time)

    cache_dir = args.cache_dir if args.cache_dir else (cache_dir / "pfam_retrieval")
    make_output_directory(cache_dir, args.force, args.nodelete_cache)

    (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    ) = get_expansion_configuration(args)

    with sql_orm.Session(bind=connection) as session:
        retrieved_data = "Pfam domain annotations from InterPro"
        log_scrape_in_db(
            time_stamp,
            config_dict,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            'InterPro',
            retrieved_data,
            session,
            args,
        )

    if args.genbank_accessions or args.uniprot_accessions:
        protein_accessions = set()
        if args.genbank_accessions:
            protein_accessions.update(
                get_acc_from_file(args.genbank_accessions, args.database)
            )
        if args.uniprot_accessions:
            uniprot_accessions = get_acc_from_file(
                args.uniprot_accessions, args.database, table="UniProt"
            )
            protein_accessions.update(
                get_ncbi_acc_for_uniprot_acc(uniprot_accessions, args.database)
            )

        if not protein_accessions:
            logger.warning("No proteins in the local db matched the accessions provided")
            return "get_pfams"

        uniprot_pairs = get_uniprot_accessions(
            class_filters,
            family_filters,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            args.database,
            protein_accessions=sorted(protein_accessions),
        )
    elif args.update:
        # retrieve for every matching protein, including ones already carrying Pfam matches
        uniprot_pairs = get_uniprot_accessions(
            class_filters,
            family_filters,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            args.database,
        )
    else:
        # only retrieve for proteins with no Pfam matches yet
        uniprot_pairs = get_uniprot_accessions(
            class_filters,
            family_filters,
            kingdom_filters,
            taxonomy_filter_dict,
            ec_filters,
            args.database,
            additional_filter="P.protein_id NOT IN (SELECT protein_id FROM Proteins_Pfams)",
        )

    if not uniprot_pairs:
        logger.warning(
            "No proteins with UniProt accessions matched the criteria provided.\n"
            "Retrieving no Pfam data from InterPro.\n"
            "Note: get_pfams looks up Pfam matches by UniProt accession, so a protein must "
            "already have UniProt data before it can be included - run get_uniprot_data first."
        )
        return "get_pfams"

    uniprot_acc_to_protein_id = {uniprot_acc: protein_id for protein_id, uniprot_acc in uniprot_pairs}

    logger.warning("Retrieving Pfam data for %d UniProt accessions", len(uniprot_acc_to_protein_id))
    get_pfam_data(uniprot_acc_to_protein_id, cache_dir, time_stamp, args)

    return "get_pfams"
