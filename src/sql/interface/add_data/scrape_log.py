#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
# (c) James Hutton Institute 2024
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


import argparse
import logging

from src.sql import sql_interface, sql_orm
from src.utilities import termcolour


logger = logging.getLogger(__name__)


def add_main_scrape_message(
    kingdom_filters: set[str],
    taxonomy_filters: set[str],
    taxonomy_filter_dict: dict,
    time_stamp: str,
    config_dict: dict,
    args: argparse.Namespace,
    connection
) -> None:
    """add information of scraping CAZy to the local CAZyme database"""
    scrape_config_message = (
        "Configuration:\n"
        f"Classes to scrape: {config_dict['classes']}\n"
        f"GH fams to scrape: {config_dict['Glycoside Hydrolases (GHs)']}\n"
        f"GT fams to scrape: {config_dict['GlycosylTransferases (GTs)']}\n"
        f"PL fams to scrape: {config_dict['Polysaccharide Lyases (PLs)']}\n"
        f"CE fams to scrape: {config_dict['Carbohydrate Esterases (CEs)']}\n"
        f"AA fams to scrape: {config_dict['Auxiliary Activities (AAs)']}\n"
        f"CBM fams to scrape: {config_dict['Carbohydrate-Binding Modules (CBMs)']}\n"
        f"Scraping subfamilies: {args.subfamilies}"
    )
    scrape_config_message += "\nTaxonomy filters applied." if len(taxonomy_filters) != 0 else ""
    scrape_config_message += f"\nScraping only tax kingdoms: {kingdom_filters}" if len(kingdom_filters) < 5 else ""
    
    logger.info(termcolour(scrape_config_message, "cyan"))

    logger.info("Adding log of scrape to the local CAZyme database")
    with sql_orm.Session(bind=connection) as session:
        sql_interface.log_scrape_in_db(
            time_stamp,
            config_dict,
            kingdom_filters,
            taxonomy_filter_dict,
            set(),  # ec_filters not applied when scraping CAZy
            'CAZy',
            'CAZy annotations',
            session,
            args,
        )


def log_scrape_in_db(
    time_stamp,
    config_dict,
    taxonomy_filters,
    kingdoms,
    ec_filter,
    db,
    retrieved_annotations,
    session,
    args,
):
    """Add a log of scraping CAZy to the local database.

    :param time_stamp: str, date and time cazy_webscraper was invoked
    :param config_dict: dict of CAZy classes and families to be scraped
    :param taxonomy_filters: dict of genera, species and strains to restrict the scrape to
    :param kingdoms: list of taxonomy Kingdoms to restrict scrape to
    :param ec_filter: set of EC numbers to restrict retrieval of data to
    :param db: str, name of the external database from which data is retrieved
    :param retrieved_annotations: str, types of annotations retrieved (e.g. UniProt accessions)
    :param session: open SQL database session
    :param args: cmd arguments

    Return nothing."""
    logger = logging.getLogger(__name__)

    logger.info("Adding log of scrape to db")

    date = time_stamp.split("_")[0]
    time = time_stamp.split("_")[1]

    new_log = sql_orm.Log(
        date=date,
        time=time,
        database=db,
        retrieved_annotations=retrieved_annotations,
    )

    classes = []

    if config_dict is not None:
        # get classes that user named to be scraped
        try:
            classes = config_dict["classes"]
            # E.G. {'classes': ['Polysaccharide Lyases (PLs)', 'Carbohydrate Esterases (CEs)']}
            if classes is not None:
                classes = ""
                for cazy_class in config_dict['classes']:
                    if len(classes) == 0:
                        classes = cazy_class
                    else:
                        classes += f", {cazy_class}"

                new_log.classes = classes
        except KeyError:
            pass

        if len(classes) != 0:
            new_log.classes = classes

        # create a list of families instructed to be scraped
        families = ""
        for key in config_dict:
            if (key != "classes") and (config_dict[key] is not None):
                for fam in config_dict[key]:
                    if len(families) == 0:
                        families = fam
                    else:
                        families += f", {fam}"

        if len(families) != 0:
            new_log.families = families

    # get taxonomy filters defined by user, and separate into genera, species and strains
    try:
        genera = ""
        if len(taxonomy_filters["genera"]) != 0:
            for genus in taxonomy_filters["genera"]:
                if len(genera) == 0:
                    genera = genus
                else:
                    genera += f", {genus}"
        if len(genera) != 0:
            new_log.genera = genera
    except (TypeError, KeyError):
        pass

    try:
        species = ""
        if len(taxonomy_filters["species"]) != 0:
            for organism in taxonomy_filters["species"]:
                if len(species) == 0:
                    species = organism
                else:
                    species += f", {organism}"

        if len(species) != 0:
            new_log.species = species
    except (TypeError, KeyError):
        pass

    try:
        strains = ""
        if len(taxonomy_filters["strains"]) != 0:
            for organism in taxonomy_filters["strains"]:
                if len(strains) == 0:
                    strains = organism
                else:
                    strains += f", {organism}"

        if len(strains) != 0:
            new_log.strains = strains
    except (TypeError, KeyError):
        pass

    # get Taxonomy Kingdoms defined by user to be scraped
    if kingdoms is not None:
        kingdoms_str = ""
        for kingdom in kingdoms:
            if len(kingdoms_str) == 0:
                kingdoms_str = kingdom
            else:
                kingdoms_str += f", {kingdom}"

        if len(kingdoms_str) != 0:
            new_log.kingdoms = kingdoms_str
        else:
            new_log.kingdoms = "ALL (Archaea, Bacteria, Eukaryota, Viruses, Unclassified)"
    else:
        new_log.kingdoms = "ALL (Archaea, Bacteria, Eukaryota, Viruses, Unclassified)"

    # retrieve commands from the command line
    cmd_line = ""
    for cmd in [
        [args.classes, " --classes '"],
        [args.families, " --families '"],
        [args.kingdoms, " --kingdoms '"],
        [args.genera, " --genera '"],
        [args.species, " --species '"],
        [args.strains, " --strains '"],
    ]:
        try:
            cmd_line = cmd_line + cmd[1] + cmd[0] + "' "
        except TypeError:
            pass

    if len(ec_filter) != 0:
        cmd_line = cmd_line + "--ec_filter '" + (args.ec_filter) + "'"
        new_log.ec_filter = ','.join(list(ec_filter))

    if len(cmd_line) != 0:
        cmd_line = cmd_line.strip()
        new_log.cmd_line = cmd_line

    session.add(new_log)
    session.commit()

    return
