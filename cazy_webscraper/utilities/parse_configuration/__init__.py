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
"""Module for handling input and output files and directories."""


import json
import logging
import re
import sys
import yaml

from cazy_webscraper.utilities.parse_configuration.cazy_class_synonym_dict import cazy_synonym_dict


def parse_configuration(args):
    """Parse configuration data, and retrieve user specified CAZy classes and families.

    If no user defined configuration is given the default behaviour to scrape the entirity of CAZy
    is enabled.

    :param args: parser arguments

    Return:
    - list of classes not to scrape
    - dict of families to scrape
    - dict of class synonoms
    """
    logger = logging.getLogger(__name__)

    # Get dictionary of accepted CAZy class synonyms
    cazy_class_synonym_dict = get_cazy_class_synonym_dict(args)

    # create dictionary to store families and classes to be scraped
    config_dict = {
        'classes': set(),
        'Glycoside Hydrolases (GHs)': set(),
        'GlycosylTransferases (GTs)': set(),
        'Polysaccharide Lyases (PLs)': set(),
        'Carbohydrate Esterases (CEs)': set(),
        'Auxiliary Activities (AAs)': set(),
        'Carbohydrate-Binding Modules (CBMs)': set(),
    }
    # create dict for storing genera, species and strains scrape is to be restricted to
    taxonomy_filter_dict = {"genera": set(), "species": set(), "strains": set()}

    kingdom_filters = set()

    # Retrieve user specified CAZy classes and families to be scraped at CAZy

    # user passed a YAML configuration file
    if args.config is not None:
        # add configuration data from YAML file yo configuration dictionary
        config_dict, taxonomy_filter_dict, kingdom_filters = get_yaml_configuration(
            config_dict,
            cazy_class_synonym_dict,
            taxonomy_filter_dict,
            kingdom_filters,
            args,
        )

    # add configuration from the cmd-line
    config_dict, taxonomy_filter_dict, kingdom_filters = get_cmd_scrape_config(
        config_dict,
        cazy_class_synonym_dict,
        taxonomy_filter_dict,
        kingdom_filters,
        args,
    )

    if len(kingdom_filters) != 0:
        kingdoms = set()
        correct_kingdoms = {"Archaea", "Bacteria", "Eukaryota", "Viruses", "Unclassified"}
        for kingdom in kingdom_filters:
            user_kingdom = f"{kingdom[0].upper()}{kingdom[1:].lower()}"
            if user_kingdom in correct_kingdoms:
                kingdoms.add(f"{kingdom[0].upper()}{kingdom[1:].lower()}")
        kingdom_filters = kingdoms

    # get list of CAZy classes not to scrape
    excluded_classes = get_excluded_classes(config_dict, cazy_class_synonym_dict)

    # convert dict into a single set of filters
    taxonomy_filter_set = get_filter_set(taxonomy_filter_dict)

    # retrieve class abbreivations
    class_filters = set([cazy_class.split(" (")[0] for cazy_class in set(config_dict['classes'])])

    family_filters = set()
    for key in config_dict:
        if key != 'classes':
            for fam in config_dict[key]:
                family_filters.add(fam)

    return (
        excluded_classes,
        config_dict,
        cazy_class_synonym_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        taxonomy_filter_set,
    )


def get_cazy_class_synonym_dict(args):
    """Retrieve dictionary of acccepted CAZy class synonym names.

    :param args: cmd args parser

    Return dictionary of acccepted CAZy class synonym names.
    """
    logger = logging.getLogger(__name__)

    if args.cazy_synonyms is None:
        logger.warning("Using default CAZy class synonyms")
        cazy_class_synonym_dict = cazy_synonym_dict()
        return cazy_class_synonym_dict

    try:
        logger.warning("Using user defined CAZy class synonyms")
        with open(args.cazy_synonyms, "r") as fh:
            cazy_class_synonym_dict = json.load(fh)
    except FileNotFoundError:
        logger.error(
            "Could not open the CAZy synonym dictionary, required for "
            "translating CAZy class abbreviations.\n"
            "Check the file cazy_class_synonym_dictionary.json is located at:\n"
            f"{args.cazy_synonyms}\n"
            "Terminating programme"
        )
        sys.exit(1)

    return cazy_class_synonym_dict


def get_yaml_configuration(
    config_dict,
    cazy_class_synonym_dict,
    taxonomy_filter_dict,
    kingdom_filters,
    args,
):
    """Parse all data from configuration YAML file.

    :param config_dict: dict, store CAZy classes and families to scrape
    :param cazy_class_synonym: dict, acceppted CAZy class name synonyms
    :param taxonomy_filter: dict, genera, species and strains scrape is to be restricted to
    :param kingdom_filters: set of tax kingdoms to restrict the scrape to
    :param args: cmd args parser

    Return dictionary containing CAZy classes and families named in YAML configuration file.
    """
    logger = logging.getLogger(__name__)

    # open configuration file
    try:
        with open(args.config) as fh:
            yaml_config_dict = yaml.full_load(fh)
    except FileNotFoundError:
        logger.error(
            "Could not find configuration file when option was enabled.\n"
            "Make sure path to the configuration file is correct\n"
            "Scrapping will not be performed becuase configuration is wrong.\n"
            "Had looked for the configuration file at:\n"
            f"{args.config}"
            "Terminating program."
        )
        sys.exit(1)

    for key in yaml_config_dict:
        if key == "classes":
            yaml_cazy_classes = get_yaml_cazy_classes(yaml_config_dict, cazy_class_synonym_dict)
            for cazy_class in yaml_cazy_classes:
                config_dict["classes"].add(cazy_class)

        elif key == "genera":
            if yaml_config_dict["genera"] is not None:
                for genus in yaml_config_dict["genera"]:
                    taxonomy_filter_dict["genera"].add(genus)

        elif key == "species":
            if yaml_config_dict["species"] is not None:
                for organism in yaml_config_dict["species"]:
                    taxonomy_filter_dict["species"].add(organism)

        elif key == "strains":
            if yaml_config_dict["strains"] is not None:
                for organism_strain in yaml_config_dict["strains"]:
                    taxonomy_filter_dict["strains"].add(organism_strain)

        elif key == "kingdoms":
            if yaml_config_dict["kingdoms"] is not None:
                for kingdom in yaml_config_dict["kingdoms"]:
                    if kingdom.lower() == 'unclassified':
                        kingdom_filters.add(kingdom.lower())
                    else:
                        new_kingdom = kingdom[0].upper() + kingdom[1:].lower()
                        if new_kingdom == 'Eukaryotes':
                            new_kingdom = 'Eukaryota'
                        kingdom_filters.add(new_kingdom)

        elif key == 'ec':
            continue

        else:  # key is a CAZy class and underneath are CAZy families to scrape
            if yaml_config_dict[key] is not None:  # CAZy families were listed
                for fam in yaml_config_dict[key]:
                    config_dict[key].add(fam)

    return config_dict, taxonomy_filter_dict, kingdom_filters


def get_yaml_cazy_classes(yaml_config_dict, cazy_class_synonym_dict):
    """Retrieve list of CAZy classes listed in the YAML file, in their standardised CAZy name.

    :param config_dict: dictionary of YAML content
    :param cazy_class_synonym_dict: dictionary of accepted CAZy class name synonyms

    Return list of CAZy classes specified by user to be scraped.
    """
    logger = logging.getLogger(__name__)

    try:
        cazy_classes = yaml_config_dict["classes"]
    except (KeyError, TypeError):
        logger.warning(
            (
                "Did not find the 'classes' tag in the configuration files.\n"
                "Did not retrieve any CAZy classes from configuration file."
            )
        )
        cazy_classes = []

    if cazy_classes is None:
        logger.info("Not items retrieved from 'classes' in the config file.")
        cazy_classes = []

    if len(cazy_classes) != 0:
        # standardise CAZy class names
        cazy_classes = parse_user_cazy_classes(
            cazy_classes,
            cazy_class_synonym_dict,
        )

    return cazy_classes


def parse_user_cazy_classes(cazy_classes, cazy_class_synonym_dict):
    """Standardise the CAZy class names listed in configuration file.

    :param cazy_classes: list, list of CAZy classes from configuration file
    :param cazy_class_synonym_dict: dict, keyed by class name, keyed by list of accepted synonoms
    :param list(cazy_class_synonym_dict.keys()): list, list of all CAZy classes

    Return list of CAZy classes listed by user in the configuration file.
    """
    logger = logging.getLogger(__name__)
    logger.info("Standardising names of class listed in configuration file")

    class_names = list(cazy_class_synonym_dict.keys()) + list(cazy_class_synonym_dict.values())
    accepted_class_names = []
    for i in class_names:
        if type(i) == list:
            for j in i:
                accepted_class_names.append(j)
        else:
            accepted_class_names.append(i)

    standardised_class_names = list(cazy_class_synonym_dict.keys())

    selected_classes = []

    for cazy_class in cazy_classes:
        # identify user defined CAZy classes not written in standardised format
        if cazy_class not in standardised_class_names:

            # not written in standardised format
            # check if accepted class name
            if cazy_class not in accepted_class_names:
                logger.warning(
                    (
                        f"'{cazy_class}' could not be standardised.\n"
                        "Please use a synonym in the file_io/cazy_class_synonym_dictionary.json.\n"
                        f"'{cazy_class}' will NOT be scraped."
                    )
                )

            else:  # written in accepted alternative format
                for standardised_class_name in cazy_class_synonym_dict:
                    if cazy_class in cazy_class_synonym_dict[standardised_class_name]:
                        selected_classes.append(standardised_class_name)

        else:  # written in standardised format
            selected_classes.appent(cazy_class)

    return list(set(cazy_classes))


def get_cmd_scrape_config(
    config_dict,
    cazy_class_synonym_dict,
    taxonomy_filter,
    kingdom_filters,
    args,
):
    """Retrieve scraping configuration data from the cmd-line args parser

    :param config_dict: dict of CAZy classes and families to be scrapped
    :param cazy_class_synonym_dict: dict of accepted CAZy class name synonyms
    :param taxonomy_filter: dict of genera, species and strains to restrict the scrape to
    :param kingdom_filters: set of tax kingdoms to restrict the scrape to
    :param args: cmd-line args parser

    Return config_dict and taxonomy_filter (dict)
    """
    if args.classes is not None:
        cmd_cazy_classes = args.classes.strip().split(",")
        cmd_cazy_classes = parse_user_cazy_classes(cmd_cazy_classes, cazy_class_synonym_dict)
        for cazy_class in cmd_cazy_classes:
            config_dict["classes"].add(cazy_class)

    if args.families is not None:
        config_dict = get_cmd_defined_families(config_dict, args)

    if args.genera is not None:
        cmd_genera = (args.genera).split(",")
        for cmd_genus in cmd_genera:
            taxonomy_filter["genera"].add(cmd_genus)

    if args.species is not None:
        cmd_species = (args.species).split(",")
        for species in cmd_species:
            taxonomy_filter["species"].add(species)

    if args.strains is not None:
        cmd_strains = (args.strains).split(",")
        for strain in cmd_strains:
            taxonomy_filter["strains"].add(strain)

    if args.kingdoms is not None:
        cmd_kingdoms = (args.kingdoms).split(",")
        for kingdom in cmd_kingdoms:
            if kingdom.lower() == 'unclassified':
                kingdom_filters.add(kingdom.lower())
            else:
                new_kingdom = kingdom[0].upper() + kingdom[1:].lower()
                kingdom_filters.add(new_kingdom)

    return config_dict, taxonomy_filter, kingdom_filters


def get_cmd_defined_families(config_dict, args):
    """Retrieve CAZy families specified at the cmd-line for scraping

    :param config_dict: dict of CAZy classes and families to scrape
    :param args: cmd-line args parser

    Return config_dict
    """
    logger = logging.getLogger(__name__)

    cmd_families = (args.families).strip().split(",")

    for fam in cmd_families:
        try:
            if fam.startswith("GH"):
                re.match(r"GH\d+", fam, re.IGNORECASE).group()
                config_dict['Glycoside Hydrolases (GHs)'].add(fam)

            elif fam.startswith("GT"):
                re.match(r"GT\d+", fam, re.IGNORECASE).group()
                config_dict['GlycosylTransferases (GTs)'].add(fam)

            elif fam.startswith("PL"):
                re.match(r"PL\d+", fam, re.IGNORECASE).group()
                config_dict['Polysaccharide Lyases (PLs)'].add(fam)

            elif fam.startswith("CE"):
                re.match(r"CE\d+", fam, re.IGNORECASE).group()
                config_dict['Carbohydrate Esterases (CEs)'].add(fam)

            elif fam.startswith("AA"):
                re.match(r"AA\d+", fam, re.IGNORECASE).group()
                config_dict['Auxiliary Activities (AAs)'].add(fam)

            elif fam.startswith("CBM"):
                re.match(r"CBM\d+", fam, re.IGNORECASE).group()
                config_dict['Carbohydrate-Binding Modules (CBMs)'].add(fam)

            else:
                logger.warning(
                    f"Do recognise class prefix for family specified at cmd line: {fam}\n"
                    "This family will not be scraped."
                )

        except AttributeError:
            logger.warning(
                f"Invalid family specified at cmd line: {fam}\n"
                "This family will not be scraped."
            )

    return config_dict


def get_excluded_classes(config_dict, cazy_class_synonym_dict):
    """Define the CAZy classes that will not be scraped.

    This includes classes for which not Families have been specified for scraping.

    :param config_dict: configuration dict defining classes and families to be scraped
    :param cazy_class_synonym_dict: dict, accepted CAZy classes synonyms

    Return list of CAZy classes not to be scraped.
    """
    cazy_classes_to_scrape = set()

    # retrieve the names of classes for which specific families to be scraped HAVE BEEN named
    for key in config_dict:
        if (key != "classes") and (len(config_dict[key]) != 0):
            # add the class of families to be scraped to the list of CAZy classes to be scraped
            cazy_classes_to_scrape.add(key)
    cazy_classes_to_scrape.union(config_dict["classes"])

    # create list of CAZy classes not to be scraped
    excluded_classes = list(cazy_class_synonym_dict.keys())
    excluded_classes = list(set(excluded_classes).difference(cazy_classes_to_scrape))

    if len(excluded_classes) != 0:
        # change names of classes into format for excluding classes during scrape
        index = 0
        for index in range(len(excluded_classes)):
            excluded_classes[index] = f"<strong>{excluded_classes[index]}</strong>"
    else:
        excluded_classes = None

    return excluded_classes


def convert_empty_sets_to_none(config_dict):
    """Convert empty lists to None type objects in configuration dictionary.

    :param config_dict: dict, CAZy classes and families to be scraped

    Return dictionary with no empty lists.
    """
    for key in config_dict:
        if config_dict[key] is not None:
            if len(config_dict[key]) == 0:
                config_dict[key] = None

    return config_dict


def get_filter_set(taxonomy_filters_dict):
    """Create a set of all taxonomy filters from a dictionary.

    :param taxonomy_filers: dict of taxonomy filters

    Return a set.
    """
    taxonomy_filters = set()

    for key in taxonomy_filters_dict:
        try:
            for tax_filter in taxonomy_filters_dict[key]:
                taxonomy_filters.add(tax_filter)
        except TypeError:
            pass

    return taxonomy_filters


def get_expansion_configuration(args):
    """Get configuration for the Expand module.

    :param args: cmd-line argument parser

    Return configuration dictionary, set of taxonomy filters and set of Taxonomy Kingdoms.
    """
    # Get dictionary of accepted CAZy class synonyms
    cazy_class_synonym_dict = get_cazy_class_synonym_dict(args)

    # create dictionary to store families and classes to be scraped
    config_dict = {
        'classes': set(),
        'Glycoside Hydrolases (GHs)': set(),
        'GlycosylTransferases (GTs)': set(),
        'Polysaccharide Lyases (PLs)': set(),
        'Carbohydrate Esterases (CEs)': set(),
        'Auxiliary Activities (AAs)': set(),
        'Carbohydrate-Binding Modules (CBMs)': set(),
    }

    taxonomy_filter_dict = {"genera": set(), "species": set(), "strains": set()}
    kingdom_filters = set()
    ec_filters = set()

    # retrieve user configuration, defining for which records data should be retrieved
    if args.config is not None:
        config_dict, taxonomy_filter_dict, kingdom_filters = get_yaml_configuration(
            config_dict,
            cazy_class_synonym_dict,
            taxonomy_filter_dict,
            kingdom_filters,
            args,
        )

    # add configuration from the cmd-line
    config_dict, taxonomy_filter_dict, kingdom_filters = get_cmd_scrape_config(
        config_dict,
        cazy_class_synonym_dict,
        taxonomy_filter_dict,
        kingdom_filters,
        args,
    )

    ec_filters = get_ec_config(ec_filters, args)

    class_filters = set(config_dict['classes'])

    family_filters = set()
    for key in config_dict:
        if key != 'classes':
            for fam in config_dict[key]:
                family_filters.add(fam)

    return (
        config_dict,
        class_filters,
        family_filters,
        kingdom_filters,
        taxonomy_filter_dict,
        ec_filters,
    )


def get_ec_config(ec_filters, args):
    """Parse EC number configuration.

    Standardise the EC numbers missing digits and remove 'EC' prefix.

    :param ec_filters: set of EC numbers
    :param args: cmd line args parser

    Return set of EC numbers.
    """
    logger = logging.getLogger(__name__)

    if args.config is not None:
        # open configuration file
        try:
            with open(args.config) as fh:
                yaml_config_dict = yaml.full_load(fh)
        except FileNotFoundError:
            logger.error(
                "Could not find configuration file when option was enabled.\n"
                "Make sure path to the configuration file is correct\n"
                "Scrapping will not be performed becuase configuration is wrong.\n"
                "Had looked for the configuration file at:\n"
                f"{args.config}"
                "Terminating program."
            )
            sys.exit(1)

        ec_filters = set(yaml_config_dict['ec'])

    if ec_filters is not None:
        try:
            ec_filters = ec_filters.union(set((args.ec_filter).split(",")))
        except AttributeError:
            pass

        ec_filters = [ec.replace("EC", "") for ec in ec_filters]
        ec_filters = [ec.replace("ec", "") for ec in ec_filters]
        ec_filters = [ec.replace("*", "-") for ec in ec_filters]

        ec_filters = set(ec_filters)

    return ec_filters
