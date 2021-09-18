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
"""Module for handling input and output files and directories."""


import json
import logging
import re
import sys
import yaml

from pathlib import Path


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
    taxonomy_filter = {"genera": set(), "species": set(), "strains": set()}

    # Retrieve user specified CAZy classes and families to be scraped at CAZy

    # user passed a YAML configuration file
    if args.config is not None:
        # add configuration data from YAML file yo configuration dictionary
        config_dict = get_yaml_configuration(config_dict, cazy_class_synonym_dict, taxonomy_filter, args)

        # add configuration from the cmd-line
        config_dict,taxonomy_filter = get_cmd_scrape_config(
            config_dict,
            cazy_class_synonym_dict,
            taxonomy_filter, 
            args,
        )

    else:  
        config_dict, taxonomy_filter = get_cmd_scrape_config(
            config_dict,
            cazy_class_synonym_dict,
            taxonomy_filter, 
            args,
        )
    
    # get list of CAZy classes not to scrape
    excluded_classes = get_excluded_classes(config_dict, cazy_class_synonym_dict)

    # convert empty sets to None type objects
    config_dict = convert_empty_sets_to_none(config_dict)

    # convert dict into a single set of filters
    taxonomy_filter = get_filter_set(taxonomy_filter)

    return excluded_classes, config_dict, cazy_class_synonym_dict, taxonomy_filter


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
            "Could not open the CAZy synonym dictionary, required for translating CAZy class abbreviations.\n"
            "Check the file cazy_class_synonym_dictionary.json is located at:\n"
            f"{args.cazy_synonyms}\n"
            "Terminating programme"
        )
        sys.exit(1)

    return cazy_class_synonym_dict


def cazy_synonym_dict():
    """Create a dictionary of accepted synonms for CAZy classes."""
    cazy_class_synonym_dict = {
        "Glycoside Hydrolases (GHs)": ["Glycoside-Hydrolases", "Glycoside-Hydrolases", "Glycoside_Hydrolases", "GlycosideHydrolases", "GLYCOSIDE-HYDROLASES", "GLYCOSIDE-HYDROLASES", "GLYCOSIDE_HYDROLASES", "GLYCOSIDEHYDROLASES", "glycoside-hydrolases", "glycoside-hydrolases", "glycoside_hydrolases", "glycosidehydrolases", "GH", "gh"],
        "GlycosylTransferases (GTs)": ["Glycosyl-Transferases", "GlycosylTransferases", "Glycosyl_Transferases", "Glycosyl Transferases", "GLYCOSYL-TRANSFERASES", "GLYCOSYLTRANSFERASES", "GLYCOSYL_TRANSFERASES", "GLYCOSYL TRANSFERASES", "glycosyl-transferases", "glycosyltransferases", "glycosyl_transferases", "glycosyl transferases", "GT", "gt"],
        "Polysaccharide Lyases (PLs)": ["Polysaccharide Lyases", "Polysaccharide-Lyases", "Polysaccharide_Lyases", "PolysaccharideLyases", "POLYSACCHARIDE LYASES", "POLYSACCHARIDE-LYASES", "POLYSACCHARIDE_LYASES", "POLYSACCHARIDELYASES", "polysaccharide lyases", "polysaccharide-lyases", "polysaccharide_lyases", "polysaccharidelyases", "PL", "pl"],
        "Carbohydrate Esterases (CEs)": ["Carbohydrate Esterases", "Carbohydrate-Esterases", "Carbohydrate_Esterases", "CarbohydrateEsterases", "CARBOHYDRATE ESTERASES", "CARBOHYDRATE-ESTERASES", "CARBOHYDRATE_ESTERASES", "CARBOHYDRATEESTERASES", "carbohydrate esterases", "carbohydrate-esterases", "carbohydrate_esterases", "carbohydrateesterases", "CE", "ce"],
        "Auxiliary Activities (AAs)": ["Auxiliary Activities", "Auxiliary-Activities", "Auxiliary_Activities", "AuxiliaryActivities", "AUXILIARY ACTIVITIES", "AUXILIARY-ACTIVITIES", "AUXILIARY_ACTIVITIES", "AUXILIARYACTIVITIES", "auxiliary activities", "auxiliary-activities", "auxiliary_activities", "auxiliaryactivities", "AA", "aa"],
        "Carbohydrate-Binding Modules (CBMs)": ["Carbohydrate-Binding-Modules", "Carbohydrate_Binding_Modules", "Carbohydrate_Binding Modules", "CarbohydrateBindingModules", "CARBOHYDRATE-BINDING-MODULES", "CARBOHYDRATE_BINDING_MODULES", "CARBOHYDRATE_BINDING MODULES", "CARBOHYDRATEBINDINGMODULES", "carbohydrate-binding-modules", "carbohydrate_binding_modules", "carbohydrate_binding modules", "carbohydratebindingmodules", "CBMs", "CBM", "cbms", "cbm"]
    }
    return cazy_class_synonym_dict


def get_yaml_configuration(config_dict, cazy_class_synonym_dict, taxonomy_filter, args):
    """Parse all data from configuration YAML file.

    :param config_dict: dict, store CAZy classes and families to scrape
    :param cazy_class_synonym: dict, acceppted CAZy class name synonyms
    :param taxonomy_filter: dict, genera, species and strains scrape is to be restricted to
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

    # retrieve CAZy classes defined in the YAML configuration file
    config_dict["classes"] = get_yaml_cazy_classes(
        yaml_config_dict,
        cazy_class_synonym_dict,
    )

    for key in yaml_config_dict:
        if key == "classes":
            yaml_cazy_classes = get_yaml_cazy_classes(yaml_config_dict, cazy_class_synonym_dict)
            for cazy_class in yaml_cazy_classes:
                config_dict["classes"].add(cazy_class)
            
        elif key == "genera":
            if yaml_config_dict["genera"] is not None:
                for genus in taxonomy_filter["genera"]:
                    config_dict["genera"].add(genus)
               
        elif key == "species":
            if yaml_config_dict["genera"] is not None:
                for organism in taxonomy_filter["genera"]:
                    config_dict["genera"].add(organism)

        elif key == "strains":
            if yaml_config_dict["genera"] is not None:
                for organism_strain in taxonomy_filter["genera"]:
                    config_dict["genera"].add(organism_strain)

        else:  # key is a CAZy class and underneath are CAZy families to scrape
            if yaml_config_dict[key] is not None:  # CAZy families were listed
                for fam in yaml_config_dict[key]:
                    config_dict[key].add(fam)

    return config_dict


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

    standardised_class_names = list(cazy_class_synonym_dict.keys())

    index = 0
    for index in range(len(cazy_classes)):
        # identify user defined CAZy classes not written in standardised format
        if cazy_classes[index] not in standardised_class_names:
            for key in cazy_class_synonym_dict:
                if cazy_classes[index] in cazy_class_synonym_dict[key]:  # if in synonyms in cazy_class_synonym_dict
                    cazy_classes[index] = key  # standardise the class name

        # check all names are standardised, remove names that could not be standardised
        if cazy_classes[index] not in standardised_class_names:
            logger.warning(
                (
                    f"'{cazy_classes[index]}' could not be standardised.\n"
                    "Please use a synonym in the file_io/cazy_class_synonym_dictionary.json.\n"
                    f"'{cazy_classes[index]}' will NOT be scraped."
                )
            )
            del cazy_classes[index]

    return cazy_classes


def get_cmd_scrape_config(
    config_dict,
    cazy_class_synonym_dict,
    taxonomy_filter, 
    args,
):
    """Retrieve scraping configuration data from the cmd-line args parser
    
    :param config_dict: dict of CAZy classes and families to be scrapped
    :param cazy_class_synonym_dict: dict of accepted CAZy class name synonyms
    :param taxonomy_filter: dict of genera, species and strains to restrict the scrape to
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
        for cmd_species in cmd_species:
            taxonomy_filter["species"].add(cmd_species)
    
    if args.strains is not None:
        cmd_strains = (args.strains).split(",")
        for cmd_genus in cmd_strains:
            taxonomy_filter["strains"].add(cmd_strains)
    
    return config_dict, taxonomy_filter


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
            if fam.find("GH") != -1:
                re.match(r"GH\d+", fam, re.IGNORECASE).group()
                config_dict['Glycoside Hydrolases (GHs)'].add(fam)

            elif fam.find("GT") != -1:
                re.match(r"GT\d+", fam, re.IGNORECASE).group()
                config_dict['GlycosylTransferases (GTs)'].add(fam)

            elif fam.find("PL") != -1:
                re.match(r"PL\d+", fam, re.IGNORECASE).group()
                config_dict['Polysaccharide Lyases (PLs)'].add(fam)

            elif fam.find("CE") != -1:
                re.match(r"CE\d+", fam, re.IGNORECASE).group()
                config_dict['Carbohydrate Esterases (CEs)'].add(fam)

            elif fam.find("AA") != -1:
                re.match(r"AA\d+", fam, re.IGNORECASE).group()
                config_dict['Auxiliary Activities (AAs)'].add(fam)

            elif fam.find("CBM") != -1:
                re.match(r"CBM\d+", fam, re.IGNORECASE).group()
                config_dict['Carbohydrate-Binding Modules (CBMs)'].ad(fam)

            else:
                logger.warning(
                    f"Invalid family specified at cmd line: {fam}\n"
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
    cazy_classes_to_scrape = []

    # retrieve the names of classes for which specific families to be scraped HAVE BEEN named
    for key in config_dict:
        if (key != "classes") and (key not in config_dict["classes"]) and (len(config_dict[key]) != 0):
            # add the class of families to be scraped to the list of CAZy classes to be scraped
            cazy_classes_to_scrape.append(key)

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
    taxonomy_filters = []

    for key in taxonomy_filters_dict:
        try:
            if len(taxonomy_filters_dict[key]) != 0:
                taxonomy_filters += taxonomy_filters_dict[key]
        except TypeError:
            pass

    if len(taxonomy_filters) == 0:
        taxonomy_filters = None

    else:
        taxonomy_filters = set(taxonomy_filters)

    return taxonomy_filters


################################################################################


def create_streamline_scraping_warning(args):
    """Creating warning message to flag 'streamlined' scraping has been enabled.

    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    streamline = (args.streamline).split(",")
    index = 0
    for index in range(len(streamline)):
        if streamline[index] == "genbank":
            streamline[index] = "GenBank primary and non-primary accessions"
        if streamline[index] == "ec":
            streamline[index] = "EC numbers"
        if streamline[index] == "uniprot":
            streamline[index] = "UniProt primary and non-primary accessions"
        if streamline[index] == "pdb":
            streamline[index] = "PDB accessions"

    streamline = ',\n'.join(streamline)
    logger.warning(
        "Enabled 'streamlined' scraping. "
        "Presuming the following data is identical for each family table a protein appears in\n"
        f"{streamline}"
    )


def get_configuration(args):
    """Get configuration for the Expand module.

    :param file_io_path: Path to file_io module
    :param args: cmd-line argument parser

    Return configuration dictionary, set of taxonomy filters and set of Taxonomy Kingdoms.
    """
    # retrieve inital parsing of configuration data
    (
        excluded_classes, config_dict, cazy_class_synonym_dict, taxonomy_filters_dict, kingdoms, ec_filter,
    ) = parse_configuration(args)
    # excluded_classes and cazy_class_synonym_dict are used in the crawler module but are not needed for the
    # the expand module

    taxonomy_filters = []

    for key in taxonomy_filters_dict:
        try:
            if len(taxonomy_filters_dict[key]) != 0:
                taxonomy_filters += taxonomy_filters_dict[key]
        except (TypeError, KeyError):
            pass

    if len(taxonomy_filters) == 0:
        taxonomy_filters = None

    else:
        taxonomy_filters = set(taxonomy_filters)

    if kingdoms == "all":
        kingdoms = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses', 'Unclassified']

    kingdoms = set(kingdoms)

    return config_dict, taxonomy_filters, kingdoms
