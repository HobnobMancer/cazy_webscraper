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


def parse_configuration(file_io_path, args):
    """Parse configuration data, and retrieve user specified CAZy classes and families.

    If no user defined configuration is given the default behaviour to scrape the entirity of CAZy
    is enabled.

    :param file_io_path: str, path to directory where file_io is installed
    :param args: parser arguments

    Return list of classes not to scrape, dict of families to scrape, dict of class synonoms,
    a dict of taxonomy filters (genera, species and strains) to restrict the scrape to, list of
    kingdoms to be scraped, a list of EC numbers to restrict the scrape to.
    """
    logger = logging.getLogger(__name__)
    # open config dict
    raw_config_dict = None
    if args.config is not None:
        try:
            with open(args.config, "r") as fh:
                raw_config_dict = yaml.full_load(fh)
        except FileNotFoundError:
            logger.error(
                "Did not find the configuration file. Check the path is correct.\n"
                "Make sure path to the configuration file is correct\n"
                "Scrapping will not be performed becuase configuration is wrong.\n"
                "Had looked for the configuration file at:\n"
                f"{args.config}"
                "Terminating programme"
            )
            sys.exit(1)

    # retrieve taxonomy filters
    taxonomy_filter = get_genera_species_strains(args, raw_config_dict)

    ec_filter = get_ec_filter(args, raw_config_dict)

    # retrieve Kingdoms to scrape
    kingdoms = get_kingdoms(args, raw_config_dict)

    # Get dictionary of accepted CAZy class synonyms
    cazy_dict, std_class_names = get_cazy_dict_std_names(file_io_path)

    # Retrieve user specified CAZy classes and families to be scraped at CAZy

    # create dictionary to store families and classes to be scraped
    config_dict = {
        'classes': [],
        'Glycoside Hydrolases (GHs)': [],
        'GlycosylTransferases (GTs)': [],
        'Polysaccharide Lyases (PLs)': [],
        'Carbohydrate Esterases (CEs)': [],
        'Auxiliary Activities (AAs)': [],
        'Carbohydrate-Binding Modules (CBMs)': [],
    }

    # user passed a YAML configuration file
    if args.config is not None:
        # add configuration data from YAML file yo configuration dictionary
        config_dict = get_yaml_configuration(config_dict, cazy_dict, std_class_names, args)

        if (args.classes is None) and (args.families is None):  # no cmd-line configuration
            # get list of CAZy classes not to scrape
            excluded_classes = get_excluded_classes(std_class_names, config_dict, cazy_dict)

            # convert empty lists to None type objects
            config_dict = convert_lists_to_none(config_dict)

            return excluded_classes, config_dict, cazy_dict, taxonomy_filter, kingdoms, ec_filter

        else:  # get cmd defined configuration
            cmd_config = get_cmd_defined_fams_classes(cazy_dict, std_class_names, args)
            cmd_config = convert_lists_to_none(cmd_config)

            # add cmd defined configuration to config_dict
            # add items from file_config to cmd_config
            for key in cmd_config:
                if cmd_config[key] is not None:
                    for item in cmd_config[key]:
                        if item not in config_dict[key]:  # do not add duplicates
                            config_dict[key].append(item)

            # get list of CAZy classes that will not be scraped
            excluded_classes = get_excluded_classes(std_class_names, config_dict, cazy_dict)

            # convert empty lists to None type objects
            config_dict = convert_lists_to_none(config_dict)

            return excluded_classes, config_dict, cazy_dict, taxonomy_filter, kingdoms, ec_filter

    else:  # user did not pass config file
        if (args.classes is None) and (args.families is None):
            # No specific families or classes specified for scraping
            return None, None, cazy_dict, taxonomy_filter, kingdoms, ec_filter

        else:  # configuration specified only via the cmd_line
            config_dict = get_cmd_defined_fams_classes(cazy_dict, std_class_names, args)

            # get list of CAZy classes that will not be scraped
            excluded_classes = get_excluded_classes(std_class_names, config_dict, cazy_dict)

            # convert empty lists to None type objects
            config_dict = convert_lists_to_none(config_dict)

            return excluded_classes, config_dict, cazy_dict, taxonomy_filter, kingdoms, ec_filter


def get_genera_species_strains(args, raw_config_dict):
    """Retrieve Genera, Species and Strains to retrict the scrape to.

    :param args: cmd-line arguments parser
    :param raw_config_dict: dictionary of content from YAML config file

    Return set contain user specificed genera, species and strains.
    """
    taxonomy_filter = {"genera": [], "species": [], "strains": []}

    if raw_config_dict is not None:
        for key in ["genera", "species", "strain"]:
            try:
                if raw_config_dict[key] is not None:
                    taxonomy_filter[key] = raw_config_dict[key]
            except KeyError:
                pass

    if args.genera is not None:
        taxonomy_filter["genera"] += (args.genera).split(",")

    if args.species is not None:
        taxonomy_filter["species"] += (args.species).split(",")

    if args.strains is not None:
        taxonomy_filter["strains"] += (args.strains).split(",")

    taxonomy_filter = convert_lists_to_none(taxonomy_filter)

    return taxonomy_filter


def get_ec_filter(args, raw_config_dict):
    """Retrieve EC numbers to retrict the scrape to.

    :param args: cmd-line arguments parser
    :param raw_config_dict: dictionary of content from YAML config file

    Return set containing user specificed EC numbers.
    """
    ecs = []

    if raw_config_dict is not None:
        try:
            if raw_config_dict["ECs"] is not None:
                ecs += raw_config_dict["ECs"]
        except KeyError:
            pass

    if args.ec is not None:
        ecs += (args.ec).split(",")
    
    i = 0
    for i in range(len(ecs)):
        ecs[i] = ecs[i].replace("EC", "")
        ecs[i] = ecs[i].replace("ec", "")

    ec_filter = set(ecs)

    return ec_filter


def get_kingdoms(args, raw_config_dict):
    """Retrieve list of Kingdoms to scrape CAZymes from.

    :param args: cmd args parser
    :param raw_config_dict: dictionary of content from YAML config file

    Return list of Kingdoms to  be scrapped.
    """
    logger = logging.getLogger(__name__)

    kingdoms = []

    if raw_config_dict is not None:
        try:
            if raw_config_dict["kingdoms"] is not None:
                kingdoms += raw_config_dict["kingdoms"]
                yaml_kingdoms = True
            else:
                yaml_kingdoms = False
        except KeyError:
            yaml_kingdoms = False
    else:
        yaml_kingdoms = False

    if args.kingdoms is not None:
        kingdoms += (args.kingdoms).split(",")

    kingdoms = [kingdom.lower() for kingdom in kingdoms]
    kingdoms = list(set(kingdoms))
    cazy_kingdoms = ['archaea', 'bacteria', 'eukaryota', 'viruses', 'unclassified']
    for kd in kingdoms:
        if kd not in cazy_kingdoms:
            logger.warning(
                f"'{kd}' written in configuration, but it not a CAZy listed taxonomy Kingdom.\n"
                f"CAZymes will not be scrapped from the user listed kingdom: {kd}"
            )
            kingdoms.remove(kd)

    if (len(kingdoms) == 0):
        if (args.kingdoms is not None) or (yaml_kingdoms is True):
            # user had specified Kingdoms to be scraped
            logger.warning(
                "None of the Kingdoms listed match the Kingdoms listed in CAZy.\n"
                "The classes in CAZy are: archaea, bacteria, eukaryota, viruses, and unclassified\n"
                "Terminating programme"
            )
            sys.exit(1)

        else:
            # user did not specify any Kingdoms
            kingdoms = 'all'  # scrape CAZymes from all Kingdoms

    return kingdoms


def get_cazy_dict_std_names(file_io_path):
    """Retrieve dictionary of acccepted CAZy class synonym names, and list of offical class.

    :param args: cmd args parser

    Return dictionary.
    """
    logger = logging.getLogger(__name__)

    # build path to the JSON file containing the cazy dict
    dict_path = file_io_path.replace("__init__.py", "cazy_dictionary.json")
    dict_path = Path(dict_path)

    try:
        with open(dict_path, "r") as fh:
            cazy_dict = json.load(fh)
            # create list of standardised CAZy classes
            std_class_names = list(cazy_dict.keys())
    except FileNotFoundError:
        logger.error(
            "Could not open the CAZy synonym dictionary, required for translating CAZy class abbreviations.\n"
            "Check the file cazy_dictionary.json is located at:\n"
            f"{dict_path}\n"
            "Terminating programme"
        )
        sys.exit(1)

    return cazy_dict, std_class_names


def get_yaml_configuration(config_dict, cazy_dict, std_class_names, args):
    """Parse data from configuration YAML file.

    :param config_dict: dict, store CAZy classes and families to scrape
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
        cazy_dict,
        std_class_names,
    )

    for key in yaml_config_dict:
        if (key != "classes") and \
           (key != "genera") and \
           (key != "species") and \
           (key != "strains") and \
           (key != "ECs"):
            if yaml_config_dict[key] is not None:
                for item in yaml_config_dict[key]:
                    if item not in config_dict[key]:  # do not add duplicates
                        config_dict[key].append(item)

    return config_dict


def get_yaml_cazy_classes(yaml_config_dict, cazy_dict, std_class_names):
    """Retrieve list of CAZy classes listed in the YAML file, in their standardised CAZy name.

    :param config_dict: dictionary of YAML content
    :param cazy_dict: dictionary of accepted CAZy class name synonyms

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
        logger.info("Did not retrieve any items under 'classes' in the config file.")
        cazy_classes = []

    if len(cazy_classes) != 0:
        # standardise CAZy class names
        cazy_classes = parse_user_cazy_classes(cazy_classes, cazy_dict, std_class_names)

    return cazy_classes


def parse_user_cazy_classes(cazy_classes, cazy_dict, std_class_names):
    """Standardise the CAZy class names listed in configuration file.

    :param cazy_classes: list, list of CAZy classes from configuration file
    :param cazy_dict: dict, keyed by class name, keyed by list of accepted synonoms
    :param std_class_names: list, list of all CAZy classes

    Return list of CAZy classes listed by user in the configuration file.
    """
    logger = logging.getLogger(__name__)
    logger.info("Standardising names of class listed in configuration file")

    index = 0
    for index in range(len(cazy_classes)):
        # identify user defined CAZy classes not written in standardised format
        if cazy_classes[index] not in std_class_names:
            for key in cazy_dict:
                if cazy_classes[index] in cazy_dict[key]:  # if in synonyms in cazy_dict
                    cazy_classes[index] = key  # standardise the class name

        # check all names are standardised, remove names that could not be standardised
        if cazy_classes[index] not in std_class_names:
            logger.warning(
                (
                    f"'{cazy_classes[index]}' could not be standardised.\n"
                    "Please use a synonym in the file_io/cazy_dictionary.json.\n"
                    f"'{cazy_classes[index]}' will NOT be scraped."
                )
            )
            del cazy_classes[index]

    return cazy_classes


def get_cmd_defined_fams_classes(cazy_dict, std_class_names, args):
    """Retrieve classes and families specified for scraping from the cmd-line args.

    :param cazy_dict: dict, accepted synonyms of CAZy class names
    :param std_class_names: list, standardised CAZy class names
    :param args: cmd args parser

    Return dictionary of CAZy classes and families to be scraped.
    """
    logger = logging.getLogger(__name__)

    # create dictionary which will store families and classes to be scraped
    config_dict = {
        'classes': [],
        'Glycoside Hydrolases (GHs)': [],
        'GlycosylTransferases (GTs)': [],
        'Polysaccharide Lyases (PLs)': [],
        'Carbohydrate Esterases (CEs)': [],
        'Auxiliary Activities (AAs)': [],
        'Carbohydrate-Binding Modules (CBMs)': [],
    }

    # add classes to config dict
    cazy_classes = args.classes
    if cazy_classes is not None:
        cazy_classes = cazy_classes.strip().split(",")
        # Standardise CAZy class names
        cazy_classes = parse_user_cazy_classes(cazy_classes, cazy_dict, std_class_names)
        config_dict["classes"] = cazy_classes

    families = args.families
    if families is not None:
        # add families to config dict
        families = families.strip().split(",")

        for fam in families:
            try:
                if fam.find("GH") != -1:
                    re.match(r"GH\d+", fam, re.IGNORECASE).group()
                    config_dict['Glycoside Hydrolases (GHs)'].append(fam)

                elif fam.find("GT") != -1:
                    re.match(r"GT\d+", fam, re.IGNORECASE).group()
                    config_dict['GlycosylTransferases (GTs)'].append(fam)

                elif fam.find("PL") != -1:
                    re.match(r"PL\d+", fam, re.IGNORECASE).group()
                    config_dict['Polysaccharide Lyases (PLs)'].append(fam)

                elif fam.find("CE") != -1:
                    re.match(r"CE\d+", fam, re.IGNORECASE).group()
                    config_dict['Carbohydrate Esterases (CEs)'].append(fam)

                elif fam.find("AA") != -1:
                    re.match(r"AA\d+", fam, re.IGNORECASE).group()
                    config_dict['Auxiliary Activities (AAs)'].append(fam)

                elif fam.find("CBM") != -1:
                    re.match(r"CBM\d+", fam, re.IGNORECASE).group()
                    config_dict['Carbohydrate-Binding Modules (CBMs)'].append(fam)

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


def get_excluded_classes(std_class_names, config_dict, cazy_dict):
    """Define the CAZy classes that will not be scraped.

    This includes classes for which not Families have been specified for scraping.

    :param std_class_names: list of standardised CAZy class names
    :param config_dict: configuration dict defining classes and families to be scraped
    :param cazy_dict: dict, accepted CAZy classes synonyms

    Return list of CAZy classes not to be scraped.
    """
    # retrieve list of CAZy classes from which all families are to be scraped
    cazy_classes = config_dict["classes"]

    # retrieve the names of classes for which specific families to be scraped have been named
    for key in config_dict:
        if (key != "classes") and (key not in cazy_classes) and (len(config_dict[key]) != 0):
            # add the class of families to be scraped to the list of CAZy classes to be scraped
            cazy_classes.append(key)

    # create list of CAZy classes not to be scraped
    excluded_classes = std_class_names
    excluded_classes = list(set(excluded_classes).difference(cazy_classes))

    if len(excluded_classes) != 0:
        # change names of classes into format for excluding classes during scrape
        index = 0
        for index in range(len(excluded_classes)):
            excluded_classes[index] = f"<strong>{excluded_classes[index]}</strong>"
    else:
        excluded_classes = None

    return excluded_classes


def convert_lists_to_none(config_dict):
    """Convert empty lists to None type objects in configuration dictionary.

    :param config_dict: dict, CAZy classes and families to be scraped

    Return dictionary with no empty lists.
    """
    for key in config_dict:
        if config_dict[key] is not None:
            if len(config_dict[key]) == 0:
                config_dict[key] = None

    return config_dict


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


def get_configuration(file_io_path, args):
    """Get configuration for the Expand module.

    :param file_io_path: Path to file_io module
    :param args: cmd-line argument parser

    Return configuration dictionary, set of taxonomy filters and set of Taxonomy Kingdoms.
    """
    # retrieve inital parsing of configuration data
    (
        excluded_classes, config_dict, cazy_dict, taxonomy_filters_dict, kingdoms, ec_filter,
    ) = parse_configuration(
        file_io_path, args,
    )
    # excluded_classes and cazy_dict are used in the crawler module but are not needed for the
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
