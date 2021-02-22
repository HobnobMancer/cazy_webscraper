#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

# The MIT License
"""Module for handling input and output files and directories."""


import json
import logging
import re
import shutil
import sys
import yaml

from pathlib import Path


def make_output_directory(output, force, nodelete):
    """Create output directory for genomic files.

    :param output: path, path of dir to be created
    :param force: bool, enable/disable creating dir if already exists
    :param nodelete: bool, enable/disable deleting content in existing dir

    Raises FileExistsError if an attempt is made to create a directory that already
    exist and force (force overwrite) is False.

    Return Nothing
    """
    logger = logging.getLogger(__name__)

    if force is True:
        logger.warning(
            "Output directory %s exists, nodelete is %s", output, nodelete,
        )

        if nodelete and output.exists():
            logger.warning("Not deleting directory %s", output)

        elif output.exists():
            logger.warning("Deleting directory %s", output)
            shutil.rmtree(output)

    logger.info("Creating directory %s", output)
    try:
        output.mkdir(exist_ok=force)
    except FileExistsError:
        logger.warning(
            (
                "Out directory already exists."
                "New directory not made, writing to existing directory."
            )
        )

    return


def parse_configuration(file_io_path, args):
    """Parse configuration data, and retrieve user specified CAZy classes and families.

    If no user defined configuration is given the default behaviour to scrape the entirity of CAZy
    is enabled.

    :param file_io_path: str, path to directory where file_io is installed
    :param args: parser arguments

    Return list of classes not to scrape, dict of families to scrape, dict of class synonoms, and
    a set of taxonomy filters (genera, species and strains) to restrict the scrape to.
    """
    # retrieve taxonomy filters
    taxonomy_filter = get_genera_species_strains(args)

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

            return excluded_classes, config_dict, cazy_dict, taxonomy_filter

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

            return excluded_classes, config_dict, cazy_dict, taxonomy_filter

    else:  # user did not pass config file
        if (args.classes is None) and (args.families is None):
            # No specific families or classes specified for scraping
            return None, None, cazy_dict, taxonomy_filter

        else:  # configuration specified only via the cmd_line
            config_dict = get_cmd_defined_fams_classes(cazy_dict, std_class_names, args)

            # get list of CAZy classes that will not be scraped
            excluded_classes = get_excluded_classes(std_class_names, config_dict, cazy_dict)

            # convert empty lists to None type objects
            config_dict = convert_lists_to_none(config_dict)

            return excluded_classes, config_dict, cazy_dict, taxonomy_filter


def get_genera_species_strains(args):
    """Retrieve Genera, Species and Strains to retrict the scrape to.

    :param args: cmd-line arguments parser

    Return set contain user specificed genera, species and strains.
    """
    taxonomy_filter = []

    # open config dict
    if args.config is not None:
        try:
            with open(args.config, "r") as fh:
                raw_config_dict = json.load(fh)
        except FileNotFoundError:
            logger.error(
                "Did not find the configuration file. Check the path is correct.\n"
                "Terminating programme"
            )
            sys.exit(1)
        
        for key in ["genera", "species", "strain"]:
            try:
                if raw_config_dict[key] is not None:
                    taxonomy_filter += raw_config_dict[key]
            except KeyError:
                pass
    
    if args.genera is not None:
        taxonomy_filter += args.genera
    
    if args.species is not None:
        taxonomy_filter += args.species
    
    if args.strains is not None:
        taxonomy_filter += args.strains
    
    if len(taxonomy_filter) == 0:
        taxonomy_filter = None
    else:
        taxonomy_filter = set(taxonomy_filter)
    
    return taxonomy_filter


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
            "Could not open the CAZy synonym dictionary.\n"
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
        logger.warning(
            "Could not find configuration file when option was enabled.\n"
            "Make sure path to the configuration file is correct\n"
            "Scrapping will not be performed becuase configuration is wrong.\n"
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
        if (key != "classes") and (key != "genera") and (key != "species") and (key != "strains"):
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
    except (KeyError, TypeError) as e:
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


def write_out_failed_scrapes(failed_urls, time_stamp, args):
    """Write out the URLs for which a connection to CAZy failed.

    :param failed_urls: list, contains the URL and reason for the failed scrape
    :param args: cmd args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    if args.output is not sys.stdout:
        output_path = args.output / f"failed_cazy_scrapes_{time_stamp}.txt"

        with open(output_path, "a") as fh:
            for url in failed_urls:
                fh.write(f"{url}\n")

    else:
        logger.error("The following items were not scraped:")
        for url in failed_urls:
            logger.error(url)

    return


def write_out_failed_proteins(sql_failures, time_stamp, args):
    """Write out the names of proteins which raised errors when being added to the local db.

    :param sql_failures: list, the names of proteins that were unsuccessfully added to the db
    :param args: cmd args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    if args.output is not sys.stdout:
        output_path = args.output / f"failed_db_protein_additions_{time_stamp}.txt"

        with open(output_path, "a") as fh:
            for fail in sql_failures:
                fh.write(f"{fail}\n")

    else:
        logger.error("The following proteins were not entered into database:")
        for fail in sql_failures:
            logger.error(fail)

    return
