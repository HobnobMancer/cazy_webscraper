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
import re
import shutil
import sys
import yaml

from pathlib import Path


def make_output_directory(output, logger, force, nodelete):
    """Create output directory for genomic files.

    :param output: path, path of dir to be created
    :param logger: logger object
    :param force: bool, enable/disable creating dir if already exists
    :param nodelete: bool, enable/disable deleting content in existing dir

    Raises FileExistsError if an attempt is made to create a directory that already
    exist and force (force overwrite) is False.

    Return Nothing
    """
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


def parse_configuration(file_io_path, args, logger):
    """Parse configuration data, and retrieve user specified CAZy classes and families.

    Return only the CAZy class synonoms dictionary if a path to a configuration file was not given.
    This results in the default behaviour of the webscraper to scrape the entirty of CAZy to be
    invoked. If no items are listed under a heading/tag in the config file, the retrieved item will
    be a None type object.

    :param file_io_path: str, path to directory where file_io is installed
    :param args: parser arguments
    :param logger: logger object

    Return list of classes not to scrape, dict of families to scrape, and dict of class synonoms.
    """
    # open dictionary of accepted CAZy class synonyms)
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

    if args.config is None:
        if (args.classes is None) and (args.families is None):
            # No specific families or classes specified for scraping
            return None, None, cazy_dict

        else:
            config_dict = get_cmd_defined_fams_clsses(args)
            return 

    try:
        # open configuration file
        with open(args.config) as fh:
            config_dict = yaml.full_load(fh)
    except FileNotFoundError:
        logger.warning(
            "Could not find configuration file.\n"
            "Make sure path to the configuration file is correct\n"
            "Scrapping will continue using default scraper configuration\n"
        )
        return None, None, cazy_dict

    # Retrieve CAZy classes listed in config file
    try:
        cazy_classes = config_dict["classes"]
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
        cazy_classes = parse_user_cazy_classes(cazy_classes, cazy_dict, std_class_names, logger)

    # retrieve classes of families/subfamilies specified in configuration file
    for key in config_dict:
        if (key != "classes") and (key not in cazy_classes) and (config_dict[key] is not None):
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

    return excluded_classes, config_dict, cazy_dict


def get_cmd_defined_fams_clsses(args, logger):
    """Retrieve classes and families specified for scraping from the cmd-line args.

    :param args: cmd args parser.

    Return two lists, one list of families, one line of classes.
    """
    # create dictionary which will store families and classes to be scraped
    config_dict = {
        'classes': None,
        'Glycoside Hydrolases (GHs)': [],
        'GlycosylTransferases (GTs)': [],
        'Polysaccharide Lyases (PLs)': [],
        'Carbohydrate Esterases (CEs)': [],
        'Auxiliary Activities (AAs)': [],
        'Carbohydrate-Binding Modules (CBMs)': [],
    }

    # add classes to config dict
    classes = args.cazy_classes
    if classes is not None:
        classes = classes.strip().split(",")
        config_dict["classes"] = classes

    families = args.family
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

    # Convert empty list to None type objects
    for key in config_dict:
        if config_dict[key] is not None:
            if len(config_dict[key]) == 0:
                config_dict[key] = None

    return config_dict


def parse_user_cazy_classes(cazy_classes, cazy_dict, std_class_names, logger):
    """Standardise the CAZy class names listed in configuration file.

    :param cazy_classes: list, list of CAZy classes from configuration file
    :param cazy_dict: dict, keyed by class name, keyed by list of accepted synonoms
    :param std_class_names: list, list of all CAZy classes
    :param logger: logger object

    Return list of CAZy classes listed by user in the configuration file.
    """
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


def write_out_df(dataframe, df_name, outdir, logger, force):
    """Write out dataframe to output directory.

    :param dataframe: pandas dataframe
    :param df_name: str, name of dataframe
    :param outdir: Path, path to output directory
    :param logger: logger object
    :param force: bool, enable/disable over writing of existing file

    Return nothing.
    """
    if outdir is not sys.stdout:
        # build output path
        output_path = outdir / f"{df_name}.csv"

        logger.info("Checking if output directory for dataframe already exists")
        if output_path.exists():
            if force is False:
                logger.warning(
                    (
                        "Specified directory for dataframe already exists.\n"
                        "Exiting writing out dataframe."
                    )
                )
                return ()
            else:
                logger.warning(
                    (
                        "Specified directory for dataframe already exists.\n"
                        "Forced overwritting enabled."
                    )
                )

        logger.info("Writing out species dataframe to directory")
        dataframe.to_csv(output_path)

    else:
        logger.info(f"Output dataframe: {df_name}")
        dataframe.to_csv(outdir)

    return


def write_out_failed_scrapes(failed_urls, time_stamp, args, logger):
    """Write out the URLs for which a connection to CAZy failed.

    :param failed_urls: list, contains the URL and reason for the failed scrape
    :param args: cmd args parser
    :param logger: logger object

    Return nothing.
    """
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
