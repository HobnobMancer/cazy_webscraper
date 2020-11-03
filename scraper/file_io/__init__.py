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
import shutil
import yaml

import pandas as pd


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


def parse_configuration(args, logger):
    """Parse configuration data, and retrieve user specified CAZy classes and families.

    Return nothing if a path to a configuration file was not given. Cases default behaviour of the
    scrapper to be performed. Return list CAZy classes not to be scraped, list CAZy families to be
    scraped and dictionary of CAZy synonoms.

    :param args: parser arguments
    :param logger: logger object

    Return two lists and a dictionary.
    """
    # open dictionary of accepted CAZy class synonyms)
    with open("cazy_dictionary.json", "r") as fh:
        cazy_dict = json.load(fh)
        # create list of standardised CAZy classes
        class_names = list(cazy_dict.keys())

    if args.config is None:
        return None, None, cazy_dict

    # open configuration file
    with open(args.config) as fh:
        config_dict = yaml.full_load(fh)

    # Retrieve CAZy classes listed in config file
    try:
        cazy_classes = config_dict["classes"]
        # standardise CAZy class names
        cazy_classes = parse_user_cazy_classes(cazy_classes, cazy_dict, class_names, logger)
    except (KeyError, TypeError) as e:
        logger.info("No CAZy classes specified in configuration file")
        cazy_classes = None

    if len(cazy_classes) == 0:
        logger.info("No CAZy classes specified in configuration file")
        cazy_classes = None

    # retrieve classes of families/subfamilies specified in configuration file
    for key in config_dict:
        if (key != "classes") and (len(config_dict[key]) > 0):
            # add the class of families to be scraped to the list of CAZy classes to be scraped
            if key not in cazy_classes:
                try:
                    cazy_classes.append(config_dict[key])
                except TypeError:
                    # raised if cazy_classes is None
                    cazy_classes = [config_dict[key]]

    # create list of CAZy classes not to be scraped
    try:
        excluded_classes = class_names
        excluded_classes = list(set(excluded_classes).difference(cazy_classes))
        # change names of classes into format for excluding classes during scrape
        index = 0
        for index in range(len(excluded_classes)):
            excluded_classes[index] = f"<strong>{excluded_classes[index]}</strong>"
    except TypeError:
        # raised if CAZy classes is None
        excluded_classes = None

    return excluded_classes, config_dict, cazy_dict


def parse_user_cazy_classes(cazy_classes, cazy_dict, class_names, logger):
    """Standardise the CAZy class names listed in configuration file.

    :param cazy_classes: list, list of CAZy classes from configuration file
    :param cazy_dict: dict, keyed by class name, keyed by list of accepted synonoms
    :param class_names: list, list of all CAZy classes
    :param logger: logger object

    Return list of CAZy classes named in configuration file.
    """
    # identify user defined CAZy classes not written in standardised format
    index = 0
    for index in range(len(cazy_classes)):
        if cazy_classes[index] not in class_names:
            for key in cazy_dict:
                if cazy_classes[index] in cazy_dict[key]:
                    cazy_classes[index] = key
        # check all names are standardised, remove names that could not be standardised
        if cazy_classes[index] not in class_names:
            logger.warning(
                (
                    f"'{cazy_classes[index]}' could not be standardised.\n"
                    "Please use a synonym in the file_io/cazy_dictionary.json.\n"
                    f"'{cazy_classes[index]}' will NOT be scraped."
                )
            )
            del cazy_classes[index]

    return class_names, cazy_classes


def write_out_df(dataframe, df_name, outdir, logger, force):
    """Write out dataframe to output directory.

    :param dataframe: pandas dataframe
    :param df_name: str, name of dataframe
    :param outdir: Path, path to output directory
    :param logger: logger object
    :param force: bool, enable/disable over writing of existing file

    Return nothing.
    """
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
    return
