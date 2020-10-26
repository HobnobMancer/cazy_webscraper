#!/usr/bin/env python
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


def write_out_df(dataframe, df_name, out_dir, logger, force):
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
