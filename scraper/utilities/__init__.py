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
"""Module containing functions for building parsers, loggers and parsing configuration."""


import logging


def config_logger(args) -> logging.Logger:
    """Configure package wide logger.

    Configure a logger at the package level, from which the module will inherit.
    If CMD-line args are provided, these are used to define output streams, and
    logging level.

    :param args: cmd-line args parser

    Return nothing
    """
    logger = logging.getLogger(__package__)

    # Set format of loglines
    log_formatter = logging.Formatter("[%(levelname)s] [%(name)s]: %(message)s")

    # define logging level
    if args.verbose is True:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    # Setup console handler to log to terminal
    console_log_handler = logging.StreamHandler()
    console_log_handler.setFormatter(log_formatter)
    logger.addHandler(console_log_handler)

    # Setup file handler to log to a file
    if args.log is not None:
        file_log_handler = logging.FileHandler(args.log)
        file_log_handler.setFormatter(log_formatter)
        logger.addHandler(file_log_handler)

    return
