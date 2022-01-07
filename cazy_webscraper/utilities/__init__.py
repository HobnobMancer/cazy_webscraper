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
"""Module containing functions for building parsers, loggers and parsing configuration."""


import logging
import os

from pathlib import Path
from typing import Optional


def config_logger(args) -> logging.Logger:
    """Configure package wide logger.

    Configure a logger at the package level, from which the module will inherit.
    If CMD-line args are provided, these are used to define output streams, and
    logging level.

    :param args: cmd-line args parser

    Return nothing
    """
    logger = logging.getLogger(__package__)
    logger.setLevel(logging.INFO)

    # define logging level
    if args.verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    # Set format of loglines
    log_formatter = logging.Formatter("[%(levelname)s] [%(name)s]: %(message)s")

    # Setup console handler to log to terminal
    console_log_handler = logging.StreamHandler()
    console_log_handler.setFormatter(log_formatter)
    if args.verbose:
        console_log_handler.setLevel(logging.INFO)
    else:
        console_log_handler.setLevel(logging.WARNING)
    logger.addHandler(console_log_handler)

    # Setup file handler to log to a file
    if args.log is not None:
        file_log_handler = logging.FileHandler(args.log)
        file_log_handler.setFormatter(log_formatter)
        if args.verbose:
            file_log_handler.setLevel(logging.INFO)
        else:
            file_log_handler.setLevel(logging.WARNING)
        logger.addHandler(file_log_handler)

    return


def build_logger(output):
    """Build loggers with pre-defined parameters for writing out errors and failed scrapes.

    :param output: Path to output dir or None

    Return logger object.
    """
    file_name = output.name.split('.')[0]
    logger = logging.getLogger(file_name)

    # Set format of loglines
    log_formatter = logging.Formatter(file_name + ": {} - {}".format("%(asctime)s", "%(message)s"))

    # Setup file handler to log to a file
    file_log_handler = logging.FileHandler(output)
    file_log_handler.setLevel(logging.WARNING)
    file_log_handler.setFormatter(log_formatter)
    logger.addHandler(file_log_handler)

    return logger


def termcolour(
    logstr: str, color: Optional[str] = None, bold: Optional[bool] = False
) -> str:
    """Return the passed logstr, wrapped in terminal colouring."""
    # For terminal colouring
    termcolours = {
        "BLACK": 0,
        "RED": 1,
        "GREEN": 2,
        "YELLOW": 3,
        "BLUE": 4,
        "MAGENTA": 5,
        "CYAN": 6,
        "WHITE": 7,
    }
    reset = "\033[0m"

    # Colour the string
    if isinstance(color, str) and color.upper() in termcolours:
        logstr = f"\033[1;{30 + termcolours[color.upper()]}m{logstr}{reset}"

    # Make the string bold
    if bold is True:
        logstr = f"\033[1m{logstr}"
        if not logstr.endswith(reset):
            logstr += reset

    return logstr
