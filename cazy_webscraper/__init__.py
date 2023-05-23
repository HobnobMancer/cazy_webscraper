#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Web scraper to scrape the CAZy website."""


import logging
import os
import sys

import numpy as np
import pandas as pd

from datetime import datetime
from pathlib import Path

import Bio
import bioservices
import bs4
import html5lib
import lxml
import mechanicalsoup
import requests
import saintBioutils
import sqlalchemy
import tqdm

from saintBioutils.utilities.file_io import make_output_directory

from cazy_webscraper.sql import sql_orm


__version__ = "2.3.0"


VERSION_INFO = f"cazy_webscraper version: {__version__}"


CITATION_INFO = (
    "If you use cazy_webscraper in your work, please cite the following publication:\n"
    "\tHobbs, E. E. M., Gloster, T. M., and Pritchard, L.\n"
    "\t(2022) 'cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets',\n"
    "\tbioRxiv\n"
    "\thttps://doi.org/10.1101/2022.12.02.518825"
)

WEBSITE = "https://hobnobmancer.github.io/cazy_webscraper/"

DOCUMENTATION = "https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest"

GITHUB_ISSUES = "https://github.com/HobnobMancer/cazy_webscraper/issues"

AUTHOR_EMAIL = "eemh1@st-andrews.ac.uk"


def closing_message(job, start_time, args, early_term=False):
    """Write closing messsage to terminal

    :param job: str, name of module run
    :param start_time: str, time run was started
    :param args: CLI arguments parser
    :param early_term: bool, True if run terminated early due to an error
    """
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    if early_term:
        termination_status = (
            "**Run terminated early due to do an error**\n"
            f"Run finished at {end_time}"
        )
    else:
        termination_status = f"Run finished at {end_time}"

    message = f"""
    ====================={job}=====================
    Run initiated at {start_time}
    {termination_status}
    Total run time: {total_time}

    Version: {VERSION_INFO}

    For help with trouble shooting and operating cazy_webscraper please see the documentation:
    README: {WEBSITE}
    Documentation (including tutorials): {DOCUMENTATION}
    GitHub Issues: {GITHUB_ISSUES}

    When publishing work that uses cazy_webscraper please cite:
    Citation: {CITATION_INFO}
    """

    if args.verbose:
        logger.info(message)
    else:
        print(message)

    if early_term:
        sys.exit(1)


def display_citation_info():
    """Display citation inforamtion.

    Return nothing
    """

    message = f"""
    =====================cazy_webscraper Citation Information=====================
    cazy_webscraper version: {VERSION_INFO}

    When publishing work that uses cazy_webscraper please cite:
    Citation: {CITATION_INFO}

    cazy_webscraper depends on a number of tools. To recognise the contributions that the 
    authors and developers have made, please also cite the following:

    When making an SQLite database:
    Hipp, R. D. (2020) SQLite, available: https://www.sqlite.org/index.html.

    Retrieving taxonomic, genomic or sequence data from NCBI:
    Cock, P.J.A., Antao, T., Chang, J.T., Chapman, B.A., Cox, C.J., Dalke, A., et al (2009) 
    Biopython: freely available Python tools for computational molecular biology and 
    bioinformatics, Bioinformatics, 25(11), 1422-1423.
    Wheeler,D.L., Benson,D.A., Bryant,S., Canese,K., Church,D.M., Edgar,R., Federhen,S.,
    Helmberg,W., Kenton,D., Khovayko,O. et al (2005) Database resources of the National Centre
    for Biotechnology Information: Update, Nucleic Acid Research, 33, D39-D45

    Retrieving data from UniProt:
    Cokelaer, T., Pultz, D., Harder, L. M., Serra-Musach, J., Saez-Rodriguez, J. (2013)
    BioServices: a common Python package to access biological Web Services programmatically,
    Bioinformatics, 19(24), 3241-3242.

    Downloading protein structure files from RSCB PDB:
    Berman, H.M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T.N., Weissig, H., et al (2022)
    The Protein Data Bank, Nucleic Acids Research, 28(1), 235-242.
    Hamelryck, T., Manderick, B. (2003), PDB parser and structure class implemented in Python.
    Bioinformatics, 19 (17), 2308–2310

    Retrieving and using taxonomic data from GTDB:
    Parks, D.H., Chuvochina, M., Rinke, C., Mussig, A.J., Chaumeil, P., Hugenholtz, P. (2022)
    GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically
    consistent, rank normalized and complete genome-based taxonomy, Nucleic Acids Research,
    50(D1), D785-D794.

    """

    print(message)


def display_version_info():
    """Display package version number information"""

    try:
        saintbio_version = saintBioutils.__version__
    except AttributeError:
        saintbio_version = 'unknown'

    message = f"""
    =====================cazy_webscraper Version Information=====================
    cazy_webscraper version: {VERSION_INFO}

    Third party tools used by cazy_webscraper:
    beautifulsoup4: {bs4.__version__}
    biopython: {Bio.__version__}
    bioservices: {bioservices.version}
    html5lib: {html5lib.__version__}
    lxml: {lxml.__version__}
    mechanicalsoup: {mechanicalsoup.__version__}
    numpy: {np.__version__}
    pandas: {pd.__version__}
    requests: {requests.__version__}
    saintBioutils: {saintbio_version}
    sqlalchemy: {sqlalchemy.__version__}
    tqdm: {tqdm.__version__}
    """

    print(message)


def connect_existing_db(args, time_stamp, start_time):
    """Coordinate connecting to an existing local CAZyme database, define logger name and cache dir

    :param args: cmd-line args parser
    :param time_stamp: str, time cazy_webscraper was invoked
    :param start_time: pd date-time obj, time cazy_webscraper was invoked

    Return connection to local CAZyme database, logger file name, and path to cache dir
    """
    logger = logging.getLogger(__name__)

    logger.info("Adding data to an existing local CAZyme database")

    if os.path.isfile(args.database) is False:
        logger.error(
            "Could not find local CAZy database.\n"
            "Check path is correct.\n"
            "Terminating programme."
        )
        closing_message("cazy_webscraper", start_time, args)
        sys.exit(1)

    try:
        connection = sql_orm.get_db_connection(args.database, args.sql_echo, new=False)
        logger.info("Opened connection to local CAZyme database")
    except Exception:
        logger.error(
            "Failed to open connection to an exiting local CAZyme database\n."
            "Terminating program\n",
            exc_info=True,
        )
        closing_message("cazy_webscraper", start_time, args)
        sys.exit(1)

    # used for naming additional log files
    logger_name = str(args.database).split('.')[0]

    # define path to cache family txt files
    cache_dir = Path(f"{str(args.database.parent)}/.cazy_webscraper_{time_stamp}")

    return connection, logger_name, cache_dir


def connect_to_new_db(args, time_stamp, start_time):
    """Build and connect to a new local CAZyme database.

    :param args: cmd-line args parser
    :param time_stamp: str, time cazy_webscraper was invoked
    :param start_time: pd date-time obj, time cazy_webscraper was invoked

    Return connection to the database, name of the logger, and path to the cache dir
    """
    logger = logging.getLogger(__name__)

    if args.db_output is not None:  # user defined target output for the NEW database

        if os.path.isfile(args.db_output):  # target file exists
            if args.force:
                logger.warning(
                    "Overwriting existing local CAZyme database at:\n"
                    f"{args.db_output}"
                )

            else:
                logger.warning(
                    "Target path for new database already exists.\n"
                    "Either enable forced overwriting (-f) or add data this data (-D).\n"
                    "Terminating program."
                )
                closing_message("cazy_webscraper", start_time, args)
                sys.exit(1)

        else:  # may need to build dirs
            logger.info(
                "Building new local CAZyme database\n"
                f"Output directory: {(args.db_output).parent}\n"
                f"Force overwriting exiting output file: {args.force}"
            )

        if str((args.db_output).parent) != '.':  # dirs defined in output put
            output_dir = (args.db_output).parent
            make_output_directory(output_dir, args.force, args.nodelete)
            cache_dir = Path(f"{str(output_dir)}/.cazy_webscraper_{time_stamp}")

        else:  # writing to cwd
            cache_dir = Path(f".cazy_webscraper_{time_stamp}")

        logger_name = str(args.db_output).split('.')[0]
        db_path = args.db_output

    else:
        logger.info("Using default database name and writing to cwd")
        db_path = Path(f"cazy_webscraper_{time_stamp}.db")
        cache_dir = Path(f".cazy_webscraper_{time_stamp}")
        logger_name = f'cazy_webscraper_{time_stamp}'

    try:
        connection = sql_orm.get_db_connection(db_path, args.sql_echo, new=True)
        logger.warning(f"Built new local CAZyme database at\n{db_path}")
    except Exception:
        logger.error(
            "Failed to build new SQL database\n."
            "Terminating program",
            exc_info=True,
        )
        closing_message("cazy_webscraper", start_time, args)
        sys.exit(1)

    return connection, logger_name, cache_dir
