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

import pandas as pd

from datetime import datetime
from pathlib import Path

from saintBioutils.utilities.file_io import make_output_directory

from cazy_webscraper.sql import sql_orm

__version__ = "2.1.3"

VERSION_INFO = f"cazy_webscraper version: {__version__}"

CITATION_INFO = (
    "If you use cazy_webscraper in your work, please cite the following publication:\n"
    "\tHobbs, E. E. M., Pritchard, L., Chapman, S., Gloster, T. M.,\n"
    "\t(2021) cazy_webscraper Microbiology Society Annual Conference 2021 poster.\n"
    "\tFigShare. Poster.\n"
    "\thttps://doi.org/10.6084/m9.figshare.14370860.v7"
)


def closing_message(job, start_time, args):
    """Write closing messsage to terminal"""
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    if args.verbose:
        logger.info(
            f"====================={job}=====================\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}"
            f"Version: {VERSION_INFO}\n"
            f"Citation: {CITATION_INFO}"
        )
    else:
        print(
            f"====================={job}=====================\n"
            f"Scrape initated at {start_time}\n"
            f"Scrape finished at {end_time}\n"
            f"Total run time: {total_time}\n"
            f"Version: {VERSION_INFO}\n"
            f"Citation: {CITATION_INFO}"
        )

    return


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
