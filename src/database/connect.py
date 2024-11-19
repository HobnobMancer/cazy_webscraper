#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
# (c) James Hutton Institute 2024
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

import logging
import os
import sys

from pathlib import Path

from saintBioutils.utilities.file_io import make_output_directory

from src.sql import sql_orm
from src import closing_message


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
    logger_name = str(args.database).split('.', maxsplit=1)[0]

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

    if args.db_output:  # user defined target output for the NEW database
        if str((args.db_output).parent) != '.':  # dirs defined in output
            output_dir = (args.db_output).parent
            if not output_dir.exists():
                output_dir.mkdir(exist_ok=True, parents=True)
            cache_dir = Path(f"{str(output_dir)}/.cazy_webscraper_{time_stamp}")

        else:  # writing to cwd
            cache_dir = Path(f".cazy_webscraper_{time_stamp}")

        logger_name = str(args.db_output).split('.', maxsplit=1)[0]
        db_path = args.db_output

    else:
        logger.info("Using default database name and writing to cwd")
        db_path = Path(f"cazy_webscraper_{time_stamp}.db")
        cache_dir = Path(f".cazy_webscraper_{time_stamp}")
        logger_name = f'cazy_webscraper_{time_stamp}'

    logger.info("Building new local CAZyme database: %s\n", args.db_output)

    try:
        connection = sql_orm.get_db_connection(db_path, args.sql_echo, new=True)
        logger.warning("Built new local CAZyme database at%s", db_path)
    except Exception:
        logger.error(
            "Failed to build new SQL database\n."
            "Terminating program",
            exc_info=True,
        )
        closing_message("cazy_webscraper", start_time, args)
        sys.exit(1)

    return connection, logger_name, cache_dir
