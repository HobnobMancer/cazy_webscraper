#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
# Author:
# Emma E. M. Hobbs
#
# Contact
# ehobbs@ebi.ac.uk
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
"""Sanity checks for user inputs"""


import argparse
import logging
import os
import sys

from pathlib import Path

from src.utilities import termcolour


logger = logging.getLogger(__name__)


def sanity_check_main_input(time_stamp: str, args: argparse.ArgumentParser) -> None:
    db = None

    if not args.email and not args.skip_ncbi_tax:
        error_message = """::ARGUMENT ERROR:: No email address provided.
        An email address is required by NCBI for retrieving the latest taxonomic classifications 
        for proteins listed with multiple source organisms in the CAZy database.
        Please provide an email address. Terminating program."""
        logger.error(termcolour(error_message, "red"))
        sys.exit(22)

    if args.database and args.db_output:
        error_message = """::ARGUMENT ERROR:: New and existing database paths specified
            A target path for a NEW database (--db_output, -d) and a path to an EXISTING database (--database, -D) were provided.
            Please provide one OR the other.
            Terminating program."""
        logger.error(termcolour(error_message, "red"))
        sys.exit(22)

    if args.db_output and args.db_output.exists() and not args.force:
        error_message = """::INPUT ERROR:: Not allowed to overwrite existing database.
            A database already exists at %s
            and --force/-f was not used, therefore, cazy_webscraper cannot overwrite this database
            Terminating program.""" % args.db_output
        logger.error(termcolour(error_message, "red"))
        sys.exit(5)

    if args.db_output and args.db_output.exists() and args.force:
        error_message = "Local db %s already exists. Force is True therefore, ovewriting existing database." % args.db_output
        logger.warning(termcolour(error_message, "yellow"))
        os.remove(args.db_output)

    if not any([args.db_output, args.database]):
        db = Path(f"cazy_webscraper_{time_stamp}.db")
        error_message = """No database name provided.
        Using default name (%s) and writing to PWD""" % db
        logger.warning(termcolour(error_message, "yellow"))
        
    else:
        if args.db_output:
            db = args.db_output
        else:
            db = args.database

    return db
