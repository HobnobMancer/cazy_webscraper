#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
# (c) James Hutton Institute 2024
# Author:
# Emma E. M. Hobbs

# Contact
# ehobbs@ebi.ac.uk

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
"""Build the main CLI"""


import sys

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, Namespace
from typing import List, Optional

from pathlib import Path
from typing import List, Optional

from src import __version__, __citation__


def build_parser(argv: Optional[List] = None) -> Namespace:
    """The main parser for cazy_webscraper"""
    # Create parser object
    parser_main = ArgumentParser(
        prog="cazy_webscraper",
        description="Web scraper to retrieve protein data catalogued by CAZy, UniProt, NCBI, GTDB and PDB.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser_main.add_subparsers(
        title="subcommands", description="Valid subcommands",
    )

    # add args to the main parser
    parser_main.add_argument(
        "--version",
        action="store_true",
        default=False,
        help="Print version number"
    )
    parser_main.add_argument(
        "--citation",
        action="store_true",
        default=False,
        help="Print citation information"
    )
    parser_main.add_argument(
        "-l",
        "--log",
        dest="log",
        action="store",
        default=None,
        type=Path,
        help="logfile location",
    )
    parser_main.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        help="report verbose progress to log",
    )

    # add subcommand parser

    # Parse arguments
    # The list comprehension is to allow PosixPaths to be defined and passed in testing
    if argv is None:
        if len(sys.argv) == 1:
            argv = ["-h"]
        else:
            argv = sys.argv[1:]
    return parser_main.parse_args([str(_) for _ in argv])
