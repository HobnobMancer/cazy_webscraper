#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
"""Single entry point for application"""


import logging
import sys

import pandas as pd

from datetime import datetime
from typing import List, Optional

from saintBioutils.utilities.logger import config_logger

from src import CITATION_INFO, VERSION_INFO, closing_message
from src.utilities.parsers.parse_cmd import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    if argv is None:
        args = build_parser()
    else:
        args = build_parser(argv)

    if args.citation:
        sys.stderr.write(f"{VERSION_INFO}\n")
        sys.stderr.write(f"{CITATION_INFO}")
        return 0
    if args.version:
        sys.stderr.write(f"{VERSION_INFO}\n")
        return 0

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    # Boilerplate for log
    logger.info("Processed arguments: %s", args)
    args.cmdline = " ".join(sys.argv)
    logger.info("command-line: %s", args.cmdline)

    # run subcommand
    returnval = args.func(args, time_stamp, start_time)
    print(returnval)
    closing_message("cazy_webscraper", start_time, args)
    return returnval


if __name__ == "__main__":
    main()
