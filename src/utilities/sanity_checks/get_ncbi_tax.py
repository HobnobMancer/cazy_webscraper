#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
# (c) James Hutton Institute 2024
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
import sys

from src.utilities import termcolour


logger = logging.getLogger(__name__)


def sanity_check_inputs(args: argparse.Namespace) -> None:
    if (args.genbank_accessions or args.uniprot_accessions) and args.file_only:
        error_message = """::ARGUMENT ERROR:: Cannot provide a list of protein accessions (GenBank and/or UniProt) 
        and disable retrieval dating from NCBI (using the --file_only flag).
        """
        logger.error(termcolour(error_message, "red"))
        sys.exit(22)
