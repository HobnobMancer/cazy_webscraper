#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2024
# (c) University of Strathclyde 2024
# (c) James Hutton Institute 2024
#
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
"""Parse files of accessions"""


import logging
import sys

from pathlib import Path

from src.sql.interface.filter_data.protein import filter_to_db_acc


logger = logging.getLogger(__name__)


def get_acc_from_file(
    acc_file: Path,
    db: Path,
    table="Protein"
) -> list[str]:
    """Get accession from user file and add to seq_acc_to_retrieve"""
    logger.warning("Getting %s accessions from file: %s", table, acc_file)
    seq_acc_to_retrieve = set()
    try:
        with open(acc_file, "r") as fh:
            for line in fh:
                seq_acc_to_retrieve.add(line.strip())
    except FileNotFoundError:
        logging.error(
            "Could not find list of GenBank accessions at: %s\n"
            "Check the path is correct\n"
            "Terminating program", acc_file
        )
        sys.exit(1)

    seq_acc_to_retrieve = filter_to_db_acc(db, seq_acc_to_retrieve, table=table)

    # return a list, not a set: downstream batching (get_chunks_list) slices its input,
    # and sets are not subscriptable. sorted() also makes batch composition reproducible.
    return sorted(seq_acc_to_retrieve)
