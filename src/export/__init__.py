#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Write data held in a local CAZyme database out to files.

This is the mirror image of ``src/databases`` - that package brings data in from external
sources, this one takes data back out of the local database. Nothing here talks to an
external service.
"""


import logging

from argparse import Namespace
from pathlib import Path


logger = logging.getLogger(__name__)


def compile_output_path(args: Namespace) -> Path:
    """Build the stem shared by every output file of a run.

    e.g. <output_dir>/<prefix>_<database name>
    """
    stem = args.database.name
    if stem.endswith(".db"):
        stem = stem[:-3]

    if args.prefix:
        stem = f"{args.prefix}_{stem}"

    output_dir = args.output_dir if args.output_dir else Path.cwd()

    return output_dir / stem


def check_existing_outputs(paths: list[Path], args: Namespace) -> bool:
    """Warn about output files that already exist.

    Returns True if it is safe to continue, False if the run should stop.
    """
    existing = [path for path in paths if path.exists()]

    if not existing:
        return True

    listed = "\n  ".join(str(path) for path in existing)

    if args.overwrite:
        logger.warning("Overwriting existing output file(s):\n  %s", listed)
        return True

    logger.error(
        "The following output file(s) already exist:\n  %s\n"
        "Use --overwrite to overwrite them, or --prefix to write under a different name.\n"
        "Terminating program",
        listed,
    )
    return False
