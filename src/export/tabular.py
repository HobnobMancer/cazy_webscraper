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
"""Write data from a local CAZyme database out as CSV and JSON."""


import csv
import json
import logging

from pathlib import Path

from tqdm import tqdm


logger = logging.getLogger(__name__)


# the order columns appear in, when the corresponding --include option is used
COLUMN_ORDER = (
    "class",
    "family",
    "subfamily",
    "kingdom",
    "genus",
    "organism",
    "ec",
    "pdb",
    "uniprot_acc",
    "uniprot_name",
    "genbank_seq",
    "uniprot_seq",
)

# fields that hold a set of values rather than a single one
MULTI_VALUE_FIELDS = {"class", "family", "subfamily", "ec", "pdb"}


def build_column_names(include: set[str]) -> list[str]:
    """Return the CSV header for the requested --include options."""
    return ["protein_accession"] + [name for name in COLUMN_ORDER if name in include]


def _flatten(value):
    """Render one field for CSV: sets become comma separated, None becomes empty."""
    if value is None:
        return ""
    if isinstance(value, (set, frozenset, list, tuple)):
        return ",".join(sorted(str(v) for v in value))
    return str(value)


def write_csv(record_batches, output_path: Path, include: set[str]) -> int:
    """Write protein records to CSV, a batch at a time.

    Returns the number of rows written.
    """
    columns = build_column_names(include)
    rows_written = 0

    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)

        for batch in tqdm(record_batches, desc="Writing CSV"):
            for accession in sorted(batch):
                record = batch[accession]
                writer.writerow(
                    [accession] + [_flatten(record.get(name)) for name in columns[1:]]
                )
                rows_written += 1

    return rows_written


def write_json(record_batches, output_path: Path, include: set[str]) -> int:
    """Write protein records to JSON, a batch at a time.

    The document is assembled on the fly rather than built in memory and dumped at the end,
    so the JSON output does not need the whole export resident to be written.

    Returns the number of records written.
    """
    records_written = 0

    with open(output_path, "w") as fh:
        fh.write("{\n")

        for batch in tqdm(record_batches, desc="Writing JSON"):
            for accession in sorted(batch):
                record = batch[accession]

                serialisable = {}
                for name in COLUMN_ORDER:
                    if name not in include:
                        continue
                    value = record.get(name)
                    # sets are not JSON serialisable
                    serialisable[name] = (
                        sorted(str(v) for v in value)
                        if isinstance(value, (set, frozenset))
                        else value
                    )

                if records_written:
                    fh.write(",\n")
                fh.write(f"  {json.dumps(accession)}: {json.dumps(serialisable)}")
                records_written += 1

        fh.write("\n}\n")

    return records_written
