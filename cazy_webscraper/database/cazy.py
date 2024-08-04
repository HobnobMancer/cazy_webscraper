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
"""Dump the CAZy db dump text file into the local CAZyme db"""


import sqlite3

from pathlib import Path
from zipfile import ZipFile

from tqdm import tqdm


def dump_cazy_txt(cazy_txt_path: Path, db: Path):
    """Dump cazy txt file into the local db

    :param cazy_txt_path: Path to the local tsv file dump of cazy
    :param db: Path to the local CAZyme db
    """
    conn = sqlite3.connect(db)
    cur = conn.cursor()

    with ZipFile(cazy_txt_path) as zip_handle:
        cazy_filepath = zip_handle.namelist()[0]

        with zip_handle.open(cazy_filepath) as fh:
            num_lines = sum(1 for _ in fh)  # Count total lines in the file
            fh.seek(0)  # Reset file pointer to the beginning

            for line_bytes in tqdm(fh, desc="Dumping CAZy data into a temp table", total=num_lines):
                data = line_bytes.decode('utf-8').strip().split()
                # e.g. GH157   Bacteria   Bacteroides cellulosilyticus BFG-250   UBD70155.1    ncbi
                fam, king, genus, protein_id, source = data[0], data[1], data[2], data[-2], data[-1]
                sp = ' '.join([_ for _ in data if _ not in [fam, king, genus, protein_id, source]])
                cur.execute(
                    """
                    INSERT INTO TempTable (family, kingdom, genus, species, protein_id, source)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    (fam, king, genus, sp, protein_id, source)
                )

    conn.commit()
    cur.close()
    conn.close()
