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
# SOFTWARE.
"""Installation file for pip"""


import setuptools

from pathlib import Path


# get long description from README.md
with Path("README.md").open("r") as long_description_handle:
    long_description = long_description_handle.read()


setuptools.setup(
    name="cazy_webscraper",
    version="2.3.0",
    # Metadata
    author="Emma E. M. Hobbs",
    author_email="eemh1@st-andrews.ac.uk",
    description="".join(
        [
            (
                "A tool to automate retrieving data from CAZy, "
                "build a local CAZyme SQL database, and throughly interrogating the data. "
                "Also, automate retrieving protein data, sequences, EC numbers and structure files "
                "for specific datasets in the CAZyme database from UniProt, GenBank and PDB."
            )
        ]
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords="bioinforamtics protein webscraper",
    platforms="Posix, MacOS X",
    url="https://github.com/HobnobMancer/cazy_webscraper",
    entry_points={
        "console_scripts": [
            "cazy_webscraper = cazy_webscraper.cazy_scraper:main",
            "cw_get_genbank_seqs = cazy_webscraper.expand.genbank.sequences.get_genbank_sequences:main",
            "cw_get_ncbi_taxs = cazy_webscraper.expand.genbank.taxonomy.get_ncbi_taxs:main",
            "cw_get_genomics = cazy_webscraper.expand.genbank.genomes.get_genome_accs:main",
            "cw_get_gtdb_taxs = cazy_webscraper.expand.gtdb.get_gtdb_tax:main",
            "cw_get_uniprot_data = cazy_webscraper.expand.uniprot.get_uniprot_data:main",
            "cw_extract_db_seqs = cazy_webscraper.expand.extract_seqs.extract_db_seqs:main",
            "cw_get_pdb_structures = cazy_webscraper.expand.pdb.get_pdb_structures:main",
            "cw_query_database = cazy_webscraper.api.cw_query_database:main",
            "cw_get_db_schema = cazy_webscraper.sql.get_schema:main",
        ]
    },
    install_requires=[
        "biopython",
        "bioservices>=1.11.0",
        "mechanicalsoup",
        "pandas",
        "pyyaml",
        "requests",
        "saintBioutils>=0.0.25",
        "sqlalchemy==1.4.20",
        "tqdm",
        "html5lib",
    ],
    packages=setuptools.find_packages(),
    package_data={
        "Conda microenvironment": ["environment.yml"],
        "CAZy dictionary": ["cazy_dictionary.json"],
        "Configuration file": ["scraper_config.yaml"],
    },
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
