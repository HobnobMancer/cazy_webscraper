#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author : Emma E. M. Hobbs
#
# Contact:
# eemh1@st-andrews.ac.uk
#
# Emma Hobbs,
# School of Biology,
# University of St Andrews,
# Biomedical Sciences Research Complex,
# St Andrews,
# Fife,
# KY16 9ST
# Scotland,
# UK
#
# MIT License

import setuptools

from pathlib import Path


# get long description from README.md
with Path("README.md").open("r") as long_description_handle:
    long_description = long_description_handle.read()


setuptools.setup(
    name="cazy_webscraper",
    version="0.1.2",
    # Metadata
    author="Emma E. M. Hobbs",
    author_email="eemh1@st-andrews.ac.uk",
    description="".join(
        [
            (
                "cazy_webscraper provides a webscraper to automate "
                "the retrieval of protein data from the CAZy database, "
                "found at http://www.cazy.org"
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
            "cazy_webscraper.py = scraper.cazy_webscraper:main",
            "get_genbank_sequences.py = scraper.expand.get_genbank_sequences:main",
            "get_pdb_structures.py = scraper.expand.get_pdb_structures:main",
        ]
    },
    install_requires=[
        "biopython>=1.76",
        "mechanicalsoup",
        "pandas>=1.0.3",
        "pyyaml",
        "requests",
        "sqlalchemy==1.3.20",
        "tqdm",
    ],
    packages=setuptools.find_packages(),
    package_data={
        "Conda microenvironment": ["environment.yml"],
        "CAZy dictionary": ["cazy_dictionary.json"],
        "Configuration file": ["scraper_config.yaml"],
    },
    include_package_data=True,
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Licence :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific :: Bioinformatics",
    ],
)
