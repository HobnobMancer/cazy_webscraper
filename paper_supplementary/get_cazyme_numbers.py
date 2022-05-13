#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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
"""Build a df with the number of CAZymes in CAZy for a list of species"""


import pandas as pd

from argparse import Namespace
from datetime import datetime
from pathlib import Path

from tqdm import tqdm

from cazy_webscraper import (
    connect_existing_db,
    closing_message,
)
from cazy_webscraper.sql.sql_orm import (
    Genbank,
    Taxonomy,
    Session,
)


SPECIES_LIST = "taxonomy_result.txt"
DB_PATH = Path("all_cazy_2022_01_13.db")
OUTPUT_DF = "fungi_cazyme_totals.csv"


def main():
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    argsdict = {
        "args": Namespace(
            database=DB_PATH,
            sql_echo=False,
            verbose=False,
        )
    }

    db_connection, logger_name, cache_dir = connect_existing_db(argsdict["args"], time_stamp, start_time)

    organisms = get_organisms()
    print(f"Parsing {len(organisms)} genera")

    cazyme_totals = get_cazyme_totals(organisms, db_connection)
    print(f"Parsed {len(cazyme_totals)} organisms")

    df = pd.DataFrame(cazyme_totals, columns=["Genus", "Species", "Num_of_CAZymes"])
    df.to_csv(OUTPUT_DF)

    closing_message("cazy_webscraper", start_time, argsdict["args"])


def get_organisms():
    """Get organisms from file and compile into a dict
    
    Separate genus and species
    
    Return dict {genus: {species: num of cazymes (int)}}
    """
    with open(SPECIES_LIST, "r") as fh:
        ncbi_results = fh.read().splitlines()

    organisms = {}
    for hit in tqdm(ncbi_results , desc="Separating genus and species names"):
        genus = hit.split(" ")[0].strip()
        species = hit.replace(genus,"").strip()
        
        try:
            organisms[genus]
            try:
                organisms[genus][species]
            except KeyError:
                organisms[genus][species] = 0
        except KeyError:
            organisms[genus] = {species: 0}

    return organisms


def get_cazyme_totals(organisms, db_connection):
    """Get the number of CAZymes listed in CAZy for each organism
    
    :param organisms: dict {genus: {species: num of cazymes (int)}}
    :param db_connection: open sqlalchemy connection to db engine
    
    Return list of nested lists, one nested list per organism"""
    cazyme_totals = []

    # load add GenBank and taxonomy records
    print("Loading GenBank and Taxnomy records from the local database")
    with Session(bind=db_connection) as session:
        gbk_tax_records = session.query(Genbank, Taxonomy).\
            join(Genbank, (Genbank.taxonomy_id == Taxonomy.taxonomy_id)).\
            all()
    print(f"Loaded {len(gbk_tax_records)} GenBank records")

    # convert organsisms into dict
    for record in tqdm(gbk_tax_records, desc="Getting CAZyme totals"):
        genus = record[1].genus
        species = record[1].species

        try:
            organisms[genus]
            try:
                organisms[genus][species] += 1
            except KeyError:
                pass
        except KeyError:
            pass

    cazyme_totals = []
    genera = list(organisms.keys())
    for genus in tqdm(genera, desc="Reorganising data for the df per genus"):
        genus_species_dict = organisms[genus]

        for species in genus_species_dict:
            cazyme_count = organisms[genus][species]
            cazyme_totals.append([genus, species, cazyme_count])

    return cazyme_totals


if __name__ == "__main__":
    main()
