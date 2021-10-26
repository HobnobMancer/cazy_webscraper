#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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
"""Parse data retrieved from CAZy."""


import logging
import re

from tqdm.notebook import tqdm
from zipfile import ZipFile


def extract_cazy_file_data(cazy_txt_path):
    """Retrieve data from the text file downloaded from CAZy.
    
    :param cazy_txt_path: Path, location where CAZy text file was downloaded
    
    Return list of lines from the CAZy text file.
    """
    with ZipFile(cazy_txt_path) as zip_handle:
        cazy_filepath = zip_handle.namelist()[0]
    
        with zip_handle.open(cazy_filepath) as fh:
            lines = fh.read().splitlines()
    
    return lines


def parse_cazy_data(
    lines,
    class_filter,
    fam_filter,
    kingdom_filter,
    tax_filter,
    cazy_fam_populations,
):
    """Extract GenBank accession, taxonomy data and CAZy (sub)family annotations from CAZy txt file.
    
    :param lines: list of str, lines from CAZy txt file, one unqiue line is one item in the list
    :param class_filter: set of CAZy class to scrape
    :param fam_filter: set of CAZy families to scrape
    :param kingdom_filter: set of tax Kingdoms to limit scrape to
    :param tax_filter: set of tax (genus, species, strains) filters to limit scrape to
    :param cazy_fam_populations: None if args.validate is False, or dict of CAZy listed CAZy 
        (sub)fam population sizes if args.validate is True
    
    Return cazy_data: dict, {gbk: {"organism": set(str), "families": {'fam': set (subfam)}}}
    """
    logger = logging.getLogger(__name__)

    # define dicts CAZy data will be stored in
    cazy_data = {} # {genbank_accession: {organism,}, families: {(fam, subfam,),} }
    taxa_data = {} # {kingdom: organism}
    
    # used for verbose logging
    gbk_accessions = set()
    fam_annotations = set()
    organisms = set()
    kingdoms = set()

    for line in tqdm(lines, 'Parsing CAZy txt file'):
        line_str = line.decode()
        line_data = line_str.split('\t')

        # retrieve CAZyme data
        cazy_fam = line_data[0]
        if cazy_fam.find("_") != -1:
            cazy_subfam = cazy_fam
            cazy_fam = cazy_fam[:cazy_fam.find("_")]
        else:
            cazy_subfam = None
        fam_annotations.add( (cazy_fam, cazy_subfam) )
       
        kingdom = line_data[1]
        kingdoms.add(kingdom)
        
        organism = line_data[2]
        organisms.add(organism)
        
        gbk_accession = line_data[3]
        gbk_accessions.add(gbk_accession)

        # Apply filters
        if (len(class_filter) == 0) and (len(fam_filter) == 0):
            cazy_data = apply_kingdom_tax_filters(
                cazy_data,
                kingdom_filter,
                tax_filter,
                gbk_accession,
                cazy_fam,
                cazy_subfam,
                organism,
                kingdom,
            )

        # check if scraping the class the family belongs to
        if len(class_filter) != 0:
            fam_class = re.match(r"\D{2,3}", cazy_fam).group()
            
            if fam_class in class_filter:
                cazy_data = apply_kingdom_tax_filters(
                    cazy_data,
                    kingdom_filter,
                    tax_filter,
                    gbk_accession,
                    cazy_fam,
                    cazy_subfam,
                    organism,
                    kingdom,
                )

                continue
            
        if (cazy_fam in fam_filter) or (cazy_subfam in fam_filter):
            cazy_data = apply_kingdom_tax_filters(
                cazy_data,
                kingdom_filter,
                tax_filter,
                gbk_accession,
                cazy_fam,
                cazy_subfam,
                organism,
                kingdom,
            )

            continue

    if cazy_fam_populations is not None:
        validate_data_retrieval(cazy_data, cazy_fam_populations)
        
    logger.info(
        "CAZy txt file contained:\n"
        f"{len(list(gbk_accessions))} unique GenBank accessions"
        f"{len(list(fam_annotations))} unique Fam-Subfam annotations"
        f"{len(list(organisms))} unique source organisms"
        f"{len(list(kingdoms))} unique tax kingdoms"
    )

    return cazy_data, taxa_data


def apply_kingdom_tax_filters(
    cazy_data,
    kingdom_filter,
    tax_filter,
    gbk_accession,
    cazy_fam,
    cazy_subfam,
    organism,
    kingdom,
):
    """Apply User defined Kingdom and Taxonomy filters to determine if protein is to be added to 
    the local CAZyme database.
    
    :param cazy_data: dict of proteins to add to the local CAZyme database
    :param kingdom_filter: set of kingdoms to limit the addition of proteins to the db to
    :param tax_filter: set of genus, species and strains to limit the addition of protein to the db to
    :param gbk_accession: str, the GenBank accession of the protein
    :param cazy_fam: str, CAZy family annotation
    :param cazy_subfam: str, CAZy subfamily annotation, or None if protein is not a CAZy subfamily
    :param kingdom: str, taxonomy kingdom of the source organism of the protein

    Return nothing.
    """
    if (len(kingdom_filter) == 0) and (len(tax_filter) == 0):  # kingdom and tax filters NOT enabled
        cazy_data = add_protein_to_dict(
            cazy_data,
            gbk_accession,
            cazy_fam,
            cazy_subfam,
            organism,
            kingdom,
        )

    if len(kingdom_filter) != 0:  # user enabled kingdom filter
        if kingdom in kingdom_filter:
            cazy_data = add_protein_to_dict(
                cazy_data,
                gbk_accession,
                cazy_fam,
                cazy_subfam,
                organism,
                kingdom,
            )

    if len(tax_filter) != 0:  # user enabled tax filter
        if any(filter in organism for filter in tax_filter):
            cazy_data = add_protein_to_dict(
                cazy_data,
                gbk_accession,
                cazy_fam,
                cazy_subfam,
                organism,
                kingdom,
            )
    
    return cazy_data


def add_protein_to_dict(cazy_data, gbk_accession, cazy_fam, cazy_subfam, organism, kingdom):
    """Add protein to dict containing all proteins to be added to the local CAZyme database.
    
    :param cazy_data: dict of proteins to add to the local CAZyme database
    :param cazy_fam: str, name of the CAZy family
    :param cazy_subfam: str, name of the CAZy subfamily or None if not in a subfamily
    :param organism: str, scientific name of the source organism
    :param kingdom: str, taxonomy kingdom of the source organism of the protein

    Return dict of CAZy data
    """
    try:
        cazy_data[gbk_accession]
        cazy_data[gbk_accession]["kingdom"].add(kingdom)
        cazy_data[gbk_accession]["organism"].add(organism)
        
        try:
            cazy_data[gbk_accession]["families"][cazy_fam].add(cazy_subfam)
        except KeyError:
            cazy_data[gbk_accession]["families"][cazy_fam] = {cazy_subfam}
        
    except KeyError:
        cazy_data[gbk_accession] = {
            "kingdom": {kingdom},
            "organism": {organism},
            "families": {cazy_fam: {cazy_subfam}},
        }
    
    return cazy_data


def validate_data_retrieval(cazy_data, cazy_fam_populations):
    """Extract GenBank accession, taxonomy data and CAZy (sub)family annotations from CAZy txt file.

    Check the number of retrieved proteins per CAZy family against the family population sizes
    previously retrieved to the CAZy website.
    
    :param lines: list of str, lines from CAZy txt file, one unqiue line is one item in the list
    :param class_filter: set of CAZy class to scrape
    :param fam_filter: set of CAZy families to scrape
    :param kingdom_filter: set of tax Kingdoms to limit scrape to
    :param tax_filter: set of tax (genus, species, strains) filters to limit scrape to
    
    Return nothing
    """
    logger = logging.getLogger(__name__)

    total_retrieved_fam_pops = {}  # {fam: population}
    post_filtering_fam_pops = {}  # {fam: population}


    for gbk_accession in cazy_data:

        for fam in cazy_data[gbk_accession]['family']:
            # if protein has is not listed under the parent family alone, add fam with None subfam value
            if None not in cazy_data[gbk_accession]['family'][fam]:
                cazy_data[gbk_accession]['family'][fam].add(None)
            
            for subfam in cazy_data[gbk_accession]['family'][fam]:
                if subfam is None:
                    try:
                        post_filtering_fam_pops[fam] += 1
                    except KeyError:
                        post_filtering_fam_pops[fam] = 1
                
                else:
                    try:
                        post_filtering_fam_pops[subfam] += 1
                    except KeyError:
                        post_filtering_fam_pops[subfam] = 1
    
    for fam in total_retrieved_fam_pops:
        try:
            post_filtering_pop = post_filtering_fam_pops[fam]
        except KeyError:
            post_filtering_pop = 0
        
        try:
            cazy_pop = cazy_fam_populations[fam]
        except KeyError:
            cazy_pop = 'Not Retrieved'
        
        if cazy_pop == 'Not Retrieved' or cazy_pop == 'Failed Retrieval':
            logger.warning(
                f"Cannot validate data retrieval for {fam} because failed to retrieve\n"
                "family population size from the CAZy website."
                f"CAZy txt file contained {total_retrieved_fam_pops[fam]} proteins from {fam}\n"
                f"{post_filtering_pop} proteins from {fam} added to db after applying filters."
            )
            continue
        
        if total_retrieved_fam_pops[fam] < cazy_pop:
            logger.warning(
                f"Fewer proteins from {fam} found in the CAZy txt file than in CAZy listed population size\n"
                f"CAZy listed {fam} as containing {total_retrieved_fam_pops[fam]} proteins\n"
                f"CAZy txt file contained {total_retrieved_fam_pops[fam]} proteins from {fam}\n"
                f"{post_filtering_pop} proteins from {fam} added to db after applying filters."
            )

        elif total_retrieved_fam_pops[fam] < cazy_pop:
            logger.warning(
                f"More proteins from {fam} found in the CAZy txt file than in CAZy listed population size\n"
                f"CAZy listed {fam} as containing {total_retrieved_fam_pops[fam]} proteins\n"
                f"CAZy txt file contained {total_retrieved_fam_pops[fam]} proteins from {fam}\n"
                f"{post_filtering_pop} proteins from {fam} added to db after applying filters."
            )

        else:
            logger.warning(
                f"The same number of proteins from {fam} found in the CAZy txt file as the CAZy listed population size\n"
                f"CAZy listed {fam} as containing {total_retrieved_fam_pops[fam]} proteins\n"
                f"CAZy txt file contained {total_retrieved_fam_pops[fam]} proteins from {fam}\n"
                f"{post_filtering_pop} proteins from {fam} added to db after applying filters."
            )

    return
