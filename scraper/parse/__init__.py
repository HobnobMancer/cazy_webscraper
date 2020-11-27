#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
"""Module for parsing the webpages from CAZy."""


import json
import os
import re
import sys
import time

import pandas as pd

from datetime import datetime

from Bio import Entrez, SeqIO
from Bio.PDB import PDBList
from tqdm import tqdm

from scraper.file_io import make_output_directory, write_out_df


def proteins_to_dataframe(families, args, logger):
    """Write protein members of CAZy families to a pandas dataframe.

    Duplicate rows are removed. These rows must be identical across the entire row, not only share
    the same protein name. Duplicates are mostly likely to arise when scraping families and
    subfamilies becuase proteins under the subfamily are also catalogued under the parent family.

    :param families: list, list of CAZy families containing protein members
    :param args: args parser
    :param logger: logger object

    Return nothing.
    """
    # Build empty dictionary to store GenBank accession number synonyms
    all_genbank_synonyms = {}

    # Build empty dataframe to add proteins to
    protein_dataframe = pd.DataFrame(
        data={},
        columns=[
            "Protein_name",
            "CAZy_family",
            "EC#",
            "Source_organism",
            "GenBank",
            "UniProt",
            "PDB/3D",
        ]
    )

    # add proteins to dataframe
    for family in families:
        proteins = family.get_proteins()
        for protein in tqdm(proteins, desc=f"Writing {family.name} proteins to df"):
            if protein is not None:
                # write out protein to df
                protein_dict = protein.get_protein_dict()
                df = pd.DataFrame(protein_dict)
                protein_dataframe = protein_dataframe.append(df, ignore_index=True)
                # add synonym GenBank accessions to all synonyms accessions
                genbank_synonyms = protein.genbank_synonyms
                if genbank_synonyms is not None:
                    all_genbank_synonyms.update(genbank_synonyms)

    time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

    if args.subfamilies is True:
        subfam = "_incld_subfams"
    else:
        subfam = ""
    # compile dataframe name
    if args.data_split == "family":
        df_name = f"cazy_{families[0].name}{subfam}_{time_stamp}"
        synonyms_dict = f"cazy_genbank_synonyms_{families[0].name}{subfam}_{time_stamp}"
    elif args.data_split == "class":
        df_name = f"cazy_{families[0].cazy_class}{subfam}_{time_stamp}"
        synonyms_dict = f"cazy_genbank_synonyms_{families[0].cazy_class}{subfam}_{time_stamp}"
    else:
        df_name = f"cazy_scrape{subfam}_{time_stamp}"
        synonyms_dict = f"cazy_genbank_synonyms{subfam}_{time_stamp}"

    # Remove duplicate proteins
    # This can arise when scraping subfamilies and families within a class
    # The proteins within the subfamilies will also be listed under their parent family
    protein_dataframe = protein_dataframe.drop_duplicates()

    # write out dataframe
    write_out_df(protein_dataframe, df_name, args.output, logger, args.force)

    # write out dictionary of genbank_synonyms
    if args.output is sys.stdout:
        json.dump(all_genbank_synonyms, sys.stdout)

    else:
        output_path = args.output / f"{synonyms_dict}.json"
        with open(output_path, "w") as fh:
            json.dump(all_genbank_synonyms, fh)

    # Additional parsing and retrieval of data: protein sequences and/or structures
    if (args.genbank is not None) and (args.pdb is not None):
        get_structures_and_sequences(protein_dataframe, df_name, args, logger)

    elif (args.genbank is not None) and (args.pdb is None):
        get_genbank_fasta(protein_dataframe, df_name, args, logger)

    elif (args.genbank is None) and (args.pdb is not None):
        get_pdb_structures(protein_dataframe, df_name, args, logger)

    return


def get_structures_and_sequences(protein_dataframe, df_name, args, logger):
    """Coordinate the retrieval of CAZyme structure and sequences from PDB and GenBank.

    :param protein_dataframe: Pandas dataframe, dataframe containing protein data from CAZy
    :param df_name: str, name of the CAZy dataframe
    :param args: cmd args parser
    :param logger: logger object

    Return nothing.
    """
    logger.info(f"Retrieving structures and sequences, for proteins in {df_name}")

    # set user email address for Entrez
    Entrez.email = args.genbank

    # create directory to write FASTA files to
    if (args.genbank_output is not sys.stdout) and (args.genbank_output != args.output):
        make_output_directory(args.genbank_output, logger, args.force, args.nodelete)
    # create directory to write strucure files to
    if (args.pdb_output is not sys.stdout) and (args.pdb_output != args.output):
        make_output_directory(args.pdb_output, logger, args.force, args.nodelete)

    index = 0
    pdbl = PDBList()
    for index in tqdm(
        range(len(protein_dataframe["Protein_name"])),
        desc="Downloading PDBs and FASTAs",
    ):
        # Retrieve accession from GenBank cell
        df_row = protein_dataframe.iloc[index]
        genbank_accession, cazy_family = get_genbank_accession(df_row, df_name, index, logger)
        pdb_accessions = get_pdb_accessions(df_row, df_name, logger)

        if genbank_accession is not None:
            # get FASTA file from GenBank
            fasta_name = f"{genbank_accession}_{cazy_family}.fasta"

            if args.genbank_output is not sys.stdout:
                fasta_name = args.genbank_output / fasta_name

            download_fasta(genbank_accession, fasta_name, args, logger)

        if pdb_accessions is not None:
            # Get structure from PDB
            if pdb_output is not None:
                for accession in pdb_accessions:
                    pdbl.retrieve_pdb_file(
                        f"{accession}",
                        file_format=args.pdb,
                        pdir=args.pdb_output,
                    )
            
            else:
                if args.output is not sys.stdout:
                    for accession in pdb_accessions:
                        pdbl.retrieve_pdb_file(
                            f"{accession}",
                            file_format=args.pdb,
                            pdir=args.output,
                        )
                else:
                    for accession in pdb_accessions:
                        pdbl.retrieve_pdb_file(
                            f"{accession}",
                            file_format=args.pdb,
                            pdir=os.getcwd,
                        )

    return


def get_genbank_accession(df_row, df_name, row_index, logger):
    """Retrieve GenBank accession for protein in the dataframe (df) row.

    Retrieve the first GenBank accession becuase if multiple are given, CAZy only links
    to the first GenBank accession.

    :param df_row: Pandas series, from protein protein dataframe
    :param df_name: str, name of the Pandas dataframe from which row was retrieved
    :param row_index: int, index of the row in the protein dataframe
    :param logger: logger object

    Return two strings, GenBank accession of the protein and the protein's CAZy family.
    """
    family = df_row[1]
    genbank_cell = df_row[4]

    # Retrieve the first accession number, and separate from assoicated HTML link
    find_result = genbank_cell.find(" ")  # finds space separating the first accession and html link
    first_accession = genbank_cell[:find_result]

    # check result is in expected genbank format
    try:
        re.match(r"\D{3}\d+.\d+", first_accession)
    except AttributeError:
        logger.warning(
            f"Could not return accession for protein in row {row_index} in\n"
            "{df_name}.\n"
            "Not retrieving FASTA file for this protein."
        )
        return None, None

    return first_accession, family


def get_pdb_accessions(df_row, df_name, logger):
    """Retrieve PDB accession for protein in the dataframe (df) row.

    Retrieves all PDB accessions becuase multiple are accession can be given. Duplicate accession
    are not retrieved.

    :param df_row: Pandas series, from protein protein dataframe
    :param df_name: str, name of the Pandas dataframe from which row was retrieved
    :param row_index: int, index of the row in the protein dataframe
    :param logger: logger object

    Return list of PDB accessions.
    """
    # Check if null value is stored in pdb cell
    if type(pdb_cell) is float:
        return None

    # Separate the PDB accession
    pdb_accessions = pdb_cell.split(",\n")
    if len(pdb_accessions) == 0:
        return None

    index = 0
    for index in range(len(pdb_accessions)):
        # remove html links
        pdb_accessions[index] = pdb_accessions[index][:pdb_accessions[index].find(" ")]
        # remove additional data in square brackest
        pdb_accessions[index] = pdb_accessions[index].split("[")[0]

    # remove duplicate accessions
    pdb_accessions = list(dict.fromkeys(pdb_accessions))

    return pdb_accessions


def get_pdb_structures(protein_dataframe, df_name, args, logger):
    """Coordinate the retrieval of CAZyme structures from PDB.

    :param protein_dataframe: Pandas dataframe, dataframe containing protein data from CAZy
    :param df_name: str, name of the CAZy dataframe
    :param args: cmd args parser
    :param logger: logger object

    Return nothing.
    """
    logger.info(f"Retrieve PDB structure files for proteins in {df_name}")

    if (args.pdb_output is not sys.stdout) and (args.pdb_output != args.output):
        make_output_directory(args.pdb_output, logger, args.force, args.nodelete)

    # build PDBList object
    pdbl = PDBList()
    index = 0
    for index in tqdm(
        range(len(protein_dataframe["Protein_name"])),
        desc="Downloading structures from RSCB",
    ):
        # Retrieve accession from GenBank cell
        df_row = protein_dataframe.iloc[index]
        pdb_accessions = get_pdb_accessions(df_row, df_name, logger)

        if pdb_accessions is None:
            continue

        # Get download structures from PDB
        for accession in pdb_accessions:
            pdbl.retrieve_pdb_file(f"{accession}", file_format=args.pdb, pdir=args.pdb_output)

    return


def get_genbank_fasta(dataframe, df_name, args, logger):
    """Retrieve GenBank accessions, coordinate downloading associated GenBank FASTA file.

    :param df_row: Pandas series, from protein protein dataframe
    :param df_name: str, name of the Pandas dataframe from which row was retrieved
    :param row_index: int, index of the row in the protein dataframe
    :param logger: logger object0345 454 111

    Return  nothing.
    """
    logger.info(f"Retrieving FASTA files from GenBank, for proteins in {df_name}")

    # set user email address for Entrez
    Entrez.email = args.genbank

    # create directory to write FASTA files to
    if (args.genbank_output is not sys.stdout) and (args.genbank_output != args.output):
        make_output_directory(args.genbank_output, logger, args.force, args.nodelete)

    index = 0
    for index in tqdm(range(len(dataframe["Protein_name"])), desc="Downloading GenBank FASTAs"):
        # Retrieve accession from GenBank cell
        df_row = dataframe.iloc[index]
        accession, cazy_family = get_genbank_accession(df_row, df_name, index, logger)

        if accession is None:
            continue

        # create file name
        file_name = f"{accession}_{cazy_family}.fasta"
        if args.genbank_output is not sys.stdout:
            file_name = args.genbank_output / file_name

        download_fasta(accession, file_name, args, logger)

    return


def download_fasta(accession, file_name, args, logger):
    """Download FASTA file from GenBank FTP server.

    :param accession: str, accession of protein
    :param file_name: str, path to output file desintation
    :param args: cmd args parser
    :param logger: logger object

    Returning nothing.
    """
    if args.genbank_output is not sys.stdout:
        if file_name.exists():
            logger.warning(f"FASTA file {file_name} aleady exists, not downloading again")
            return

    handle = entrez_retry(
        logger,
        Entrez.efetch,
        db="protein",
        id=accession,
        rettype="fasta",
        retmode="text",
    )

    if handle is None:
        logger.warning(f"Failed to download FASTA file for {accession}")
        return

    for record in SeqIO.parse(handle, "fasta"):
        if args.genbank_output is not sys.stdout:
            with open(file_name, "w") as fh:
                SeqIO.write(record, fh, "fasta")
        else:
            SeqIO.write(record, sys.stdout, "fasta")

    return


def entrez_retry(logger, entrez_func, *func_args, **func_kwargs):
    """Call to NCBI using Entrez.
    Maximum number of retries is 10, retry initated when network error encountered.
    :param logger: logger object
    :param retries: parser argument, maximum number of retries excepted if network error encountered
    :param entrez_func: function, call method to NCBI
    :param *func_args: tuple, arguments passed to Entrez function
    :param ** func_kwargs: dictionary, keyword arguments passed to Entrez function
    Returns record.
    """
    record, retries, tries = None, 10, 0

    while record is None and tries < retries:
        try:
            record = entrez_func(*func_args, **func_kwargs)

        except IOError:
            # log retry attempt
            if tries < retries:
                logger.warning(
                    f"Network error encountered during try no.{tries}.\nRetrying in 10s",
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered too many times. Exiting attempt to call to NCBI"
        )
        return

    return record
