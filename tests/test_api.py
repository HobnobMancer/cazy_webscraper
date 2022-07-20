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
"""Tests api mod

These test are intened to be run from the root of the repository using:
pytest -v
"""


import json
import pandas as pd

from argparse import Namespace, ArgumentParser
from datetime import datetime
from pathlib import Path

import pytest

from sqlalchemy.exc import IntegrityError
from saintBioutils.utilities import logger as saint_logger
from saintBioutils.utilities import file_io as saint_fileIO

from cazy_webscraper import cazy_scraper
from cazy_webscraper.api import cw_query_database
from cazy_webscraper.utilities.parsers import api_parser
from cazy_webscraper.utilities import parse_configuration
from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks, get_api_data


@pytest.fixture
def config_dict():
    configuration_dict = {
        'classes': [],
        "Glycoside Hydrolases (GHs)": ["GH3"],
        'GlycosylTransferases (GTs)': [],
        "Polysaccharide Lyases (PLs)": None,
        'Carbohydrate Esterases (CEs)': [],
        'Auxiliary Activities (AAs)': [],
        'Carbohydrate-Binding Modules (CBMs)': [],
    }
    return configuration_dict


@pytest.fixture
def argsdict(db_path):
    data = {
        "args": Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
            cazy_synonyms=None,
            classes=None,
            config=None,
            citation=False,
            db_output=Path(db_path),  ###
            database=None,
            delete_old_relationships=False,
            force=False,
            families=None,
            genera=None,
            kingdoms=None,
            log=None,
            nodelete=False,
            nodelete_cache=False,
            nodelete_log=False,
            retries=10,
            sql_echo=False,
            subfamilies=False,
            species=None,
            strains=None,
            timeout=45,
            validate=True,
            verbose=False,
            version=False,
        )
    }
    return data


@pytest.fixture
def mock_building_parser(*args, **kwargs):
    parser_args = ArgumentParser(
        prog="cw_query_database",
        usage=None,
        description="Interrogate the CAZy database",
        conflict_handler="error",
        add_help=True,
    )
    return parser_args


@pytest.fixture
def argsdict_all():
    data = {
        "args": Namespace(
            email="dummy@domain.com",
            cache_dir=None,
            cazy_data=None,
            cazy_synonyms=None,
            classes=None,
            config=None,
            citation=False,
            db_output=Path("db_path"),
            database=Path("db_path"),
            output_dir=Path("output_dir/output"),
            delete_old_relationships=False,
            force=False,
            families=None,
            genera=None,
            kingdoms=None,
            log=None,
            nodelete=False,
            nodelete_cache=False,
            nodelete_log=False,
            retries=10,
            sql_echo=False,
            subfamilies=False,
            species=None,
            strains=None,
            timeout=45,
            validate=True,
            verbose=False,
            version=False,
            prefix="test",
            include=['class', 'family', 'subfamily', 'kingdom', 'genus', 'organism', 'genbank_seq', 'uniprot_acc', 'uniprot_name', 'ec', 'pdb', 'uniprot_seq']
        )
    }
    return data


@pytest.fixture
def query_data():
    query_data = {
        "genbank_accession": {
            'class': {"CAZy classes"},
            'family': {"CAZy families"},
            'subfamily': {"CAZy subfamilies"},
            'kingdom': "kingdom",
            'genus': "genus",
            'organism': "genus species strain",
            'ec_numbers': {"EC number annotations"},
            'pdb_accessions': {"PDB accessions"},
            'uniprot_accession': "UniProt protein accession",
            'uniprot_name': "Name of the protein from UniProt",
            'uniprot_sequence': "Protein Aa seq from UniProt",
            'uniprot_sequence_date': "Date the seq was last modified in UniProt",
            'gbk_sequence': "Protein Aa seq from GenBank",
            'gbk_sequence_date': "Date the seq was last modified in Gbk",
    },}
    return query_data


def test_main_new_output(config_dict, db_path, mock_config_logger, mock_building_parser, monkeypatch):
    """Test main()"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            file_types=["json", "csv"],
            output_paths="output",
            email="dummy@domain.com",
            cache_dir='cache',
            cazy_data=None,
            cazy_synonyms=None,
            classes=None,
            config=None,
            citation=False,
            db_output=None,
            database=Path(db_path),  ###
            delete_old_relationships=False,
            force=False,
            families=None,
            genera=None,
            kingdoms=None,
            log=Path('log_dir'),
            nodelete=False,
            nodelete_cache=False,
            nodelete_log=False,
            retries=10,
            sql_echo=False,
            subfamilies=False,
            species=None,
            strains=None,
            timeout=45,
            validate=False,
            verbose=False,
            version=False,
            output_dir="tests/test_outputs/test_api",
            prefix="test_",
        )
        return parser

    def mock_parse_config(*args, **kwards):
        return (
            config_dict,
            {},
            {'GH'},
            {'CE1'},
            {'Bacteria'},
            {'genus': 'a genus'},
            {'ec number'},
        )

    def mock_connect_existing_db(*args, **kwards):
        return (
            'connection',
            'logger_name',
            'cache_dir',
        )

    def mock_return_none(*args, **kwargs):
        return

    def mock_query_data(*args, **kwargs):
        return config_dict

    def mock_names(*args, **kwargs):
        return "output_path"

    monkeypatch.setattr(api_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger,  "config_logger", mock_config_logger)
    monkeypatch.setattr(parse_configuration, "get_expansion_configuration", mock_parse_config)
    # monkeypatch.setattr("cazy_webscraper.api.cw_query_database.make_output_directory", mock_return_none)
    monkeypatch.setattr(saint_fileIO, "make_output_directory", mock_return_none)
    monkeypatch.setattr(cazy_scraper, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr(get_selected_gbks, "get_genbank_accessions", mock_return_none)
    monkeypatch.setattr(cw_query_database, "get_query_data", mock_query_data)
    monkeypatch.setattr(cw_query_database, "write_json_output", mock_return_none)
    monkeypatch.setattr(cw_query_database, "write_csv_output", mock_return_none)
    monkeypatch.setattr(cw_query_database, "closing_message", mock_return_none)
    monkeypatch.setattr(cw_query_database, "compile_output_name",mock_names)

    cw_query_database.main()


def test_main_existing_output(config_dict, db_path, mock_config_logger, mock_building_parser, monkeypatch):
    """Test main() when output files already exist"""
    def mock_parser(*args, **kwargs):
        parser = Namespace(
            file_types=["json", "csv"],
            output_paths="output",
            email="dummy@domain.com",
            cache_dir='cache',
            cazy_data=None,
            cazy_synonyms=None,
            classes=None,
            config=None,
            citation=False,
            db_output=None,
            database=Path("db_path"),  ###
            delete_old_relationships=False,
            force=False,
            families=None,
            genera=None,
            kingdoms=None,
            log=Path('log_dir'),
            nodelete=False,
            nodelete_cache=False,
            nodelete_log=False,
            retries=10,
            sql_echo=False,
            subfamilies=False,
            species=None,
            strains=None,
            timeout=45,
            validate=False,
            verbose=False,
            version=False,
            output_dir=Path("tests/test_outputs/test_api"),
            prefix="test_",
            include="families",
            overwrite=False,
        )
        return parser

    def mock_parse_config(*args, **kwards):
        return (
            config_dict,
            {},
            {'GH'},
            {'CE1'},
            {'Bacteria'},
            {'genus': 'a genus'},
            {'ec number'},
        )

    def mock_connect_existing_db(*args, **kwards):
        return (
            'connection',
            'logger_name',
            'cache_dir',
        )

    def mock_return_none(*args, **kwargs):
        return

    def mock_query_data(*args, **kwargs):
        return config_dict

    def mock_names(*args, **kwargs):
        return "tests/test_outputs/test_api/test__db_path_gbkAcc"

    monkeypatch.setattr(api_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger,  "config_logger", mock_config_logger)
    monkeypatch.setattr(parse_configuration, "get_expansion_configuration", mock_parse_config)
    # monkeypatch.setattr("cazy_webscraper.api.cw_query_database.make_output_directory", mock_return_none)
    monkeypatch.setattr(saint_fileIO, "make_output_directory", mock_return_none)
    monkeypatch.setattr(cazy_scraper, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr(get_selected_gbks, "get_genbank_accessions", mock_return_none)
    monkeypatch.setattr(cw_query_database, "get_query_data", mock_query_data)
    monkeypatch.setattr(cw_query_database, "write_json_output", mock_return_none)
    monkeypatch.setattr(cw_query_database, "write_csv_output", mock_return_none)
    monkeypatch.setattr(cw_query_database, "closing_message", mock_return_none)
    monkeypatch.setattr(cw_query_database, "compile_output_name",mock_names)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cw_query_database.main()
    assert pytest_wrapped_e.type == SystemExit


def test_main_existing_output_overwrite(config_dict, db_path, mock_config_logger, mock_building_parser, monkeypatch):
    """Test main() when output files already exist"""
    def mock_parser(*args, **kwargs):
        parser = Namespace(
            file_types=["json", "csv"],
            output_paths="output",
            email="dummy@domain.com",
            cache_dir='cache',
            cazy_data=None,
            cazy_synonyms=None,
            classes=None,
            config=None,
            citation=False,
            db_output=None,
            database=Path("db_path"),  ###
            delete_old_relationships=False,
            force=False,
            families=None,
            genera=None,
            kingdoms=None,
            log=Path('log_dir'),
            nodelete=False,
            nodelete_cache=False,
            nodelete_log=False,
            retries=10,
            sql_echo=False,
            subfamilies=False,
            species=None,
            strains=None,
            timeout=45,
            validate=False,
            verbose=False,
            version=False,
            output_dir=Path("tests/test_outputs/test_api"),
            prefix="test_",
            include="families",
            overwrite=True,
        )
        return parser

    def mock_parse_config(*args, **kwards):
        return (
            config_dict,
            {},
            {'GH'},
            {'CE1'},
            {'Bacteria'},
            {'genus': 'a genus'},
            {'ec number'},
        )

    def mock_connect_existing_db(*args, **kwards):
        return (
            'connection',
            'logger_name',
            'cache_dir',
        )

    def mock_return_none(*args, **kwargs):
        return

    def mock_query_data(*args, **kwargs):
        return config_dict

    def mock_names(*args, **kwargs):
        return "tests/test_outputs/test_api/test__db_path_gbkAcc"

    monkeypatch.setattr(api_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger,  "config_logger", mock_config_logger)
    monkeypatch.setattr(parse_configuration, "get_expansion_configuration", mock_parse_config)
    # monkeypatch.setattr("cazy_webscraper.api.cw_query_database.make_output_directory", mock_return_none)
    monkeypatch.setattr(saint_fileIO, "make_output_directory", mock_return_none)
    monkeypatch.setattr(cazy_scraper, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr(get_selected_gbks, "get_genbank_accessions", mock_return_none)
    monkeypatch.setattr(cw_query_database, "get_query_data", mock_query_data)
    monkeypatch.setattr(cw_query_database, "write_json_output", mock_return_none)
    monkeypatch.setattr(cw_query_database, "write_csv_output", mock_return_none)
    monkeypatch.setattr(cw_query_database, "closing_message", mock_return_none)
    monkeypatch.setattr(cw_query_database, "compile_output_name",mock_names)

    cw_query_database.main()


def test_compile_names(argsdict_all):
    """Test compile_output_name()"""

    assert cw_query_database.compile_output_name(argsdict_all["args"]) == Path("output_dir/output/test_db_path_gbkAcc_classes_fams_subfams_kngdm_genus_orgnsm_gbkSeq_uni_acc_uni_name_ec_pdb_uniprotSeq")


def test_query_data(argsdict_all, monkeypatch):
    """Test get_query_data()"""
    gbk_dict = {"acc": "id"}

    def mock_get_data(*args, **kwards):
        return gbk_dict
    
    monkeypatch.setattr(get_api_data, "get_class_fam_annotations", mock_get_data)
    monkeypatch.setattr(get_api_data, "get_tax_annotations", mock_get_data)
    monkeypatch.setattr(get_api_data, "get_ec_annotations", mock_get_data)
    monkeypatch.setattr(get_api_data, "get_pdb_accessions", mock_get_data)
    monkeypatch.setattr(get_api_data, "get_uniprot_data", mock_get_data)
    monkeypatch.setattr(get_api_data, "get_gbk_seq", mock_get_data)

    cw_query_database.get_query_data(gbk_dict, "db", argsdict_all["args"])


def test_column_names(argsdict_all):
    """Tests get_column_names()"""
    assert cw_query_database.get_column_names(argsdict_all["args"]) == ['genbank_accession', 'class', 'family', 'subfamily', 'kingdom', 'genus', 'source_organism', 'genbank_sequence', 'genbank_sequence_date', 'uniprot_accession', 'uniprot_name', 'ec_number', 'pdb_accession', 'uniprot_sequence', 'uniprot_sequence_date']


def test_single_values():
    """Test add_single_value_to_rows()"""
    new_rows = [[1], [2]]
    key = "family"
    query_data = {"gbk_acc": {key: [1, 2]}}
    
    assert cw_query_database.add_single_value_to_rows(query_data, "gbk_acc", key, new_rows) == [[1, [1, 2]], [2, [1, 2]]]


def test_json_output(query_data, monkeypatch, argsdict_all):
    """Test write_json_output()"""

    def mock_return_none(*args, **kwards):
        return None

    monkeypatch.setattr(json, "dump", mock_return_none)

    cw_query_database.write_json_output("tests/test_outputs/test_api/test_json.json", query_data, argsdict_all["args"])


def test_csv_output(query_data, monkeypatch, argsdict_all):
    """Test write_csv_output()"""

    def mock_columns(*args, **kwards):
        return ['genbank_accession', 'class', 'family', 'subfamily', 'kingdom', 'genus', 'source_organism', 'genbank_sequence', 'genbank_sequence_date', 'uniprot_accession', 'uniprot_name', 'ec_number', 'pdb_accession', 'uniprot_sequence', 'uniprot_sequence_date']
    
    def mock_class(*args, **kwards):
        return [[1,2,3,4,5,6,7,8,9,10,11]]

    def mock_return_none(*args, **kwards):
        return None

    monkeypatch.setattr(cw_query_database, "get_column_names", mock_columns)
    monkeypatch.setattr(cw_query_database, "get_class_fam_relationships", mock_class)
    monkeypatch.setattr(cw_query_database, "add_single_value_to_rows", mock_class)

    cw_query_database.write_csv_output(query_data, argsdict_all["args"], "tests/test_outputs/test_api/test.csv")


def test_get_relationships_all(query_data, argsdict_all):
    """Test get_class_fam_relationships()"""
    cw_query_database.get_class_fam_relationships("genbank_accession", query_data, argsdict_all["args"])


def test_get_relationships_class_only(query_data):
    """Test get_class_fam_relationships() when retrieving only class annotations"""
    data = {
        "args": Namespace(
            include=['class'],
        )
    }
    cw_query_database.get_class_fam_relationships("genbank_accession", query_data, data["args"])


def test_get_relationships_family_only(query_data):
    """Test get_class_fam_relationships() when retrieving only family annotations"""
    data = {
        "args": Namespace(
            include=['family'],
        )
    }
    cw_query_database.get_class_fam_relationships("genbank_accession", query_data, data["args"])


def test_get_relationships_subfamily_only(query_data):
    """Test get_class_fam_relationships() when retrieving only subfamily annotations"""
    data = {
        "args": Namespace(
            include=['subfamily'],
        )
    }
    cw_query_database.get_class_fam_relationships("genbank_accession", query_data, data["args"])


def test_get_relationships_class_fam(query_data):
    """Test get_class_fam_relationships() when retrieving class and family annotations"""
    data = {
        "args": Namespace(
            include=['class', 'family'],
        )
    }
    cw_query_database.get_class_fam_relationships("genbank_accession", query_data, data["args"])


def test_get_relationships_class_subfam(query_data):
    """Test get_class_fam_relationships() when retrieving class and subfamily annotations"""
    data = {
        "args": Namespace(
            include=['class', 'subfamily'],
        )
    }
    cw_query_database.get_class_fam_relationships("genbank_accession", query_data, data["args"])


def test_get_relationships_fam_subfam(query_data):
    """Test get_class_fam_relationships() when retrieving family and subfamily annotations"""
    data = {
        "args": Namespace(
            include=['family','subfamily'],
        )
    }
    cw_query_database.get_class_fam_relationships("genbank_accession", query_data, data["args"])

