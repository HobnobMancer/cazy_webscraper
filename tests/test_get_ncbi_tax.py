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
"""Tests expand.taxonomy.get_ncbi_taxs.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


from argparse import Namespace, ArgumentParser
from datetime import datetime
from pathlib import Path

import pytest

from saintBioutils.utilities import logger as saint_logger
from sqlalchemy.exc import IntegrityError

from cazy_webscraper.expand.genbank.taxonomy import get_ncbi_taxs
from cazy_webscraper.ncbi.taxonomy import lineage
from cazy_webscraper.sql import sql_interface
from cazy_webscraper.sql.sql_interface.get_data import get_selected_gbks
from cazy_webscraper.sql.sql_interface.add_data import add_ncbi_tax_data
from cazy_webscraper.utilities.parsers import tax_ncbi_parser
from cazy_webscraper.sql.sql_interface.get_data import get_table_dicts
from cazy_webscraper.sql.sql_interface.get_data import get_records


@pytest.fixture
def mock_building_parser(*args, **kwargs):
    parser_args = ArgumentParser(
        prog="get_ncbi_taxonomy.py",
        usage=None,
        description="Get lineages from NCBI",
        conflict_handler="error",
        add_help=True,
    )
    return parser_args


def test_main(
    mock_building_parser,
    mock_return_logger,
    config_dict,
    db_connection,
    monkeypatch,
    test_dir,
):
    """Test main()"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            cache_dir=(test_dir / "test_outputs" / "test_outputs_ncbi"),
            email="dummyemail",
            nodelete_cache=False,
            batch_size=150,
            config=None,
            classes=None,
            database="fake_database_path",
            ec_filter=None,
            force=True,
            families=None,
            genbank_accessions=None,
            genera=None,
            get_pages=True,
            kingdoms=None,
            log=None,
            nodelete=True,
            output=None,
            retries=10,
            sequence=True,
            seq_update=True,
            subfamilies=True,
            species=None,
            strains=None,
            uniprot_accessions=None,
            use_lineage_cache=None,
            verbose=False,
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return

    def mock_connect_existing_db(*args, **kwards):
        return db_connection, None, "cache_dir"

    def mock_get_expansion_configuration(*args, **kwards):
        return config_dict, set(), set(), set(), dict(), set()

    def mock_get_genbank_accessions(*args, **kwards):
        return {1: 1, 2: 2, 3: 3}

    def mock_get_lineage(*args, **kwards):
        return {'tax_id': {'linaege info': 'kingdom', 'proteins': {'local db protein ids'}}}

    def mock_get_ncbi_data(*args, **kwards):
        return (
            {
                1: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                2: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                3: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}}
            },
            {1, 2, 3},
        )

    monkeypatch.setattr(tax_ncbi_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger, "config_logger", mock_return_logger)
    monkeypatch.setattr(get_ncbi_taxs, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr("cazy_webscraper.expand.genbank.taxonomy.get_ncbi_taxs.make_output_directory", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "get_expansion_configuration", mock_get_expansion_configuration)
    monkeypatch.setattr(sql_interface, "log_scrape_in_db", mock_return_none)
    # not using cached lineages
    monkeypatch.setattr(get_ncbi_taxs, "get_db_proteins", mock_get_genbank_accessions)
    monkeypatch.setattr(get_ncbi_taxs, "get_ncbi_ids", mock_get_ncbi_data)
    monkeypatch.setattr(get_ncbi_taxs, "get_lineage_protein_data", mock_get_lineage)
    # mock adding data to the local CAZyme database
    monkeypatch.setattr(get_ncbi_taxs, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "closing_message", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "cache_taxonomy", mock_return_none)

    get_ncbi_taxs.main()


def test_main_using_lineage_cache(
    mock_building_parser,
    mock_return_logger,
    config_dict,
    db_connection,
    monkeypatch,
    test_dir,
):
    """Test main()"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            cache_dir=(test_dir / "test_outputs" / "test_outputs_ncbi"),
            email="dummyemail",
            nodelete_cache=False,
            batch_size=150,
            config=None,
            classes=None,
            database="fake_database_path",
            ec_filter=None,
            force=True,
            families=None,
            genbank_accessions=None,
            genera=None,
            get_pages=True,
            kingdoms=None,
            log=None,
            nodelete=True,
            output=None,
            retries=10,
            sequence=True,
            seq_update=True,
            subfamilies=True,
            species=None,
            strains=None,
            uniprot_accessions=None,
            use_lineage_cache=Path("tests/test_inputs/test_inputs_ncbi_tax/test_lineage_cache.json"),
            verbose=False,
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return

    def mock_connect_existing_db(*args, **kwards):
        return db_connection, None, "cache_dir"

    def mock_get_expansion_configuration(*args, **kwards):
        return config_dict, set(), set(), set(), dict(), set()

    def mock_get_genbank_accessions(*args, **kwards):
        return {1: 1, 2: 2, 3: 3}

    def mock_get_ncbi_data(*args, **kwards):
        return (
            {
                1: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                2: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                3: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}}
            },
            {1, 2, 3},
        )

    def mock_get_lineage(*args, **kwards):
        return {'tax_id': {'linaege info': 'kingdom', 'proteins': {'local db protein ids'}}}

    monkeypatch.setattr(tax_ncbi_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger, "config_logger", mock_return_logger)
    monkeypatch.setattr(get_ncbi_taxs, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr("cazy_webscraper.expand.genbank.taxonomy.get_ncbi_taxs.make_output_directory", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "get_expansion_configuration", mock_get_expansion_configuration)
    monkeypatch.setattr(sql_interface, "log_scrape_in_db", mock_return_none)
    # not using cached lineages
    monkeypatch.setattr(get_ncbi_taxs, "get_db_proteins", mock_get_genbank_accessions)
    monkeypatch.setattr(get_ncbi_taxs, "get_ncbi_ids", mock_get_ncbi_data)
    monkeypatch.setattr(get_ncbi_taxs, "get_lineage_protein_data", mock_get_lineage)
    # mock adding data to the local CAZyme database
    monkeypatch.setattr(get_ncbi_taxs, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "closing_message", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "cache_taxonomy", mock_return_none)

    get_ncbi_taxs.main()


def test_main_using_lineage_cache_fails(
    mock_building_parser,
    mock_return_logger,
    config_dict,
    db_connection,
    monkeypatch,
    test_dir,
):
    """Test main() when cannot find linaege cache file"""

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            cache_dir=(test_dir / "test_outputs" / "test_outputs_ncbi"),
            email="dummyemail",
            nodelete_cache=False,
            batch_size=150,
            config=None,
            classes=None,
            database="fake_database_path",
            ec_filter=None,
            force=True,
            families=None,
            genbank_accessions=None,
            genera=None,
            get_pages=True,
            kingdoms=None,
            log=None,
            nodelete=True,
            output=None,
            retries=10,
            sequence=True,
            seq_update=True,
            subfamilies=True,
            species=None,
            strains=None,
            uniprot_accessions=None,
            use_lineage_cache=Path("tests/test_inputs/test_inputs_ncbi_tax/test_lineage_cache_NOT_EXIST.json"),
            verbose=False,
        )
        return parser

    def mock_return_none(*args, **kwargs):
        return

    def mock_connect_existing_db(*args, **kwards):
        return db_connection, None, "cache_dir"

    def mock_get_expansion_configuration(*args, **kwards):
        return config_dict, set(), set(), set(), dict(), set()

    def mock_get_genbank_accessions(*args, **kwards):
        return {1: 1, 2: 2, 3: 3}

    def mock_get_ncbi_data(*args, **kwards):
        return (
            {
                1: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                2: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}},
                3: {'ec': {1, 2, 3}, 'pdb': {1, 2, 3}}
            },
            {1, 2, 3},
        )

    def mock_get_lineage(*args, **kwards):
        return {'tax_id': {'linaege info': 'kingdom', 'proteins': {'local db protein ids'}}}

    monkeypatch.setattr(tax_ncbi_parser, "build_parser", mock_building_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(saint_logger, "config_logger", mock_return_logger)
    monkeypatch.setattr(get_ncbi_taxs, "connect_existing_db", mock_connect_existing_db)
    monkeypatch.setattr("cazy_webscraper.expand.genbank.taxonomy.get_ncbi_taxs.make_output_directory", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "get_expansion_configuration", mock_get_expansion_configuration)
    monkeypatch.setattr(sql_interface, "log_scrape_in_db", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "cache_taxonomy", mock_return_none)
    # not using cached lineages
    monkeypatch.setattr(get_ncbi_taxs, "get_db_proteins", mock_get_genbank_accessions)
    monkeypatch.setattr(get_ncbi_taxs, "get_ncbi_ids", mock_get_ncbi_data)
    monkeypatch.setattr(get_ncbi_taxs, "get_lineage_protein_data", mock_get_lineage)
    # mock adding data to the local CAZyme database
    monkeypatch.setattr(get_ncbi_taxs, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "add_ncbi_taxonomies", mock_return_none)
    monkeypatch.setattr(add_ncbi_tax_data, "update_genbank_ncbi_tax", mock_return_none)
    monkeypatch.setattr(get_ncbi_taxs, "closing_message", mock_return_none)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_ncbi_taxs.main()
    assert pytest_wrapped_e.type == SystemExit


def test_get_db_proteins_usr_acc(monkeypatch):
    argsdict = {"args": Namespace(
        genbank_accessions="tests/test_inputs/test_inputs_ncbi_tax/test_accs.txt",
        uniprot_accessions="tests/test_inputs/test_inputs_ncbi_tax/test_accs.txt",
        update_gbk=False,
    )}

    def mock_gbk_dict(*args, **kwards):
        return {'gbk_acc': 'id', 'gbk_acc_1': 'id'}

    monkeypatch.setattr(get_ncbi_taxs, "get_ids", mock_gbk_dict)
    monkeypatch.setattr(get_selected_gbks, "get_ids", mock_gbk_dict)

    monkeypatch.setattr(get_ncbi_taxs, "get_gbk_table_dict", mock_gbk_dict)
    monkeypatch.setattr(get_table_dicts, "get_gbk_table_dict", mock_gbk_dict)

    monkeypatch.setattr(get_ncbi_taxs, "get_uniprot_table_dict", mock_gbk_dict)
    monkeypatch.setattr(get_table_dicts, "get_uniprot_table_dict", mock_gbk_dict)

    monkeypatch.setattr(get_ncbi_taxs, "get_user_uniprot_sequences", mock_gbk_dict)
    monkeypatch.setattr(get_records, "get_user_uniprot_sequences", mock_gbk_dict)

    monkeypatch.setattr(get_ncbi_taxs, "get_no_tax_gbk_table_dict", mock_gbk_dict)
    monkeypatch.setattr(get_table_dicts, "get_no_tax_gbk_table_dict", mock_gbk_dict)

    get_ncbi_taxs.get_db_proteins(set(), set(), set(), {}, set(), 'db_connection', argsdict['args'])


def test_get_db_proteins(monkeypatch):
    argsdict = {"args": Namespace(
        genbank_accessions=None,
        uniprot_accessions=None,
        update_gbk=True,
    )}

    def mock_gbk_dict(*args, **kwards):
        return {'gbk_acc': 'id', 'gbk_acc_1': 'id'}

    monkeypatch.setattr(get_ncbi_taxs, "get_genbank_accessions", mock_gbk_dict)
    monkeypatch.setattr(get_selected_gbks, "get_genbank_accessions", mock_gbk_dict)

    get_ncbi_taxs.get_db_proteins(set(), set(), set(), {}, set(), 'db_connection', argsdict['args'])


def test_get_lineage_fails(monkeypatch):
    argsdict = {"args": Namespace(
        retries=1,
    )}

    def mock_entrez_tax_call(*args, **kwargs):
        """Mocks call to Entrez to retrieve taxonomy record."""
        return None

    monkeypatch.setattr(get_ncbi_taxs, "entrez_retry", mock_entrez_tax_call)

    output = lineage.fetch_lineages(['2700054'], 'QK', 'Webenv', argsdict['args'])
    assert output == None


def test_get_lineage(monkeypatch):
    """Retrieve mock result with https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=2"""
    argsdict = {"args": Namespace(
        retries=1,
    )}

    efetch_result = "tests/test_inputs/test_inputs_ncbi_tax/efetchTaxLineage.xml"

    with open(efetch_result, "rb") as fh:
        result = fh

        def mock_entrez_tax_call(*args, **kwargs):
            """Mocks call to Entrez to retrieve taxonomy record."""
            return result

        monkeypatch.setattr(get_ncbi_taxs, "entrez_retry", mock_entrez_tax_call)

        output = lineage.fetch_lineages(['2700054'], 'queryKey', 'Webenv', argsdict['args'])
        assert output == None


def test_get_linked_proteins(monkeypatch):
    """Get mocked output at https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?db=protein&dbfrom=taxonomy&id=51453&linkname=taxonomy_protein"""
    argsdict = {"args": Namespace(
        retries=10,
    )}

    efetch_result = "tests/test_inputs/test_inputs_ncbi_tax/efetchLinkedProteins.xml"

    with open(efetch_result, "rb") as fh:
        result = fh

        def mock_entrez_tax_call(*args, **kwargs):
            """Mocks call to Entrez to retrieve taxonomy record."""
            return result

        monkeypatch.setattr(get_ncbi_taxs, "entrez_retry", mock_entrez_tax_call)

        output = get_ncbi_taxs.get_tax_proteins(
            '51453',
            {'51453': {1}},
            {
                '2206269991': 'gbk_acc_1',
                '2206269987': 'gbk_acc_2',
            },
            {'gbk_acc_1': 'db_id1', 'gbk_acc_2': 'db_id2'},
            "tests/test_outputs/test_ncbi_tax",
            argsdict['args'],
        )
        assert output == (
                {'51453': {1, 'db_id1', 'db_id2'}},
                True,
            )


def test_get_linked_proteins_fail(monkeypatch):
    """Get mocked output at https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?db=protein&dbfrom=taxonomy&id=51453&linkname=taxonomy_protein"""
    argsdict = {"args": Namespace(
        retries=10,
    )}

    def mock_entrez_tax_call(*args, **kwargs):
        """Mocks call to Entrez to retrieve taxonomy record."""
        return None

    monkeypatch.setattr(get_ncbi_taxs, "entrez_retry", mock_entrez_tax_call)

    output = get_ncbi_taxs.get_tax_proteins(
        '51453',
        {'51453': {1}},
        {
            '2206269991': 'gbk_acc_1',
            '2206269987': 'gbk_acc_2',
        },
        {'gbk_acc_1': 'db_id1', 'gbk_acc_2': 'db_id2'},
        "tests/test_outputs/test_ncbi_tax",
        argsdict['args'],
    )
    assert output == ({'51453': {1}}, False)


def test_link_prot_to_tax(monkeypatch):
    """Get mock entrez output from https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=taxonomy&id=1995578961&linkname=protein_taxonomy"""
    argsdict = {"args": Namespace(
        retries=10,
    )}

    entrez_result = "tests/test_inputs/test_inputs_ncbi_tax/elinkProtTax.xml"

    with open(entrez_result, "rb") as fh:
        result = fh

        def mock_entrez_tax_call(*args, **kwargs):
            """Mocks call to Entrez."""
            return result

        monkeypatch.setattr(get_ncbi_taxs, "entrez_retry", mock_entrez_tax_call)

        output = get_ncbi_taxs.link_prot_taxs(
            'query key',
            'web env',
            argsdict['args'],
            ['accession'],
        )
        
        assert output == {'2810347'}


def test_get_ncbi_ids_both_no_files(monkeypatch):
    argsdict = {"args": Namespace(
        use_tax_ids="",
        use_protein_ids="",
    )}

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_ncbi_taxs.get_ncbi_ids({}, Path("tests/test_outputs/test_ncbi_tax"), argsdict['args'])
    assert pytest_wrapped_e.type == SystemExit


def test_get_ncbi_ids_both_success(monkeypatch):
    argsdict = {"args": Namespace(
        use_tax_ids="tests/test_inputs/test_inputs_ncbi_tax/test_accs.txt",
        use_protein_ids="tests/test_inputs/test_inputs_ncbi_tax/prot_ids.out",
    )}

    output = get_ncbi_taxs.get_ncbi_ids({}, Path("tests/test_outputs/test_ncbi_tax"), argsdict['args'])
    assert output == (['test_gbk', 'test_gbk'], None)


def test_get_ncbi_ids_only_tax(monkeypatch):
    argsdict = {"args": Namespace(
        use_tax_ids="tests/test_inputs/test_inputs_ncbi_tax/test_accs.txt",
        use_protein_ids=None,
    )}

    def mock_get_ncbi_tax_prot_ids(*args, **kwards):
        return ([], {})

    monkeypatch.setattr(get_ncbi_taxs, "get_ncbi_tax_prot_ids", mock_get_ncbi_tax_prot_ids)

    output = get_ncbi_taxs.get_ncbi_ids({}, Path("tests/test_outputs/test_ncbi_tax"), argsdict['args'])
    assert output == (
        ['test_gbk', 'test_gbk'],
        {},
    )


def test_get_ncbi_ids_only_tax_fails(monkeypatch):
    argsdict = {"args": Namespace(
        use_tax_ids="test_FAILs/test_inputs/test_inputs_ncbi_tax/test_accs_FAIL.txt",
        use_protein_ids=None,
    )}

    def mock_get_ncbi_tax_prot_ids(*args, **kwards):
        return ([], {})

    monkeypatch.setattr(get_ncbi_taxs, "get_ncbi_tax_prot_ids", mock_get_ncbi_tax_prot_ids)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_ncbi_taxs.get_ncbi_ids({}, Path("tests/test_outputs/test_ncbi_tax"), argsdict['args'])
    assert pytest_wrapped_e.type == SystemExit


def test_get_ncbi_ids_only_prots(monkeypatch):
    argsdict = {"args": Namespace(
        use_tax_ids=None,
        use_protein_ids="tests/test_inputs/test_inputs_ncbi_tax/prot_ids.out",
    )}

    def mock_get_ncbi_tax_prot_ids(*args, **kwards):
        return ([], {})

    monkeypatch.setattr(get_ncbi_taxs, "get_ncbi_tax_prot_ids", mock_get_ncbi_tax_prot_ids)

    output = get_ncbi_taxs.get_ncbi_ids({}, Path("tests/test_outputs/test_ncbi_tax"), argsdict['args'])
    assert output == (([], None))


def test_get_ncbi_ids_only_prots_fails(monkeypatch):
    argsdict = {"args": Namespace(
        use_tax_ids=None,
        use_protein_ids="tests/test_inputs/test_inputs_ncbi_tax/prot_ids_FAILS.out",
    )}

    def mock_get_ncbi_tax_prot_ids(*args, **kwards):
        return ([], {})

    monkeypatch.setattr(get_ncbi_taxs, "get_ncbi_tax_prot_ids", mock_get_ncbi_tax_prot_ids)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_ncbi_taxs.get_ncbi_ids({}, Path("tests/test_outputs/test_ncbi_tax"), argsdict['args'])
    assert pytest_wrapped_e.type == SystemExit
