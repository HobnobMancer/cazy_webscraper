#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2026
# (c) University of Strathclyde 2026
# (c) James Hutton Institute 2026
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# ehobbs@ebi.ac.uk
#
# The MIT License
"""Unit tests for the get_ncbi_taxs subcommand: src.databases.ncbi.taxonomies.

Network access (Bio.Entrez) is always mocked. Where possible, real XML captured from NCBI
(tests/test_inputs/test_inputs_ncbi_tax/) is fed through Bio.Entrez's own parser rather than
hand-built dicts, so the tests exercise the real XML-to-dict shape, not an assumed one.
"""


import sqlite3

from argparse import Namespace
from unittest.mock import patch

import pytest

from src.databases.ncbi.taxonomies import NcbiProtein, get_linked_proteins_batch, get_ncbi_lineages
from src.sql.interface.temp_tables import create_temp_tax_tables


def _args(**overrides):
    defaults = {"retries": 2}
    defaults.update(overrides)
    return Namespace(**defaults)


@pytest.fixture
def ncbi_tax_fixture_dir(test_input_dir):
    return test_input_dir / "test_inputs_ncbi_tax"


class TestNcbiProtein:
    def test_extracts_accession_genus_species_kingdom_from_a_genbank_record(self):
        record = {
            "GBSeq_accession-version": "WP_012345.1",
            "GBSeq_organism": "Trichoderma achlamydosporum",
            "GBSeq_taxonomy": "Eukaryota; Fungi; Ascomycota",
        }

        protein = NcbiProtein()
        protein.get_ncbi_data(record)

        assert protein.protein_id == "WP_012345.1"
        assert protein.genus == "Trichoderma"
        assert protein.species == "achlamydosporum"
        assert protein.kingdom == "Eukaryota"

    def test_multi_word_species_name(self):
        record = {
            "GBSeq_accession-version": "WP_099999.1",
            "GBSeq_organism": "Bacillus subtilis subsp. spizizenii",
            "GBSeq_taxonomy": "Bacteria; Bacillota",
        }

        protein = NcbiProtein()
        protein.get_ncbi_data(record)

        assert protein.genus == "Bacillus"
        assert protein.species == "subtilis subsp. spizizenii"


class TestGetNcbiLineages:
    @patch("src.databases.ncbi.taxonomies.Entrez.efetch")
    @patch("src.databases.ncbi.taxonomies.Entrez.epost")
    def test_parses_a_real_captured_lineage_record(self, mock_epost, mock_efetch, db_path, ncbi_tax_fixture_dir):
        mock_epost.return_value = open(ncbi_tax_fixture_dir / "epost.xml", "rb")
        mock_efetch.return_value = open(ncbi_tax_fixture_dir / "efetchTaxLineage.xml", "rb")

        conn = sqlite3.connect(db_path)
        create_temp_tax_tables(conn)
        query = """
            INSERT INTO NcbiTaxs
            (domain, kingdom, phylum, tax_class, tax_order, family, genus, species, strain, ncbi_tax_id)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """

        err_occurred = get_ncbi_lineages(["2700054"], query, conn, _args())

        # the function falls through to an implicit `return None` on success rather than
        # `return False`, but its only caller (get_ncbi_taxa) just does `if err_occurred:`,
        # so both are equivalent in practice - assert falsy rather than `is False`
        assert not err_occurred
        row = conn.execute(
            "SELECT domain, kingdom, phylum, tax_class, tax_order, family, genus, species, strain, ncbi_tax_id "
            "FROM NcbiTaxs"
        ).fetchone()
        assert row == (
            None,  # no rank in this lineage is literally named "domain" - Eukaryota is a superkingdom
            "Fungi", "Ascomycota", "Sordariomycetes", "Hypocreales", "Hypocreaceae",
            "Trichoderma", "achlamydosporum", None, 2700054,
        )
        conn.close()

    @patch("src.databases.ncbi.taxonomies.Entrez.epost")
    def test_epost_connection_failure_is_reported_without_raising(self, mock_epost, db_path):
        mock_epost.side_effect = RuntimeError("network down")

        conn = sqlite3.connect(db_path)
        create_temp_tax_tables(conn)

        err_occurred = get_ncbi_lineages(["123"], "INSERT INTO NcbiTaxs (ncbi_tax_id) VALUES (?)", conn, _args())

        assert err_occurred is True


class TestGetLinkedProteinsBatch:
    @patch("src.databases.ncbi.taxonomies.Entrez.read")
    @patch("src.databases.ncbi.taxonomies.Entrez.elink")
    def test_maps_taxonomy_id_to_protein_accession_via_the_prot_id_lookup(self, mock_elink, mock_read):
        mock_read.return_value = [
            {
                "IdList": ["2810347"],
                "LinkSetDb": [{"Link": [{"Id": "1995578961"}]}],
            }
        ]

        linked, err_occurred = get_linked_proteins_batch(
            ["2810347"], {"1995578961": "WP_012345.1"}, _args(),
        )

        assert err_occurred is False
        assert linked == [("2810347", "WP_012345.1")]

    @patch("src.databases.ncbi.taxonomies.Entrez.read")
    @patch("src.databases.ncbi.taxonomies.Entrez.elink")
    def test_protein_ids_with_no_local_match_are_silently_skipped(self, mock_elink, mock_read):
        mock_read.return_value = [
            {"IdList": ["2810347"], "LinkSetDb": [{"Link": [{"Id": "UNKNOWN_PROT_ID"}]}]}
        ]

        linked, err_occurred = get_linked_proteins_batch(["2810347"], {}, _args())

        assert linked == []
        assert err_occurred is False

    @patch("src.databases.ncbi.taxonomies.Entrez.elink")
    def test_connection_failure_is_reported_without_raising(self, mock_elink):
        mock_elink.side_effect = RuntimeError("network down")

        linked, err_occurred = get_linked_proteins_batch(["2810347"], {}, _args())

        assert err_occurred is True
        assert linked == []
