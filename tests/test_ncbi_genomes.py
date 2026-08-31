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
"""Unit tests for the get_ncbi_genomes subcommand: src.databases.ncbi.genomes.

Network access (Bio.Entrez and requests) is always mocked. Where a captured fixture exists
(tests/test_inputs/test_inputs_ncbi_genomes/) it is fed through Bio.Entrez's own parser.
"""


import gzip
import sqlite3

from argparse import Namespace
from unittest.mock import MagicMock, patch

import pytest
from requests.exceptions import RequestException

from src.databases.ncbi.genomes import (
    NcbiGenome, download_feature_table, is_refseq_protein, link_proteins_to_assemblies,
)


def _args(**overrides):
    defaults = {"retries": 2}
    defaults.update(overrides)
    return Namespace(**defaults)


class FakeDocSummary(dict):
    """Stands in for the Bio.Entrez DictionaryElement esummary returns: dict-like, plus an
    `.attributes` mapping (NCBI puts the UID there, not as a regular key)."""
    def __init__(self, *args, attributes=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.attributes = attributes or {}


class TestIsRefseqProtein:
    @pytest.mark.parametrize("accession", ["WP_012345678.1", "NP_012345.1", "XP_012345.1", "YP_012345.1"])
    def test_refseq_prefixes(self, accession):
        assert is_refseq_protein(accession) is True

    @pytest.mark.parametrize("accession", ["AAA12345.1", "QRS12345.1", "CAA12345.1"])
    def test_genbank_prefixes(self, accession):
        assert is_refseq_protein(accession) is False


class TestParseGenomeRecord:
    def test_extracts_uid_name_and_both_accessions_from_synonym(self):
        record = FakeDocSummary(
            {"AssemblyName": "ASM123v1", "Synonym": {"Genbank": "GCA_000001.1", "RefSeq": "GCF_000001.1"}},
            attributes={"uid": "12345"},
        )

        genome = NcbiGenome()
        genome.parse_genome_record(record)

        assert genome.genome_id == "12345"
        assert genome.assembly_name == "ASM123v1"
        assert genome.gbk_accession == "GCA_000001.1"
        assert genome.refseq_accession == "GCF_000001.1"

    def test_falls_back_to_assembly_accession_when_synonym_is_missing(self):
        record = FakeDocSummary(
            {"AssemblyName": "ASM123v1", "AssemblyAccession": "GCF_000001.1"},
            attributes={"uid": "12345"},
        )

        genome = NcbiGenome()
        genome.parse_genome_record(record)

        assert genome.refseq_accession == "GCF_000001.1"
        assert genome.gbk_accession is None

    def test_gca_only_fallback_accession_is_treated_as_genbank(self):
        record = FakeDocSummary(
            {"AssemblyName": "ASM123v1", "AssemblyAccession": "GCA_000001.1"},
            attributes={"uid": "12345"},
        )

        genome = NcbiGenome()
        genome.parse_genome_record(record)

        assert genome.gbk_accession == "GCA_000001.1"
        assert genome.refseq_accession is None

    def test_malformed_record_is_handled_without_raising(self):
        genome = NcbiGenome()
        genome.parse_genome_record(None)  # .get()/.attributes on None raises AttributeError
        assert genome.genome_id is None


class TestCompileUrl:
    def test_builds_the_ncbi_ftp_path_from_the_accession_digits(self):
        genome = NcbiGenome()
        genome.assembly_name = "ASM123v1"

        url = genome.compile_url("GCF_000123456.1")

        assert url == (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/123/456/"
            "GCF_000123456.1_ASM123v1/GCF_000123456.1_ASM123v1_feature_table.txt.gz"
        )


class TestParseFeatureTableToDb:
    def _make_feature_table(self, tmp_path, rows):
        """rows: list of (product_accession, non_redundant_refseq) tuples."""
        path = tmp_path / "features.txt.gz"
        with gzip.open(path, "wt") as fh:
            for product_acc, nr_refseq in rows:
                # a real feature table has many columns; only 10 (product_accession) and
                # 11 (non_redundant_refseq), 0-indexed, are read - pad the rest with "-"
                cols = ["-"] * 12
                cols[0] = "CDS"
                cols[10] = product_acc
                cols[11] = nr_refseq
                fh.write("\t".join(cols) + "\n")
            fh.write("\t".join(["gene"] + ["-"] * 11) + "\n")  # non-CDS rows are skipped
        return path

    def test_maps_wanted_proteins_to_the_genome_and_deletes_the_file_afterwards(self, db_path, tmp_path):
        conn = sqlite3.connect(db_path)
        conn.execute("INSERT INTO Proteins (protein_accession, source) VALUES ('WP_000001.1', 'ncbi')")
        protein_id = conn.execute("SELECT protein_id FROM Proteins").fetchone()[0]
        conn.execute(
            "CREATE TABLE TEMP_GENOME2PROTEIN (ncbi_genome_id INTEGER, protein_id INTEGER, "
            "genome_id INTEGER NULL, PRIMARY KEY (ncbi_genome_id, protein_id))"
        )
        conn.commit()

        table_path = self._make_feature_table(tmp_path, [
            ("WP_000001.1", "-"),
            ("UNWANTED.1", "-"),
        ])

        genome = NcbiGenome()
        genome.genome_id = "999"
        genome.gbk_accession = "GCA_000001.1"

        written = genome.parse_feature_table_to_db(table_path, {"WP_000001.1"}, conn)

        assert written == 1
        row = conn.execute("SELECT ncbi_genome_id, protein_id FROM TEMP_GENOME2PROTEIN").fetchone()
        # ncbi_genome_id has INTEGER affinity, so sqlite coerces the "999" string on insert
        assert row == (999, protein_id)
        assert not table_path.exists()  # cleaned up after parsing
        conn.close()

    def test_checks_the_non_redundant_refseq_column_too(self, db_path, tmp_path):
        # some GCA tables populate col 11 with the equivalent WP_ accession instead of col 10
        conn = sqlite3.connect(db_path)
        conn.execute("INSERT INTO Proteins (protein_accession, source) VALUES ('WP_000002.1', 'ncbi')")
        protein_id = conn.execute("SELECT protein_id FROM Proteins").fetchone()[0]
        conn.execute(
            "CREATE TABLE TEMP_GENOME2PROTEIN (ncbi_genome_id INTEGER, protein_id INTEGER, "
            "genome_id INTEGER NULL, PRIMARY KEY (ncbi_genome_id, protein_id))"
        )
        conn.commit()

        table_path = self._make_feature_table(tmp_path, [("GBK_LOCAL_ACC.1", "WP_000002.1")])

        genome = NcbiGenome()
        genome.genome_id = "999"

        written = genome.parse_feature_table_to_db(table_path, {"WP_000002.1"}, conn)

        assert written == 1
        conn.close()

    def test_no_wanted_proteins_writes_nothing(self, db_path):
        conn = sqlite3.connect(db_path)
        genome = NcbiGenome()
        assert genome.parse_feature_table_to_db(None, set(), conn) == 0


class TestLinkProteinsToAssemblies:
    @patch("src.databases.ncbi.genomes.Entrez.efetch")
    def test_hop1_then_hop2_returns_assembly_ids(self, mock_efetch, test_input_dir):
        real_hop1_fixture = test_input_dir / "test_inputs_ncbi_genomes" / "elink_prot_nuccore.xml"
        with patch("src.databases.ncbi.genomes.Entrez.elink") as mock_elink:
            mock_elink.side_effect = [
                open(real_hop1_fixture, "rb"),  # hop 1: protein -> nuccore
                MagicMock(),                     # hop 2: nuccore -> assembly (content faked below)
            ]
            with patch("src.databases.ncbi.genomes.Entrez.read") as mock_read:
                mock_read.side_effect = [
                    # hop 1's real parsed shape (matches the fixture used above)
                    [{"LinkSetDb": [{"Link": [{"Id": "2131947417"}]}]}],
                    # hop 2: nuccore 2131947417 -> assembly 555
                    [{"LinkSetDb": [{"Link": [{"Id": "555"}]}]}],
                ]

                assembly_ids, connection_err = link_proteins_to_assemblies("env", "1", _args())

        assert connection_err is False
        assert assembly_ids == ["555"]

    @patch("src.databases.ncbi.genomes.Entrez.elink")
    @patch("src.databases.ncbi.genomes.Entrez.read")
    def test_no_nuccore_link_short_circuits_before_hop2(self, mock_read, mock_elink):
        mock_read.return_value = [{"LinkSetDb": []}]  # no Link entries at all

        assembly_ids, connection_err = link_proteins_to_assemblies("env", "1", _args())

        assert assembly_ids == []
        assert connection_err is False
        mock_elink.assert_called_once()  # only hop 1 was attempted

    @patch("src.databases.ncbi.genomes.Entrez.elink")
    def test_connection_error_is_reported_without_raising(self, mock_elink):
        mock_elink.side_effect = RuntimeError("network down")

        assembly_ids, connection_err = link_proteins_to_assemblies("env", "1", _args())

        assert connection_err is True
        assert assembly_ids == []


class TestDownloadFeatureTable:
    def _genome(self):
        genome = NcbiGenome()
        genome.assembly_name = "ASM1v1"
        return genome

    @patch("src.databases.ncbi.genomes.requests.head")
    def test_404_head_response_returns_none_without_attempting_download(self, mock_head, tmp_path):
        mock_head.return_value = MagicMock(status_code=404)

        result = download_feature_table(self._genome(), "GCF_000001.1", tmp_path, _args())

        assert result is None

    @patch("src.databases.ncbi.genomes.requests.get")
    @patch("src.databases.ncbi.genomes.requests.head")
    def test_successful_download_writes_the_file(self, mock_head, mock_get, tmp_path):
        mock_head.return_value = MagicMock(status_code=200)
        response = MagicMock()
        response.raise_for_status.return_value = None
        response.iter_content.return_value = [b"chunk1", b"chunk2"]
        mock_get.return_value = response

        result = download_feature_table(self._genome(), "GCF_000001.1", tmp_path, _args())

        assert result == tmp_path / "GCF_000001.1.feature_table.txt.gz"
        assert result.read_bytes() == b"chunk1chunk2"

    @patch("src.databases.ncbi.genomes.requests.head")
    def test_already_downloaded_file_is_reused_not_refetched(self, mock_head, tmp_path):
        mock_head.return_value = MagicMock(status_code=200)
        existing = tmp_path / "GCF_000001.1.feature_table.txt.gz"
        existing.write_bytes(b"already here")

        with patch("src.databases.ncbi.genomes.requests.get") as mock_get:
            result = download_feature_table(self._genome(), "GCF_000001.1", tmp_path, _args())
            mock_get.assert_not_called()

        assert result == existing

    @patch("src.databases.ncbi.genomes.time.sleep")
    @patch("src.databases.ncbi.genomes.requests.get")
    @patch("src.databases.ncbi.genomes.requests.head")
    def test_exhausting_retries_on_download_returns_none(self, mock_head, mock_get, mock_sleep, tmp_path):
        mock_head.return_value = MagicMock(status_code=200)
        mock_get.side_effect = RequestException("boom")

        result = download_feature_table(self._genome(), "GCF_000001.1", tmp_path, _args(retries=2))

        assert result is None
        assert mock_get.call_count == 2
