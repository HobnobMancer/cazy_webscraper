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
"""Unit tests for the get_gtdb_taxs subcommand: src.databases.gtdb.

Network access (requests) is always mocked. NOTE: the bundled fixture
tests/test_inputs/test_inputs_gtdb/ar53_taxonomy.tsv.gz is stale - it wraps every field in
literal double quotes, but a live pull of the same file (confirmed 2026-08-31 against
https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_taxonomy.tsv.gz) has no
quoting at all. Using the bundled fixture as-is would test parse_lineages_to_db against a
format it doesn't actually need to handle, so the fixtures below are hand-built to match the
confirmed-live unquoted format instead.
"""


import gzip
import sqlite3

from argparse import Namespace
from unittest.mock import MagicMock, patch

from requests.exceptions import RequestException

from src.databases.gtdb import download_gtdb_file, get_gtdb_release, parse_lineages_to_db
from src.sql.interface.temp_tables import create_temp_gtdb_table


def _args(**overrides):
    defaults = {"retries": 2, "timeout": 10}
    defaults.update(overrides)
    return Namespace(**defaults)


LIVE_FORMAT_ROWS = [
    "RS_GCF_964307495.1\td__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;"
    "f__Methanobacteriaceae;g__Methanocatella;s__Methanocatella smithii",
    "GB_GCA_048462655.1\td__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;"
    "f__Methanobacteriaceae;g__Methanocatella;s__Methanocatella smithii",
    # not in accession_to_genome_id - must be skipped, not matched by accident
    "RS_GCF_999999999.1\td__Archaea;p__Other;c__Other;o__Other;f__Other;g__Other;s__Other sp.",
]


def _write_gtdb_file(path, rows):
    with gzip.open(path, "wt") as fh:
        for row in rows:
            fh.write(row + "\n")
    return path


class TestParseLineagesToDb:
    def test_strips_rs_and_gb_prefixes_before_matching(self, db_path, tmp_path):
        conn = sqlite3.connect(db_path)
        create_temp_gtdb_table(conn)
        gtdb_file = _write_gtdb_file(tmp_path / "ar53.tsv.gz", LIVE_FORMAT_ROWS)

        written = parse_lineages_to_db(
            gtdb_file,
            {"GCF_964307495.1": 1, "GCA_048462655.1": 2},
            "v232",
            conn,
        )

        assert written == 2
        rows = conn.execute(
            "SELECT assembly_accession, genome_id, kingdom, phylum, tax_class, tax_order, "
            "family, genus, species, release FROM TEMP_GTDB ORDER BY assembly_accession"
        ).fetchall()
        assert rows == [
            (
                "GCA_048462655.1", 2, "Archaea", "Methanobacteriota", "Methanobacteria",
                "Methanobacteriales", "Methanobacteriaceae", "Methanocatella",
                "Methanocatella smithii", "v232",
            ),
            (
                "GCF_964307495.1", 1, "Archaea", "Methanobacteriota", "Methanobacteria",
                "Methanobacteriales", "Methanobacteriaceae", "Methanocatella",
                "Methanocatella smithii", "v232",
            ),
        ]
        conn.close()

    def test_genomes_not_of_interest_are_skipped(self, db_path, tmp_path):
        conn = sqlite3.connect(db_path)
        create_temp_gtdb_table(conn)
        gtdb_file = _write_gtdb_file(tmp_path / "ar53.tsv.gz", LIVE_FORMAT_ROWS)

        written = parse_lineages_to_db(gtdb_file, {"GCF_964307495.1": 1}, "v232", conn)

        assert written == 1
        conn.close()

    def test_batches_flush_before_the_end_of_the_file(self, db_path, tmp_path):
        # batch_size=1 forces a flush after every row, not just the final flush() call
        conn = sqlite3.connect(db_path)
        create_temp_gtdb_table(conn)
        gtdb_file = _write_gtdb_file(tmp_path / "ar53.tsv.gz", LIVE_FORMAT_ROWS)

        written = parse_lineages_to_db(
            gtdb_file, {"GCF_964307495.1": 1, "GCA_048462655.1": 2}, "v232", conn, batch_size=1,
        )

        assert written == 2
        assert conn.execute("SELECT COUNT(*) FROM TEMP_GTDB").fetchone()[0] == 2
        conn.close()


class TestGetGtdbRelease:
    @patch("src.databases.gtdb.requests.get")
    def test_reads_the_version_from_the_first_line(self, mock_get):
        response = MagicMock(status_code=200, text="v232\n\nReleased Apr 15, 2026\n")
        response.raise_for_status.return_value = None
        mock_get.return_value = response

        assert get_gtdb_release() == "v232"

    @patch("src.databases.gtdb.requests.get")
    def test_connection_failure_returns_none_rather_than_raising(self, mock_get):
        mock_get.side_effect = RequestException("boom")

        assert get_gtdb_release() is None


class TestDownloadGtdbFile:
    @patch("src.databases.gtdb.requests.get")
    def test_already_downloaded_file_is_reused(self, mock_get, tmp_path):
        existing = tmp_path / "ar53_taxonomy.tsv.gz"
        existing.write_bytes(b"already here")

        result = download_gtdb_file("ar53_taxonomy.tsv.gz", existing, _args())

        assert result == existing
        mock_get.assert_not_called()

    @patch("src.databases.gtdb.requests.get")
    def test_streams_the_response_to_disk(self, mock_get, tmp_path):
        response = MagicMock()
        response.raise_for_status.return_value = None
        response.headers = {}
        response.iter_content.return_value = [b"chunk1", b"chunk2"]
        response.__enter__.return_value = response
        response.__exit__.return_value = False
        mock_get.return_value = response

        out_path = tmp_path / "ar53_taxonomy.tsv.gz"
        result = download_gtdb_file("ar53_taxonomy.tsv.gz", out_path, _args())

        assert result == out_path
        assert out_path.read_bytes() == b"chunk1chunk2"

    @patch("src.databases.gtdb.time.sleep")
    @patch("src.databases.gtdb.requests.get")
    def test_no_file_is_left_behind_when_every_attempt_fails(self, mock_get, mock_sleep, tmp_path):
        mock_get.side_effect = RequestException("boom")
        out_path = tmp_path / "ar53_taxonomy.tsv.gz"

        result = download_gtdb_file("ar53_taxonomy.tsv.gz", out_path, _args(retries=2))

        assert result is None
        assert not out_path.exists()
        assert mock_get.call_count == 2
