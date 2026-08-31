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
"""Unit tests for src.databases.interpro (Pfam domain retrieval from InterPro).

These target the v3 code under ``src/`` and are self-contained: no network access, and no
dependency on the (v2-only, unrelated) fixtures in tests/conftest.py.
"""


import sqlite3

from argparse import Namespace
from unittest.mock import MagicMock, patch

import pytest
from requests.exceptions import ConnectionError as ReqConnectionError

from src.databases import interpro
from src.sql import sql_orm


# Trimmed but real payload shape, captured live from
# https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/P00734/
SAMPLE_RESULTS = [
    {
        "metadata": {
            "accession": "PF00051",
            "name": "Kringle domain",
            "source_database": "pfam",
            "type": "domain",
            "integrated": "IPR000001",
        },
        "proteins": [
            {
                "accession": "p00734",
                "entry_protein_locations": [
                    {"fragments": [{"start": 108, "end": 186}]},
                    {"fragments": [{"start": 213, "end": 291}]},
                ],
            }
        ],
    },
    {
        "metadata": {
            "accession": "PF00089",
            "name": "Trypsin",
            "source_database": "pfam",
            "type": "domain",
            "integrated": "IPR001254",
        },
        "proteins": [
            {
                "accession": "p00734",
                "entry_protein_locations": [
                    {"fragments": [{"start": 364, "end": 613}]},
                ],
            }
        ],
    },
]


def _args(**overrides):
    defaults = {"retries": 2, "timeout": 10}
    defaults.update(overrides)
    return Namespace(**defaults)


def _mock_response(status_code=200, json_body=None, headers=None):
    response = MagicMock(status_code=status_code)
    response.json.return_value = json_body or {}
    response.headers = headers or {}
    response.raise_for_status.return_value = None
    return response


class TestParsePfamResults:
    def test_flattens_one_record_per_match_fragment(self):
        records = interpro.parse_pfam_results("P00734", SAMPLE_RESULTS, "109.0")

        assert len(records) == 3
        assert {r.pfam_accession for r in records} == {"PF00051", "PF00089"}
        assert sorted((r.match_start, r.match_end) for r in records) == [
            (108, 186), (213, 291), (364, 613),
        ]

    def test_carries_metadata_onto_each_record(self):
        records = interpro.parse_pfam_results("P00734", SAMPLE_RESULTS, "109.0")
        kringle = next(r for r in records if r.pfam_accession == "PF00051")

        assert kringle.uniprot_acc == "P00734"
        assert kringle.name == "Kringle domain"
        assert kringle.annotation_type == "domain"
        assert kringle.interpro_accession == "IPR000001"
        assert kringle.release == "109.0"
        assert kringle.protein_id is None  # only set later, once mapped to a local protein

    def test_ignores_entries_from_other_member_databases(self):
        results = [{
            "metadata": {
                "accession": "G3DSA:1.2.3", "source_database": "cathgene3d", "type": "domain",
            },
            "proteins": [{"entry_protein_locations": [{"fragments": [{"start": 1, "end": 10}]}]}],
        }]

        assert interpro.parse_pfam_results("P00734", results, "109.0") == []

    def test_empty_results_yields_no_records(self):
        assert interpro.parse_pfam_results("P00734", [], "109.0") == []


class TestFetchPfamMatches:
    @patch("src.databases.interpro.requests.get")
    def test_success_returns_results_and_release(self, mock_get):
        mock_get.return_value = _mock_response(
            json_body={"results": SAMPLE_RESULTS, "next": None},
            headers={"InterPro-Version": "109.0"},
        )

        results, release, connection_err = interpro.fetch_pfam_matches("P00734", _args())

        assert connection_err is False
        assert release == "109.0"
        assert results == SAMPLE_RESULTS
        mock_get.assert_called_once()

    @patch("src.databases.interpro.requests.get")
    def test_204_is_reported_as_no_results_not_a_connection_error(self, mock_get):
        # InterPro returns 204 both for "no Pfam matches" and "unrecognised accession" -
        # confirmed live against the real API - so this must not be treated as a failure.
        mock_get.return_value = _mock_response(status_code=204)

        results, release, connection_err = interpro.fetch_pfam_matches("BOGUS123", _args())

        assert results == []
        assert connection_err is False

    @patch("src.databases.interpro.requests.get")
    def test_follows_cursor_pagination(self, mock_get):
        page1 = _mock_response(json_body={
            "results": [SAMPLE_RESULTS[0]],
            "next": "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/P00734/?cursor=x",
        })
        page2 = _mock_response(json_body={"results": [SAMPLE_RESULTS[1]], "next": None})
        mock_get.side_effect = [page1, page2]

        results, release, connection_err = interpro.fetch_pfam_matches("P00734", _args())

        assert connection_err is False
        assert results == SAMPLE_RESULTS
        assert mock_get.call_count == 2

    @patch("src.databases.interpro.time.sleep")  # skip the real exponential backoff delay
    @patch("src.databases.interpro.requests.get")
    def test_retries_after_a_connection_error_then_succeeds(self, mock_get, mock_sleep):
        mock_get.side_effect = [ReqConnectionError("boom"), _mock_response(json_body={"results": [], "next": None})]

        results, release, connection_err = interpro.fetch_pfam_matches("P00734", _args(retries=3))

        assert connection_err is False
        assert mock_get.call_count == 2
        mock_sleep.assert_called_once()

    @patch("src.databases.interpro.time.sleep")
    @patch("src.databases.interpro.requests.get")
    def test_reports_connection_error_once_retries_are_exhausted(self, mock_get, mock_sleep):
        mock_get.side_effect = ReqConnectionError("boom")

        results, release, connection_err = interpro.fetch_pfam_matches("P00734", _args(retries=2))

        assert connection_err is True
        assert results == []
        assert mock_get.call_count == 2


@pytest.fixture
def pfam_db(tmp_path):
    """A fresh v3-schema db with one protein that already has a UniProt accession."""
    db_path = tmp_path / "test.db"
    sql_orm.get_db_connection(db_path, False, new=True)

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("INSERT INTO Uniprots (uniprot_accession) VALUES ('P00734')")
    uniprot_id = cur.lastrowid
    cur.execute(
        "INSERT INTO Proteins (protein_accession, source, uniprot_id) VALUES ('TEST_ACC', 'ncbi', ?)",
        (uniprot_id,),
    )
    protein_id = cur.lastrowid
    conn.commit()
    conn.close()

    return db_path, protein_id


class TestGetPfamData:
    def test_writes_pfam_families_and_matches_to_the_db(self, tmp_path, pfam_db, monkeypatch):
        db_path, protein_id = pfam_db
        monkeypatch.setattr(
            interpro, "fetch_pfam_matches",
            lambda uniprot_acc, args: (SAMPLE_RESULTS, "109.0", False),
        )

        stats = interpro.get_pfam_data(
            {"P00734": protein_id}, tmp_path, "testrun", _args(database=db_path, batch_size=10),
        )

        assert stats["new pfam families"] == 2
        assert stats["new matches"] == 3
        assert stats["failed accessions"] == 0

        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute("SELECT accession FROM Pfams ORDER BY accession")
        assert [row[0] for row in cur.fetchall()] == ["PF00051", "PF00089"]

        cur.execute(
            "SELECT interpro_accession, match_start, match_end FROM Proteins_Pfams "
            "WHERE protein_id = ? ORDER BY match_start", (protein_id,),
        )
        assert cur.fetchall() == [
            ("IPR000001", 108, 186),
            ("IPR000001", 213, 291),
            ("IPR001254", 364, 613),
        ]

        # the run's scratch table must not be left behind
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'TEMP%'")
        assert cur.fetchall() == []
        conn.close()

    def test_rerunning_does_not_duplicate_rows(self, tmp_path, pfam_db, monkeypatch):
        db_path, protein_id = pfam_db
        monkeypatch.setattr(
            interpro, "fetch_pfam_matches",
            lambda uniprot_acc, args: (SAMPLE_RESULTS, "109.0", False),
        )

        args = _args(database=db_path, batch_size=10)
        interpro.get_pfam_data({"P00734": protein_id}, tmp_path, "run1", args)
        interpro.get_pfam_data({"P00734": protein_id}, tmp_path, "run2", args)

        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Pfams")
        assert cur.fetchone()[0] == 2
        cur.execute("SELECT COUNT(*) FROM Proteins_Pfams")
        assert cur.fetchone()[0] == 3
        conn.close()

    def test_accessions_with_no_matches_are_cached_not_treated_as_failures(self, tmp_path, pfam_db, monkeypatch):
        db_path, protein_id = pfam_db
        monkeypatch.setattr(
            interpro, "fetch_pfam_matches",
            lambda uniprot_acc, args: ([], None, False),
        )

        stats = interpro.get_pfam_data(
            {"P00734": protein_id}, tmp_path, "nomatch", _args(database=db_path, batch_size=10),
        )

        assert stats["accessions with no matches"] == 1
        assert stats["failed accessions"] == 0
        cache_file = tmp_path / "pfam_no_matches_nomatch.txt"
        assert cache_file.exists()
        assert cache_file.read_text().strip() == "P00734"

    def test_failed_accessions_are_cached_separately_from_no_matches(self, tmp_path, pfam_db, monkeypatch):
        db_path, protein_id = pfam_db
        monkeypatch.setattr(
            interpro, "fetch_pfam_matches",
            lambda uniprot_acc, args: ([], None, True),
        )

        stats = interpro.get_pfam_data(
            {"P00734": protein_id}, tmp_path, "failrun", _args(database=db_path, batch_size=10),
        )

        assert stats["failed accessions"] == 1
        assert stats["accessions with no matches"] == 0
        cache_file = tmp_path / "pfam_connection_errors_failrun.txt"
        assert cache_file.exists()
        assert cache_file.read_text().strip() == "P00734"
