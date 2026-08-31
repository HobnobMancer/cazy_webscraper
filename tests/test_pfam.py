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
"""Unit tests for the get_pfams subcommand: src.databases.interpro (the InterPro API client),
src.sql.sql_orm's ProteinPfam model, src.sql.interface.get_data/add_data.*pfam_data, and the
UniProt-accession selector added to src.sql.interface.get_data.get_proteins for it.

Self-contained: builds its own v3-schema db per test via tmp_path, independent of tests/conftest.py
(which targets the unrelated v2 code). Network access is always mocked.
"""


import sqlite3

from argparse import Namespace
from unittest.mock import MagicMock, patch

import pytest
from requests.exceptions import ConnectionError as ReqConnectionError

from src.databases import interpro
from src.databases.interpro import PfamMatchRecord
from src.sql.interface.add_data.add_pfam_data import add_pfam_matches, merge_temp_pfam_relationships
from src.sql.interface.get_data.get_pfam_data import get_db_pfams
from src.sql.interface.get_data.get_proteins import get_ncbi_acc_for_uniprot_acc, get_uniprot_accessions
from src.sql.interface.temp_tables import create_temp_pfam_protein_table, drop_temp_pfam_protein_table


NO_FILTER = set()
NO_TAXONOMY_FILTER = {"genera": set(), "species": set(), "strains": set()}


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


def _add_protein(db_path, protein_accession, uniprot_accession=None):
    """Insert a minimal protein row, optionally with a linked UniProt record."""
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    uniprot_id = None
    if uniprot_accession:
        cur.execute("INSERT INTO Uniprots (uniprot_accession) VALUES (?)", (uniprot_accession,))
        uniprot_id = cur.lastrowid

    cur.execute(
        "INSERT INTO Proteins (protein_accession, source, uniprot_id) VALUES (?, 'ncbi', ?)",
        (protein_accession, uniprot_id),
    )
    protein_id = cur.lastrowid

    conn.commit()
    conn.close()
    return protein_id


# --- parsing the raw InterPro response ---


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


# --- talking to the InterPro API (network mocked) ---


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


# --- the schema Pfam data is written into ---


class TestProteinPfamSchema:
    def test_proteins_pfams_has_its_own_primary_key_and_a_composite_unique_constraint(self, db_path):
        # a plain many-to-many linker table (like Proteins_Ecs) could not represent a protein
        # matching the same Pfam family twice at different positions - this must be an
        # association object instead, with match_start/match_end as part of the identity
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute("SELECT sql FROM sqlite_master WHERE name = 'Proteins_Pfams'")
        schema_sql = cur.fetchone()[0]
        conn.close()

        assert "protein_pfam_id INTEGER NOT NULL" in schema_sql
        assert "PRIMARY KEY (protein_pfam_id)" in schema_sql
        assert "UNIQUE (protein_id, pfam_id, match_start, match_end)" in schema_sql


# --- the get_data / add_data SQL layer ---


class TestGetDbPfams:
    def test_empty_db_returns_empty_map(self, db_path):
        conn = sqlite3.connect(db_path)
        assert get_db_pfams(conn) == {}
        conn.close()

    def test_maps_accession_to_pfam_id(self, db_path):
        conn = sqlite3.connect(db_path)
        conn.execute(
            "INSERT INTO Pfams (accession, name, annotation_type, release) "
            "VALUES ('PF00051', 'Kringle domain', 'domain', '109.0')"
        )
        conn.commit()

        pfam2id = get_db_pfams(conn)

        assert list(pfam2id) == ["PF00051"]
        assert isinstance(pfam2id["PF00051"], int)
        conn.close()


class TestAddPfamMatches:
    def _record(self, protein_id, accession="PF00051", start=1, end=10, interpro_acc="IPR000001"):
        record = PfamMatchRecord(
            "P00734", accession, "Kringle domain", "domain", interpro_acc, start, end, "109.0",
        )
        record.protein_id = protein_id
        return record

    def test_new_pfam_family_and_match_are_staged(self, db_path):
        protein_id = _add_protein(db_path, "TEST_ACC", "P00734")
        conn = sqlite3.connect(db_path)

        new_families, new_matches = add_pfam_matches([self._record(protein_id)], conn)

        assert new_families == 1
        assert new_matches == 1

        cur = conn.cursor()
        cur.execute("SELECT accession FROM Pfams")
        assert cur.fetchall() == [("PF00051",)]

        # staged in the temp table, not yet merged into Proteins_Pfams
        cur.execute("SELECT protein_id, pfam_id, interpro_accession, match_start, match_end FROM TEMP_PFAM_PROTEIN")
        row = cur.fetchone()
        assert row[0] == protein_id
        assert row[2:] == ("IPR000001", 1, 10)
        conn.close()

    def test_reuses_an_existing_pfam_family_rather_than_duplicating_it(self, db_path):
        protein_id = _add_protein(db_path, "TEST_ACC", "P00734")
        conn = sqlite3.connect(db_path)
        conn.execute(
            "INSERT INTO Pfams (accession, name, annotation_type, release) "
            "VALUES ('PF00051', 'Kringle domain', 'domain', '109.0')"
        )
        conn.commit()

        new_families, new_matches = add_pfam_matches([self._record(protein_id)], conn)

        assert new_families == 0
        assert new_matches == 1
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Pfams")
        assert cur.fetchone()[0] == 1
        conn.close()

    def test_two_matches_of_the_same_family_at_different_positions_both_survive(self, db_path):
        protein_id = _add_protein(db_path, "TEST_ACC", "P00734")
        conn = sqlite3.connect(db_path)
        records = [
            self._record(protein_id, start=1, end=10),
            self._record(protein_id, start=50, end=60),
        ]

        new_families, new_matches = add_pfam_matches(records, conn)

        assert new_families == 1  # same Pfam family both times
        assert new_matches == 2   # but two distinct match locations
        conn.close()


class TestMergeTempPfamRelationships:
    def test_merges_staged_matches_into_proteins_pfams(self, db_path):
        protein_id = _add_protein(db_path, "TEST_ACC", "P00734")
        conn = sqlite3.connect(db_path)
        create_temp_pfam_protein_table(conn)
        conn.execute(
            "INSERT INTO Pfams (accession, name, annotation_type, release) "
            "VALUES ('PF00051', 'Kringle domain', 'domain', '109.0')"
        )
        pfam_id = conn.execute("SELECT pfam_id FROM Pfams").fetchone()[0]
        conn.execute(
            "INSERT INTO TEMP_PFAM_PROTEIN (protein_id, pfam_id, interpro_accession, match_start, match_end) "
            "VALUES (?, ?, 'IPR000001', 1, 10)",
            (protein_id, pfam_id),
        )
        conn.commit()
        conn.close()

        merge_temp_pfam_relationships(db_path)

        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute("SELECT protein_id, pfam_id, match_start, match_end FROM Proteins_Pfams")
        assert cur.fetchall() == [(protein_id, pfam_id, 1, 10)]
        conn.close()

    def test_is_idempotent_across_repeated_merges(self, db_path):
        protein_id = _add_protein(db_path, "TEST_ACC", "P00734")
        conn = sqlite3.connect(db_path)
        conn.execute(
            "INSERT INTO Pfams (accession, name, annotation_type, release) "
            "VALUES ('PF00051', 'Kringle domain', 'domain', '109.0')"
        )
        pfam_id = conn.execute("SELECT pfam_id FROM Pfams").fetchone()[0]
        conn.commit()
        conn.close()

        for _ in range(2):
            conn = sqlite3.connect(db_path)
            create_temp_pfam_protein_table(conn)
            conn.execute(
                "INSERT INTO TEMP_PFAM_PROTEIN (protein_id, pfam_id, interpro_accession, match_start, match_end) "
                "VALUES (?, ?, 'IPR000001', 1, 10)",
                (protein_id, pfam_id),
            )
            conn.commit()
            conn.close()
            merge_temp_pfam_relationships(db_path)

        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Proteins_Pfams")
        assert cur.fetchone()[0] == 1
        conn.close()


class TestTempPfamProteinTable:
    def test_create_then_drop_leaves_no_table_behind(self, db_path):
        conn = sqlite3.connect(db_path)
        create_temp_pfam_protein_table(conn)

        cur = conn.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE name = 'TEMP_PFAM_PROTEIN'")
        assert cur.fetchone() is not None

        drop_temp_pfam_protein_table(conn)
        cur.execute("SELECT name FROM sqlite_master WHERE name = 'TEMP_PFAM_PROTEIN'")
        assert cur.fetchone() is None
        conn.close()


# --- selecting which proteins to query InterPro for ---


class TestGetUniprotAccessions:
    def test_only_proteins_with_a_uniprot_accession_are_returned(self, db_path):
        with_uniprot = _add_protein(db_path, "HAS_UNIPROT", "P00734")
        _add_protein(db_path, "NO_UNIPROT")  # get_uniprot_data has not run for this one yet

        result = get_uniprot_accessions(
            NO_FILTER, NO_FILTER, NO_FILTER, NO_TAXONOMY_FILTER, NO_FILTER, db_path,
        )

        assert result == [(with_uniprot, "P00734")]

    def test_additional_filter_excludes_proteins_that_already_have_pfam_matches(self, db_path):
        protein_id = _add_protein(db_path, "HAS_UNIPROT", "P00734")
        conn = sqlite3.connect(db_path)
        conn.execute(
            "INSERT INTO Pfams (accession, name, annotation_type, release) "
            "VALUES ('PF00051', 'Kringle domain', 'domain', '109.0')"
        )
        pfam_id = conn.execute("SELECT pfam_id FROM Pfams").fetchone()[0]
        conn.execute(
            "INSERT INTO Proteins_Pfams (protein_id, pfam_id, interpro_accession, match_start, match_end) "
            "VALUES (?, ?, 'IPR000001', 1, 10)",
            (protein_id, pfam_id),
        )
        conn.commit()
        conn.close()

        result = get_uniprot_accessions(
            NO_FILTER, NO_FILTER, NO_FILTER, NO_TAXONOMY_FILTER, NO_FILTER, db_path,
            additional_filter="P.protein_id NOT IN (SELECT protein_id FROM Proteins_Pfams)",
        )

        assert result == []

    def test_protein_accessions_restricts_the_selection(self, db_path):
        wanted = _add_protein(db_path, "WANTED", "P00734")
        _add_protein(db_path, "NOT_WANTED", "Q99999")

        result = get_uniprot_accessions(
            NO_FILTER, NO_FILTER, NO_FILTER, NO_TAXONOMY_FILTER, NO_FILTER, db_path,
            protein_accessions=["WANTED"],
        )

        assert result == [(wanted, "P00734")]


class TestGetNcbiAccForUniprotAcc:
    # regression test: this function referenced `logger` without importing/defining it
    # anywhere in the module, so any call whose result was empty raised NameError instead of
    # returning an empty set - found and fixed while adding Pfam support (2026-08-31)
    def test_no_match_returns_empty_set_rather_than_raising(self, db_path):
        _add_protein(db_path, "SOME_ACC", "P00734")

        result = get_ncbi_acc_for_uniprot_acc(["Q_DOES_NOT_EXIST"], db_path)

        assert result == set()

    def test_returns_the_ncbi_accession_for_a_matching_uniprot_accession(self, db_path):
        _add_protein(db_path, "SOME_ACC", "P00734")

        result = get_ncbi_acc_for_uniprot_acc(["P00734"], db_path)

        assert result == {"SOME_ACC"}


# --- end-to-end: get_pfam_data() driving the whole pipeline (network mocked) ---


@pytest.fixture
def pfam_db(tmp_path, db_path):
    """A fresh v3-schema db with one protein that already has a UniProt accession."""
    protein_id = _add_protein(db_path, "TEST_ACC", "P00734")
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
