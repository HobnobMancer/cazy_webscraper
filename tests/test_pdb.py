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
"""Unit tests for the get_pdb_structures subcommand: src.databases.pdb.

Network access (requests, Bio.PDB) is always mocked.
"""


import sqlite3

from argparse import Namespace
from unittest.mock import MagicMock, patch

from requests.exceptions import RequestException

from src.databases.pdb import dump_pdb_metadata, download_structure_files, fetch_pdb_metadata, get_pdb_data
from src.sql.interface.temp_tables import create_temp_pdb_structure_table


def _args(**overrides):
    defaults = {"retries": 2, "timeout": 10, "skip_download": True, "file_formats": None}
    defaults.update(overrides)
    return Namespace(**defaults)


def _graphql_response(entries=None, errors=None, status_code=200):
    response = MagicMock(status_code=status_code)
    response.raise_for_status.return_value = None
    response.json.return_value = {"data": {"entries": entries or []}, "errors": errors}
    return response


class TestFetchPdbMetadata:
    @patch("src.databases.pdb.requests.post")
    def test_extracts_method_and_resolution_for_each_entry(self, mock_post):
        mock_post.return_value = _graphql_response(entries=[
            {
                "rcsb_id": "1ABC",
                "exptl": [{"method": "X-RAY DIFFRACTION"}],
                "rcsb_entry_info": {"resolution_combined": [1.5]},
            },
        ])

        metadata, connection_err = fetch_pdb_metadata(["1ABC"], _args())

        assert connection_err is False
        assert metadata == {"1ABC": ("X-RAY DIFFRACTION", 1.5)}

    @patch("src.databases.pdb.requests.post")
    def test_rcsb_omits_accessions_it_does_not_recognise_rather_than_erroring(self, mock_post):
        # only 1ABC is in the response even though 1ABC and BOGUS were both requested -
        # the caller (get_pdb_data) is expected to diff the request against the result
        mock_post.return_value = _graphql_response(entries=[
            {"rcsb_id": "1ABC", "exptl": [{"method": "X-RAY DIFFRACTION"}], "rcsb_entry_info": {}},
        ])

        metadata, connection_err = fetch_pdb_metadata(["1ABC", "BOGUS"], _args())

        assert set(metadata) == {"1ABC"}

    @patch("src.databases.pdb.requests.post")
    def test_method_with_no_resolution_e_g_solution_nmr(self, mock_post):
        mock_post.return_value = _graphql_response(entries=[
            {"rcsb_id": "2XYZ", "exptl": [{"method": "SOLUTION NMR"}], "rcsb_entry_info": {}},
        ])

        metadata, connection_err = fetch_pdb_metadata(["2XYZ"], _args())

        assert metadata == {"2XYZ": ("SOLUTION NMR", None)}

    @patch("src.databases.pdb.requests.post")
    def test_multiple_experimental_methods_are_joined(self, mock_post):
        mock_post.return_value = _graphql_response(entries=[
            {
                "rcsb_id": "3DEF",
                "exptl": [{"method": "X-RAY DIFFRACTION"}, {"method": "NEUTRON DIFFRACTION"}],
                "rcsb_entry_info": {},
            },
        ])

        metadata, connection_err = fetch_pdb_metadata(["3DEF"], _args())

        assert metadata["3DEF"][0] == "X-RAY DIFFRACTION, NEUTRON DIFFRACTION"

    @patch("src.databases.pdb.requests.post")
    def test_graphql_errors_field_is_treated_as_a_connection_error(self, mock_post):
        mock_post.return_value = _graphql_response(errors=[{"message": "bad query"}])

        metadata, connection_err = fetch_pdb_metadata(["1ABC"], _args())

        assert connection_err is True
        assert metadata == {}

    @patch("src.databases.pdb.time.sleep")
    @patch("src.databases.pdb.requests.post")
    def test_exhausting_retries_after_request_exceptions_reports_connection_error(self, mock_post, mock_sleep):
        mock_post.side_effect = RequestException("boom")

        metadata, connection_err = fetch_pdb_metadata(["1ABC"], _args(retries=2))

        assert connection_err is True
        assert mock_post.call_count == 2


class TestDumpPdbMetadata:
    def test_writes_into_the_temp_structure_table(self, db_path):
        conn = sqlite3.connect(db_path)
        create_temp_pdb_structure_table(conn)

        dump_pdb_metadata({"1ABC": ("X-ray", 1.5)}, conn)

        row = conn.execute("SELECT pdb_accession, method, resolution FROM TEMP_PDB_STRUCTURE").fetchone()
        assert row == ("1ABC", "X-ray", 1.5)
        conn.close()

    def test_empty_metadata_is_a_no_op(self, db_path):
        conn = sqlite3.connect(db_path)
        create_temp_pdb_structure_table(conn)
        dump_pdb_metadata({}, conn)  # must not raise
        conn.close()


class TestDownloadStructureFiles:
    def test_skip_download_returns_zero_without_calling_pdblist(self):
        with patch("src.databases.pdb.Bio.PDB.PDBList") as mock_pdblist_cls:
            written = download_structure_files(["1ABC"], _args(skip_download=True, file_formats=["mmCif"]))
            mock_pdblist_cls.assert_not_called()
        assert written == 0

    def test_downloads_one_file_per_format_per_accession(self, tmp_path):
        with patch("src.databases.pdb.Bio.PDB.PDBList") as mock_pdblist_cls:
            instance = mock_pdblist_cls.return_value
            written_file = tmp_path / "1abc.cif"
            written_file.write_text("data")
            instance.retrieve_pdb_file.return_value = str(written_file)

            written = download_structure_files(
                ["1ABC"], _args(skip_download=False, file_formats=["mmCif"], outdir=tmp_path, overwrite=False),
            )

        assert written == 1
        instance.retrieve_pdb_file.assert_called_once_with(
            "1ABC", file_format="mmCif", overwrite=False, pdir=str(tmp_path),
        )

    def test_a_missing_structure_file_does_not_take_the_rest_of_the_batch_with_it(self, tmp_path):
        with patch("src.databases.pdb.Bio.PDB.PDBList") as mock_pdblist_cls:
            instance = mock_pdblist_cls.return_value
            good_file = tmp_path / "2xyz.cif"
            good_file.write_text("data")
            instance.retrieve_pdb_file.side_effect = [None, str(good_file)]

            written = download_structure_files(
                ["1ABC", "2XYZ"], _args(skip_download=False, file_formats=["mmCif"], outdir=tmp_path, overwrite=False),
            )

        assert written == 1

    def test_an_exception_for_one_accession_does_not_abort_the_batch(self, tmp_path):
        with patch("src.databases.pdb.Bio.PDB.PDBList") as mock_pdblist_cls:
            instance = mock_pdblist_cls.return_value
            good_file = tmp_path / "2xyz.cif"
            good_file.write_text("data")
            instance.retrieve_pdb_file.side_effect = [OSError("disk full"), str(good_file)]

            written = download_structure_files(
                ["1ABC", "2XYZ"], _args(skip_download=False, file_formats=["mmCif"], outdir=tmp_path, overwrite=False),
            )

        assert written == 1


class TestGetPdbData:
    def test_retrieves_metadata_and_flags_invalid_accessions(self, db_path, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "src.databases.pdb.fetch_pdb_metadata",
            lambda batch, args: ({"1ABC": ("X-ray", 1.5)}, False),
        )
        monkeypatch.setattr("src.databases.pdb.download_structure_files", lambda accs, args: 0)

        stats = get_pdb_data(
            ["1ABC", "BOGUS"], tmp_path, "testrun",
            _args(database=db_path, batch_size=10, skip_download=True),
        )

        assert stats["metadata retrieved"] == 1
        assert stats["invalid accessions"] == 1
        cache_file = tmp_path / "pdb_invalid_ids_testrun.txt"
        assert cache_file.read_text().strip() == "BOGUS"

        conn = sqlite3.connect(db_path)
        row = conn.execute("SELECT pdb_accession, method, resolution FROM TEMP_PDB_STRUCTURE").fetchone()
        assert row == ("1ABC", "X-ray", 1.5)
        conn.close()

    def test_connection_error_is_cached_and_counted_as_a_failed_batch(self, db_path, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "src.databases.pdb.fetch_pdb_metadata",
            lambda batch, args: ({}, True),
        )

        stats = get_pdb_data(
            ["1ABC"], tmp_path, "failrun", _args(database=db_path, batch_size=10, skip_download=True),
        )

        assert stats["failed batches"] == 1
        assert (tmp_path / "pdb_connection_errors_failrun.txt").read_text().strip() == "1ABC"

    def test_only_accessions_confirmed_present_are_offered_to_the_downloader(self, db_path, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "src.databases.pdb.fetch_pdb_metadata",
            lambda batch, args: ({"1ABC": ("X-ray", 1.5)}, False),
        )
        seen = {}

        def fake_download(accessions, args):
            seen["accessions"] = accessions
            return len(accessions)

        monkeypatch.setattr("src.databases.pdb.download_structure_files", fake_download)

        get_pdb_data(
            ["1ABC", "BOGUS"], tmp_path, "run3",
            _args(database=db_path, batch_size=10, skip_download=False, file_formats=["mmCif"]),
        )

        assert seen["accessions"] == ["1ABC"]
