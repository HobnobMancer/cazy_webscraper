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
"""Unit tests for the download_cazy subcommand: src.databases.cazy.* (dump_data, filter_data,
multi_taxa, download).

v3 downloads one bulk CAZy database dump (a zipped tsv) rather than crawling per-family HTML
pages the way v2 did, so the HTML fixtures under tests/test_inputs/test_inputs_crawler/ are
for a retrieval method v3 no longer uses and are not exercised here.

Network access (urlopen, Bio.Entrez) is always mocked.
"""


import sqlite3

from argparse import Namespace
from pathlib import Path
from unittest.mock import MagicMock, patch
from zipfile import ZipFile

import pytest

from src.databases.cazy.download import download_cazy, get_cazy_db_dump
from src.databases.cazy.dump_data import dump_cazy_txt
from src.databases.cazy.filter_data import apply_class_and_family_filters, apply_tax_filters, drop_subfamilies
from src.databases.cazy.multi_taxa import keep_first_taxa, process_multi_taxa
from src.sql.interface.temp_tables import build_temp_table


def _args(**overrides):
    defaults = {"retries": 2, "timeout": 30}
    defaults.update(overrides)
    return Namespace(**defaults)


@pytest.fixture
def temp_table_db(db_path):
    build_temp_table(db_path)
    return db_path


def _insert_rows(db_path, rows):
    """rows: list of (family, kingdom, genus, species, protein_id, source) tuples."""
    conn = sqlite3.connect(db_path)
    conn.executemany(
        "INSERT INTO TEMP_TABLE (family, kingdom, genus, species, protein_id, source) "
        "VALUES (?, ?, ?, ?, ?, ?)",
        rows,
    )
    conn.commit()
    conn.close()


def _all_protein_ids(db_path):
    conn = sqlite3.connect(db_path)
    ids = {row[0] for row in conn.execute("SELECT protein_id FROM TEMP_TABLE")}
    conn.close()
    return ids


class TestDumpCazyTxt:
    def _make_cazy_zip(self, tmp_path, lines):
        txt_path = tmp_path / "cazy_data.txt"
        txt_path.write_text("\n".join(lines) + "\n")
        zip_path = tmp_path / "cazy_data.zip"
        with ZipFile(zip_path, "w") as zf:
            zf.write(txt_path, arcname="cazy_data.txt")
        return zip_path

    def test_parses_a_realistic_dump_line(self, temp_table_db, tmp_path):
        zip_path = self._make_cazy_zip(tmp_path, [
            "GH157\tBacteria\tBacteroides cellulosilyticus BFG-250\tUBD70155.1\tncbi",
        ])

        dump_cazy_txt(zip_path, temp_table_db)

        conn = sqlite3.connect(temp_table_db)
        row = conn.execute(
            "SELECT family, kingdom, genus, species, protein_id, source FROM TEMP_TABLE"
        ).fetchone()
        conn.close()

        # genus is split off into its own column, so species holds only what remains
        # (species epithet + strain), not the full binomial
        assert row == ("GH157", "Bacteria", "Bacteroides", "cellulosilyticus BFG-250", "UBD70155.1", "ncbi")

    def test_species_with_no_strain(self, temp_table_db, tmp_path):
        zip_path = self._make_cazy_zip(tmp_path, [
            "GT2\tBacteria\tEscherichia coli\tWP_000001.1\tncbi",
        ])

        dump_cazy_txt(zip_path, temp_table_db)

        conn = sqlite3.connect(temp_table_db)
        row = conn.execute("SELECT genus, species FROM TEMP_TABLE").fetchone()
        conn.close()
        assert row == ("Escherichia", "coli")

    def test_multiple_lines_are_all_dumped(self, temp_table_db, tmp_path):
        zip_path = self._make_cazy_zip(tmp_path, [
            "GH157\tBacteria\tBacteroides cellulosilyticus\tUBD70155.1\tncbi",
            "GT2\tArchaea\tMethanocatella smithii\tWP_000001.1\tncbi",
        ])

        dump_cazy_txt(zip_path, temp_table_db)

        assert _all_protein_ids(temp_table_db) == {"UBD70155.1", "WP_000001.1"}


class TestApplyTaxFilters:
    def _seed(self, db_path):
        _insert_rows(db_path, [
            ("GH1", "Bacteria", "Escherichia", "coli", "P1", "ncbi"),
            ("GH1", "Archaea", "Methanocatella", "smithii", "P2", "ncbi"),
            ("GH1", "Eukaryota", "Trichoderma", "reesei", "P3", "ncbi"),
        ])

    def test_kingdom_filter_keeps_only_matching_kingdoms(self, temp_table_db):
        self._seed(temp_table_db)

        apply_tax_filters({"Bacteria"}, set(), set(), set(), temp_table_db)

        assert _all_protein_ids(temp_table_db) == {"P1"}

    def test_genus_filter_keeps_only_matching_genera(self, temp_table_db):
        self._seed(temp_table_db)

        apply_tax_filters(set(), {"Trichoderma"}, set(), set(), temp_table_db)

        assert _all_protein_ids(temp_table_db) == {"P3"}

    def test_species_filter_is_a_prefix_match(self, temp_table_db):
        _insert_rows(temp_table_db, [
            ("GH1", "Bacteria", "Escherichia", "coli K12", "P1", "ncbi"),
            ("GH1", "Bacteria", "Escherichia", "albertii", "P2", "ncbi"),
        ])

        apply_tax_filters(set(), set(), {"coli"}, set(), temp_table_db)

        assert _all_protein_ids(temp_table_db) == {"P1"}

    def test_strain_filter_is_an_exact_match(self, temp_table_db):
        _insert_rows(temp_table_db, [
            ("GH1", "Bacteria", "Escherichia", "coli K12", "P1", "ncbi"),
            ("GH1", "Bacteria", "Escherichia", "coli K12 substrain", "P2", "ncbi"),
        ])

        apply_tax_filters(set(), set(), set(), {"coli K12"}, temp_table_db)

        assert _all_protein_ids(temp_table_db) == {"P1"}

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "BUG (found 2026-08-31, not yet fixed): apply_tax_filters combines multiple "
            "taxonomy filters with 'kingdom NOT IN (...) AND genus NOT IN (...)'. That only "
            "deletes a row when it fails BOTH filters at once; De Morgan's law means the "
            "correct combination is OR ('NOT (A AND B)' == 'NOT A OR NOT B'). Confirmed live: "
            "a row that satisfies only one of two combined filters is wrongly kept. This is "
            "reachable directly from the CLI any time two of --kingdoms/--genera/--species/"
            "--strains are given together, e.g. `download_cazy ... --kingdoms bacteria "
            "--genera Escherichia` silently returns organisms that are neither."
        ),
    )
    def test_kingdom_and_genus_filters_combine_with_and(self, temp_table_db):
        _insert_rows(temp_table_db, [
            ("GH1", "Bacteria", "Escherichia", "coli", "P1", "ncbi"),
            ("GH1", "Bacteria", "Bacillus", "subtilis", "P2", "ncbi"),
            ("GH1", "Archaea", "Escherichia", "coli", "P3", "ncbi"),  # nonsense combo but exercises AND
        ])

        apply_tax_filters({"Bacteria"}, {"Escherichia"}, set(), set(), temp_table_db)

        assert _all_protein_ids(temp_table_db) == {"P1"}

    def test_no_filters_raises_rather_than_being_called(self, temp_table_db):
        # apply_tax_filters builds "DELETE FROM TEMP_TABLE WHERE" with no condition when every
        # filter is empty, which is invalid SQL - not reachable via the CLI, since
        # src/scripts/scrape_cazy.py only calls this function when at least one taxonomy
        # filter is set (`if any([kingdom_filters, ...genera, species, strains]):`), but the
        # function itself has no such guard if called directly
        self._seed(temp_table_db)

        with pytest.raises(sqlite3.OperationalError):
            apply_tax_filters(set(), set(), set(), set(), temp_table_db)


class TestApplyClassAndFamilyFilters:
    def _seed(self, db_path):
        _insert_rows(db_path, [
            ("GH13", "Bacteria", "E", "coli", "P1", "ncbi"),
            ("GH13_1", "Bacteria", "E", "coli", "P2", "ncbi"),
            ("GH1", "Bacteria", "E", "coli", "P3", "ncbi"),
            ("CE1", "Bacteria", "E", "coli", "P4", "ncbi"),
        ])

    def test_class_filter_keeps_matching_families_and_their_subfamilies(self, temp_table_db):
        self._seed(temp_table_db)

        apply_class_and_family_filters(["GH"], set(), True, temp_table_db)

        assert _all_protein_ids(temp_table_db) == {"P1", "P2", "P3"}

    def test_family_filter_with_subfamilies_kept_includes_the_subfamily_rows(self, temp_table_db):
        self._seed(temp_table_db)

        apply_class_and_family_filters([], {"GH13"}, True, temp_table_db)

        assert _all_protein_ids(temp_table_db) == {"P1", "P2"}

    def test_family_filter_without_keeping_subfamilies_drops_the_subfamily_rows(self, temp_table_db):
        self._seed(temp_table_db)

        apply_class_and_family_filters([], {"GH13"}, False, temp_table_db)

        assert _all_protein_ids(temp_table_db) == {"P1"}

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "BUG (found 2026-08-31, not yet fixed): a subfamily-shaped --families entry given "
            "without --subfamilies logs a warning and `continue`s, but the loop already wrote "
            "an opening '(' before the continue and nothing closes it, so the query becomes "
            "'DELETE FROM TEMP_TABLE WHERE ()' - confirmed live to raise sqlite3.OperationalError "
            "('near \")\": syntax error'). Confirmed src/utilities/sanity_checks/scrape_cazy.py "
            "has no check rejecting this combination first, so it's directly reachable from the "
            "CLI: `download_cazy ... --families GH13_1` (without --subfamilies) crashes."
        ),
    )
    def test_explicit_subfamily_filter_without_the_keep_subfamilies_flag_has_no_effect_and_is_not_an_error(self, temp_table_db):
        self._seed(temp_table_db)

        apply_class_and_family_filters([], {"GH13_1"}, False, temp_table_db)

        assert _all_protein_ids(temp_table_db) == {"P1", "P2", "P3", "P4"}


class TestDropSubfamilies:
    def test_removes_only_rows_whose_family_has_an_underscore(self, temp_table_db):
        _insert_rows(temp_table_db, [
            ("GH13", "Bacteria", "E", "coli", "P1", "ncbi"),
            ("GH13_1", "Bacteria", "E", "coli", "P2", "ncbi"),
        ])

        drop_subfamilies(temp_table_db)

        assert _all_protein_ids(temp_table_db) == {"P1"}


class TestKeepFirstTaxa:
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "BUG (found 2026-08-31, not yet fixed): keep_first_taxa reads first_row[1] and "
            "first_row[2] expecting genus and species, but TEMP_TABLE's real column order is "
            "(record_id, family, kingdom, genus, species, protein_id, source) - indices 1 and 2 "
            "are actually family and kingdom. Confirmed live: since every row for one protein_id "
            "shares the same family and kingdom, 'genus != family_value OR species != "
            "kingdom_value' is true for every row including the first one, so this DELETEs ALL "
            "rows for the protein - silent data loss - instead of keeping one. Reachable via "
            "`download_cazy ... --skip_ncbi_tax` for any protein CAZy lists under more than one "
            "organism (the very case this function exists to resolve), and as the NCBI-lookup "
            "failure fallback inside use_latest_taxa even without --skip_ncbi_tax."
        ),
    )
    def test_deletes_rows_with_a_different_taxonomy_to_the_first_one_found(self, temp_table_db):
        _insert_rows(temp_table_db, [
            ("GH1", "Bacteria", "Escherichia", "coli", "SAME_ID", "ncbi"),
            ("GH1", "Bacteria", "Bacillus", "subtilis", "SAME_ID", "ncbi"),
        ])
        conn = sqlite3.connect(temp_table_db)

        keep_first_taxa("SAME_ID", conn)

        remaining = conn.execute("SELECT genus, species FROM TEMP_TABLE WHERE protein_id = 'SAME_ID'").fetchall()
        conn.close()
        assert remaining == [("Escherichia", "coli")]


class TestProcessMultiTaxa:
    @pytest.mark.xfail(
        strict=True,
        reason="Downstream consequence of the keep_first_taxa column-index bug documented on "
        "TestKeepFirstTaxa - --skip_ncbi_tax deletes the multi-taxa protein entirely instead "
        "of keeping its first-seen taxonomy.",
    )
    def test_skip_ncbi_tax_keeps_the_first_taxonomy_and_logs_the_protein_id(self, temp_table_db, tmp_path):
        _insert_rows(temp_table_db, [
            ("GH1", "Bacteria", "Escherichia", "coli", "MULTI_ID", "ncbi"),
            ("GH1", "Bacteria", "Bacillus", "subtilis", "MULTI_ID", "ncbi"),
            ("GH1", "Bacteria", "Escherichia", "coli", "SINGLE_ID", "ncbi"),
        ])
        multi_log = tmp_path / "multi_taxa.log"
        replaced_log = tmp_path / "replaced_taxa.log"

        process_multi_taxa(temp_table_db, multi_log, replaced_log, _args(skip_ncbi_tax=True))

        assert multi_log.read_text().strip() == "MULTI_ID"
        remaining = _all_protein_ids(temp_table_db)
        assert remaining == {"MULTI_ID", "SINGLE_ID"}
        conn = sqlite3.connect(temp_table_db)
        assert conn.execute(
            "SELECT COUNT(*) FROM TEMP_TABLE WHERE protein_id = 'MULTI_ID'"
        ).fetchone()[0] == 1

    def test_jgi_proteins_are_excluded_from_multi_taxa_handling(self, temp_table_db, tmp_path):
        # JGI assigns a bare integer id that cannot be tracked back to a specific JGI entry,
        # so ambiguous JGI taxonomy is left alone rather than collapsed
        _insert_rows(temp_table_db, [
            ("GH1", "Bacteria", "Escherichia", "coli", "999", "jgi"),
            ("GH1", "Bacteria", "Bacillus", "subtilis", "999", "jgi"),
        ])
        multi_log = tmp_path / "multi_taxa.log"
        replaced_log = tmp_path / "replaced_taxa.log"

        process_multi_taxa(temp_table_db, multi_log, replaced_log, _args(skip_ncbi_tax=True))

        # the log file is only opened once there is at least one multi-taxa id to write, so
        # with zero (every candidate excluded as JGI) it is never created at all
        assert not multi_log.exists()
        assert _all_protein_ids(temp_table_db) == {"999"}


class TestDownloadCazy:
    def _fake_response(self, chunks):
        response = MagicMock()
        response.info.return_value = {"Content-length": str(sum(len(c) for c in chunks))}
        response.read.side_effect = chunks + [b""]
        response.__enter__.return_value = response
        response.__exit__.return_value = False
        return response

    @patch("src.databases.cazy.download.urlopen")
    def test_writes_the_downloaded_bytes_to_the_output_path(self, mock_urlopen, tmp_path):
        mock_urlopen.return_value = self._fake_response([b"chunk1", b"chunk2"])
        out_path = tmp_path / "cazy_db.zip"

        err = download_cazy(out_path, _args(timeout=10), max_tries=1)

        assert err is None
        assert out_path.read_bytes() == b"chunk1chunk2"

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "BUG (found 2026-08-31, not yet fixed): download_file_decorator's wrapper() "
            "resets `err_message` to None at the top of each retry loop iteration, but never "
            "resets the separate `err` variable it actually checks ('if err is None: success = "
            "True'). Confirmed live: once one attempt fails, `err` keeps holding that "
            "exception forever, so every later attempt is reported as failed even when func() "
            "raises nothing that time. The file itself IS still written correctly on the "
            "successful attempt (func()'s side effects aren't affected), but the wrapper's "
            "return value stays wrong, which makes get_cazy_db_dump's own outer retry loop "
            "re-invoke download_cazy needlessly after any transient failure."
        ),
    )
    def test_retries_on_connection_failure_then_succeeds(self, mock_urlopen, mock_sleep, tmp_path):
        mock_urlopen.side_effect = [IOError("boom"), self._fake_response([b"data"])]
        out_path = tmp_path / "cazy_db.zip"

        err = download_cazy(out_path, _args(timeout=10), max_tries=2)

        assert err is None
        assert mock_urlopen.call_count == 2

    @patch("src.databases.cazy.download.time.sleep")
    @patch("src.databases.cazy.download.urlopen")
    def test_exhausting_retries_returns_the_error_rather_than_raising(self, mock_urlopen, mock_sleep, tmp_path):
        mock_urlopen.side_effect = IOError("boom")
        out_path = tmp_path / "cazy_db.zip"

        err = download_cazy(out_path, _args(timeout=10), max_tries=2)

        assert err is not None
        assert mock_urlopen.call_count == 2


class TestGetCazyDbDump:
    def test_a_predownloaded_file_is_used_without_attempting_a_download(self, tmp_path):
        predownloaded = tmp_path / "already_have_this.zip"

        with patch("src.databases.cazy.download.download_cazy") as mock_download:
            result = get_cazy_db_dump(tmp_path, "2026-08-31_12-00-00", _args(cazy_data=predownloaded))
            mock_download.assert_not_called()

        assert result == predownloaded

    def test_no_predownloaded_file_triggers_a_download_to_the_cache_dir(self, tmp_path):
        with patch("src.databases.cazy.download.download_cazy") as mock_download:
            mock_download.return_value = None  # None == success, per download_file_decorator

            result = get_cazy_db_dump(tmp_path, "2026-08-31_12-00-00", _args(cazy_data=None))

        assert result == tmp_path / "cazy_db_2026-08-31_12-00-00.zip"
        mock_download.assert_called_once()
