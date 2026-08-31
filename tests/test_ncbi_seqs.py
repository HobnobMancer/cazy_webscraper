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
"""Unit tests for the get_ncbi_seqs subcommand: src.databases.ncbi.sequences (retrieval) and
src.sql.interface.add_data.add_ncbi_seqs (the db layer).

Network access (Bio.Entrez) is always mocked.
"""


import io
import sqlite3

from argparse import Namespace
from unittest.mock import patch

import pytest
from Bio.Seq import Seq

from src.databases.ncbi.sequences import fetch_seqs_from_entrez
from src.sql.interface.add_data.add_ncbi_seqs import update_ncbi_seqs


def _args(**overrides):
    defaults = {"retries": 2}
    defaults.update(overrides)
    return Namespace(**defaults)


class TestFetchSeqsFromEntrez:
    @patch("src.databases.ncbi.sequences.Entrez.efetch")
    def test_parses_fasta_records_keyed_by_accession(self, mock_efetch):
        mock_efetch.return_value = io.StringIO(">WP_012345.1 some protein\nMKVLAA\n")

        seqs, invalid_err, connection_err = fetch_seqs_from_entrez(
            "env", "1", ["WP_012345.1"], _args(),
        )

        assert connection_err is False
        assert invalid_err is False
        assert seqs == {"WP_012345.1": Seq("MKVLAA")}

    @patch("src.databases.ncbi.sequences.Entrez.efetch")
    def test_a_record_whose_accession_is_not_in_the_requested_batch_is_dropped(self, mock_efetch):
        # NCBI can return records for accessions that were not asked for (e.g. a merged
        # record); these must not be attributed to some other accession in the batch
        mock_efetch.return_value = io.StringIO(">UNEXPECTED_ACC.1\nMKVLAA\n")

        seqs, invalid_err, connection_err = fetch_seqs_from_entrez(
            "env", "1", ["WP_012345.1"], _args(),
        )

        assert seqs == {}

    @patch("src.databases.ncbi.sequences.Entrez.efetch")
    def test_incomplete_read_is_a_connection_error(self, mock_efetch):
        from http.client import IncompleteRead
        mock_efetch.side_effect = IncompleteRead(b"")

        seqs, invalid_err, connection_err = fetch_seqs_from_entrez(
            "env", "1", ["WP_012345.1"], _args(),
        )

        assert connection_err is True
        assert seqs == {}


class TestUpdateNcbiSeqs:
    def _add_protein(self, db_path, accession, sequence=None):
        conn = sqlite3.connect(db_path)
        conn.execute(
            "INSERT INTO Proteins (protein_accession, source, sequence) VALUES (?, 'ncbi', ?)",
            (accession, sequence),
        )
        conn.commit()
        conn.close()

    def test_sequence_is_added_when_the_protein_has_none_yet(self, db_path):
        self._add_protein(db_path, "WP_012345.1")

        stats = update_ncbi_seqs({"WP_012345.1": Seq("MKVLAA")}, db_path, "2026-08-31_12-00-00", update=False)

        assert stats == {"added": 1, "updated": 0}
        conn = sqlite3.connect(db_path)
        assert conn.execute("SELECT sequence FROM Proteins WHERE protein_accession = 'WP_012345.1'").fetchone()[0] == "MKVLAA"
        conn.close()

    def test_existing_sequence_is_untouched_without_update_flag(self, db_path):
        self._add_protein(db_path, "WP_012345.1", sequence="OLDSEQ")

        stats = update_ncbi_seqs({"WP_012345.1": Seq("NEWSEQ")}, db_path, "2026-08-31_12-00-00", update=False)

        assert stats == {"added": 0, "updated": 0}
        conn = sqlite3.connect(db_path)
        assert conn.execute("SELECT sequence FROM Proteins WHERE protein_accession = 'WP_012345.1'").fetchone()[0] == "OLDSEQ"
        conn.close()

    def test_existing_sequence_is_overwritten_with_update_flag(self, db_path):
        self._add_protein(db_path, "WP_012345.1", sequence="OLDSEQ")

        stats = update_ncbi_seqs({"WP_012345.1": Seq("NEWSEQ")}, db_path, "2026-08-31_12-00-00", update=True)

        assert stats == {"added": 0, "updated": 1}
        conn = sqlite3.connect(db_path)
        assert conn.execute("SELECT sequence FROM Proteins WHERE protein_accession = 'WP_012345.1'").fetchone()[0] == "NEWSEQ"
        conn.close()

    def test_accession_not_in_the_local_db_is_skipped_without_error(self, db_path):
        stats = update_ncbi_seqs({"NOT_IN_DB.1": Seq("MKVLAA")}, db_path, "2026-08-31_12-00-00", update=False)

        assert stats == {"added": 0, "updated": 0}
