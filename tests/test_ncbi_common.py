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
"""Unit tests for src.databases.ncbi (the Entrez helpers shared by get_ncbi_seqs, get_ncbi_taxs
and get_ncbi_genomes: call_entrez, post_acc_to_entrez, validate_batch, get_protein_accession).

Network access (Bio.Entrez) is always mocked.
"""


from argparse import Namespace
from unittest.mock import MagicMock, patch

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.databases import ncbi


def _args(**overrides):
    defaults = {"retries": 2}
    defaults.update(overrides)
    return Namespace(**defaults)


class TestGetProteinAccession:
    def test_refseq_style_accession(self):
        record = SeqRecord(Seq("MKV"), id="WP_012345678.1")
        assert ncbi.get_protein_accession(record) == "WP_012345678.1"

    def test_genbank_cds_style_accession(self):
        record = SeqRecord(Seq("MKV"), id="AAA12345.1")
        assert ncbi.get_protein_accession(record) == "AAA12345.1"

    def test_pipe_delimited_id_with_sp_prefix(self):
        # e.g. "sp|P00734|THRB_HUMAN"
        record = SeqRecord(Seq("MKV"), id="sp|P00734|THRB_HUMAN")
        assert ncbi.get_protein_accession(record) == "P00734"

    def test_falls_back_to_a_known_accession_found_inside_the_id(self):
        record = SeqRecord(Seq("MKV"), id="lcl|weird_wrapper_WP_099999999.1_suffix")
        result = ncbi.get_protein_accession(record, acc_to_retrieve=["WP_099999999.1"])
        assert result == "WP_099999999.1"

    def test_returns_none_when_nothing_matches(self):
        record = SeqRecord(Seq("MKV"), id="###")
        assert ncbi.get_protein_accession(record) is None


class TestCallEntrez:
    @patch("src.databases.ncbi.entrez_retry")
    @patch("src.databases.ncbi.Entrez.read")
    def test_success_returns_result_with_no_errors(self, mock_read, mock_retry):
        mock_read.return_value = {"WebEnv": "env", "QueryKey": "1"}

        result, invalid_err, connection_err = ncbi.call_entrez(MagicMock(), retries=2, db="Protein", id="X")

        assert result == {"WebEnv": "env", "QueryKey": "1"}
        assert invalid_err is False
        assert connection_err is False

    @patch("src.databases.ncbi.entrez_retry")
    @patch("src.databases.ncbi.Entrez.read")
    def test_invalid_id_runtime_error_is_flagged_invalid_not_connection(self, mock_read, mock_retry):
        mock_read.side_effect = RuntimeError(
            "Some IDs have invalid value and were omitted. Only 998 out of 999"
        )

        result, invalid_err, connection_err = ncbi.call_entrez(MagicMock(), retries=2, db="Protein", id="X")

        assert invalid_err is True
        assert connection_err is False

    @patch("src.databases.ncbi.entrez_retry")
    @patch("src.databases.ncbi.Entrez.read")
    def test_incomplete_read_is_flagged_as_a_connection_error(self, mock_read, mock_retry):
        from http.client import IncompleteRead
        mock_read.side_effect = IncompleteRead(b"")

        result, invalid_err, connection_err = ncbi.call_entrez(MagicMock(), retries=2, db="Protein", id="X")

        assert connection_err is True
        assert invalid_err is False


class TestPostAccToEntrez:
    @patch("src.databases.ncbi.call_entrez")
    def test_extracts_webenv_and_query_key_on_success(self, mock_call_entrez):
        mock_call_entrez.return_value = ({"WebEnv": "env123", "QueryKey": "1"}, False, False)

        webenv, query_key, invalid_err, connection_err = ncbi.post_acc_to_entrez(["ACC1"], _args())

        assert webenv == "env123"
        assert query_key == "1"
        assert invalid_err is False
        assert connection_err is False

    @patch("src.databases.ncbi.call_entrez")
    def test_propagates_connection_error(self, mock_call_entrez):
        mock_call_entrez.return_value = (None, False, True)

        webenv, query_key, invalid_err, connection_err = ncbi.post_acc_to_entrez(["ACC1"], _args())

        assert connection_err is True
        assert webenv is None


class TestValidateBatch:
    @patch("src.databases.ncbi.post_acc_to_entrez")
    def test_whole_batch_valid_is_returned_unchanged(self, mock_post):
        mock_post.return_value = ("env", "1", False, False)

        valid, errs = ncbi.validate_batch(["A", "B", "C"], [], _args())

        assert valid == ["A", "B", "C"]
        assert errs == []

    @patch("src.databases.ncbi.post_acc_to_entrez")
    def test_connection_error_caches_the_batch_and_returns_nothing_valid(self, mock_post):
        mock_post.return_value = (None, None, False, True)

        valid, errs = ncbi.validate_batch(["A", "B"], [], _args())

        assert valid == []
        assert errs == [["A", "B"]]

    @patch("src.databases.ncbi.post_acc_to_entrez")
    def test_single_invalid_id_is_dropped_not_retried(self, mock_post):
        mock_post.return_value = (None, None, True, False)

        valid, errs = ncbi.validate_batch(["BAD_ACC"], [], _args())

        assert valid == []

    @patch("src.databases.ncbi.post_acc_to_entrez")
    def test_multi_id_batch_with_an_invalid_id_is_split_and_retried(self, mock_post):
        # first call (whole batch) reports invalid; subsequent calls (the two halves) succeed -
        # this is how a single bad accession is isolated out of an otherwise-valid batch
        mock_post.side_effect = [
            (None, None, True, False),
            ("env", "1", False, False),
            ("env", "1", False, False),
        ]

        valid, errs = ncbi.validate_batch(["A", "B"], [], _args())

        assert valid == ["A", "B"]
