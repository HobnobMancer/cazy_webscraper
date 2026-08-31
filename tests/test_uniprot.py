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
"""Unit tests for the get_uniprot_data subcommand: src.databases.uniprot (the bioservices
UniProt ID-mapping client) and src.sql.interface.add_data.add_uniprot_data (the db layer).

Network access (bioservices) is always mocked.
"""


import sqlite3

from unittest.mock import MagicMock

import pytest

from src.databases.uniprot import UniProtKB, UniProtRecord
from src.sql.interface.add_data.add_uniprot_data import (
    add_ec_numbers, add_go_terms, add_pdbs, add_uniprot_records,
    merge_temp_ec_relationships, merge_temp_go_relationships, merge_temp_pdb_relationships,
)


# Shape matches a real UniProt ID-mapping /idmapping/results response, as parsed by
# UniProtKB.parse_mappings.
def _mapping(ncbi_acc="WP_012345.1", uniprot_id="P00734", swissprot=True):
    return {
        "from": ncbi_acc,
        "to": {
            "uniProtkbId": uniprot_id,
            "extraAttributes": {"uniParcId": "UPI000000001"},
            "proteinDescription": {
                "recommendedName": {
                    "fullName": {"value": "Prothrombin"},
                    "ecNumbers": [{"value": "3.4.21.5"}],
                },
            },
            "comments": [
                {"commentType": "CATALYTIC ACTIVITY", "reaction": {"ecNumber": "3.4.21.5"}},
            ],
            "entryType": "UniProtKB reviewed (Swiss-Prot)" if swissprot else "UniProtKB unreviewed (TrEMBL)",
            "sequence": {"value": "MAAAGGG", "md5": "abcd1234", "molWeight": 1234.5},
            "entryAudit": {"lastSequenceUpdateDate": "2020-01-01"},
            "uniProtKBCrossReferences": [
                {
                    "database": "PDB",
                    "id": "1ABC",
                    "properties": [
                        {"key": "Method", "value": "X-ray"},
                        {"key": "Resolution", "value": "1.50 A"},
                    ],
                },
                {
                    "database": "GO",
                    "id": "GO:0004252",
                    "properties": [{"key": "GoTerm", "value": "F:serine-type endopeptidase activity"}],
                },
            ],
        },
    }


class TestParseMappings:
    def test_extracts_the_core_fields(self):
        records = UniProtKB(retries=1).parse_mappings([_mapping()])

        assert len(records) == 1
        record = records[0]
        assert record.ncbi_acc == "WP_012345.1"
        assert record.uniprot_id == "P00734"
        assert record.uniparc_id == "UPI000000001"
        assert record.name == "Prothrombin"
        assert record.swissprot is True
        assert record.sequence == "MAAAGGG"
        assert record.md5 == "abcd1234"
        assert record.mol_weight == 1234.5
        assert record.sequence_date == "2020-01-01"

    def test_ec_numbers_come_from_both_recommended_name_and_catalytic_activity_comments(self):
        record = UniProtKB(retries=1).parse_mappings([_mapping()])[0]
        assert record.ec_nums == {"3.4.21.5"}

    def test_pdb_cross_references_include_method_and_parsed_resolution(self):
        record = UniProtKB(retries=1).parse_mappings([_mapping()])[0]
        assert record.pdbs == {("1ABC", "X-ray", 1.5)}

    def test_go_cross_references_are_keyed_by_go_id(self):
        record = UniProtKB(retries=1).parse_mappings([_mapping()])[0]
        assert record.go_terms == {"GO:0004252": "F:serine-type endopeptidase activity"}

    def test_swissprot_only_filters_out_trembl_entries(self):
        records = UniProtKB(retries=1).parse_mappings(
            [_mapping(swissprot=False)], swissprot_only=True,
        )
        assert records == []

    def test_get_flags_control_what_is_parsed(self):
        record = UniProtKB(retries=1).parse_mappings(
            [_mapping()], get_ecs=False, get_pdbs=False, get_gos=False,
        )[0]
        assert record.ec_nums == set()
        assert record.pdbs == set()
        assert record.go_terms == {}

    def test_malformed_resolution_is_dropped_not_raised(self):
        mapping = _mapping()
        mapping["to"]["uniProtKBCrossReferences"][0]["properties"][1]["value"] = "not-a-number"
        record = UniProtKB(retries=1).parse_mappings([mapping])[0]
        method, resolution = next(iter(record.pdbs))[1:]
        assert method == "X-ray"
        assert resolution is None


class TestMapBatch:
    def test_returns_results_and_failed_ids_from_a_single_successful_call(self):
        client = UniProtKB(retries=2)
        client.service = MagicMock()
        client.service.mapping.return_value = {"results": [_mapping()], "failedIds": ["BAD.1"]}

        results, failed = client._map_batch(["WP_012345.1", "BAD.1"], "EMBL-GenBank-DDBJ_CDS")

        assert results == [_mapping()]
        assert failed == ["BAD.1"]

    def test_retries_on_a_falsy_response_then_succeeds(self, monkeypatch):
        monkeypatch.setattr("src.databases.uniprot.time.sleep", lambda *_: None)
        client = UniProtKB(retries=2)
        client.service = MagicMock()
        client.service.mapping.side_effect = [None, {"results": [_mapping()], "failedIds": []}]

        results, failed = client._map_batch(["WP_012345.1"], "EMBL-GenBank-DDBJ_CDS")

        assert len(results) == 1
        assert client.service.mapping.call_count == 2

    def test_exhausting_retries_returns_the_whole_batch_as_failed(self, monkeypatch):
        monkeypatch.setattr("src.databases.uniprot.time.sleep", lambda *_: None)
        client = UniProtKB(retries=1)
        client.service = MagicMock()
        client.service.mapping.return_value = None

        results, failed = client._map_batch(["WP_012345.1", "WP_999999.1"], "EMBL-GenBank-DDBJ_CDS")

        assert results == []
        assert failed == ["WP_012345.1", "WP_999999.1"]


class TestGetUniprotRecords:
    def test_ids_that_fail_as_genbank_are_retried_as_refseq(self, monkeypatch):
        # UniProt indexes RefSeq accessions (NP_/XP_/WP_) separately from classic
        # GenBank/EMBL/DDBJ CDS ids - a WP_ accession is expected to miss on the first pass
        client = UniProtKB(retries=1)
        calls = []

        def fake_map_batch(batch, fr):
            calls.append(fr)
            if fr == "EMBL-GenBank-DDBJ_CDS":
                return [], batch
            return [_mapping(ncbi_acc="WP_012345.1")], []

        monkeypatch.setattr(client, "_map_batch", fake_map_batch)

        records, failed_ids = client.get_uniprot_records(["WP_012345.1"], False, True, True, True)

        assert calls == ["EMBL-GenBank-DDBJ_CDS", "RefSeq_Protein"]
        assert len(records) == 1
        assert failed_ids == []


class TestAddUniprotRecords:
    def _record(self, ncbi_acc="TEST_ACC", uniprot_id="P00734", sequence="MAAA"):
        return UniProtRecord(
            ncbi_acc, uniprot_id, "UPI0001", "Prothrombin", True, sequence, "md5hash", 123.4,
            "2020-01-01", set(), set(), {},
        )

    def _add_protein(self, db_path, accession):
        conn = sqlite3.connect(db_path)
        conn.execute("INSERT INTO Proteins (protein_accession, source) VALUES (?, 'ncbi')", (accession,))
        conn.commit()
        conn.close()

    def test_new_uniprot_record_is_added_and_linked_to_the_protein(self, db_path):
        self._add_protein(db_path, "TEST_ACC")

        new_count, updated_count = add_uniprot_records([self._record()], db_path, update_seq=False)

        assert new_count == 1
        assert updated_count == 0
        conn = sqlite3.connect(db_path)
        row = conn.execute(
            "SELECT U.uniprot_accession FROM Proteins P JOIN Uniprots U ON P.uniprot_id = U.uniprot_id "
            "WHERE P.protein_accession = 'TEST_ACC'"
        ).fetchone()
        assert row == ("P00734",)
        conn.close()

    def test_existing_record_without_a_sequence_is_backfilled_even_without_update_flag(self, db_path):
        self._add_protein(db_path, "TEST_ACC")
        conn = sqlite3.connect(db_path)
        conn.execute("INSERT INTO Uniprots (uniprot_accession, sequence) VALUES ('P00734', NULL)")
        conn.commit()
        conn.close()

        new_count, updated_count = add_uniprot_records([self._record()], db_path, update_seq=False)

        assert new_count == 0
        assert updated_count == 1

    def test_existing_sequence_is_untouched_without_update_flag(self, db_path):
        self._add_protein(db_path, "TEST_ACC")
        conn = sqlite3.connect(db_path)
        conn.execute("INSERT INTO Uniprots (uniprot_accession, sequence) VALUES ('P00734', 'OLDSEQ')")
        conn.commit()
        conn.close()

        new_count, updated_count = add_uniprot_records([self._record()], db_path, update_seq=False)

        assert new_count == 0
        assert updated_count == 0


class TestAddEcNumbersPdbsGoTerms:
    def _linked_record(self, db_path):
        conn = sqlite3.connect(db_path)
        conn.execute("INSERT INTO Proteins (protein_accession, source) VALUES ('TEST_ACC', 'ncbi')")
        protein_id = conn.execute("SELECT protein_id FROM Proteins").fetchone()[0]
        conn.commit()
        conn.close()

        record = UniProtRecord(
            "TEST_ACC", "P00734", "UPI0001", "Prothrombin", True, "MAAA", "md5", 123.4,
            "2020-01-01", {"3.4.21.5"}, {("1ABC", "X-ray", 1.5)}, {"GO:0004252": "peptidase activity"},
        )
        record.protein_id = protein_id
        return record

    def test_add_ec_numbers_stages_the_relationship_and_merge_commits_it(self, db_path):
        record = self._linked_record(db_path)
        conn = sqlite3.connect(db_path)

        new_ecs = add_ec_numbers([record], conn)
        assert new_ecs == 1
        conn.close()

        merge_temp_ec_relationships(db_path)

        conn = sqlite3.connect(db_path)
        row = conn.execute(
            "SELECT E.ec_number FROM Proteins_Ecs PE JOIN Ecs E ON PE.ec_id = E.ec_id"
        ).fetchone()
        assert row == ("3.4.21.5",)
        conn.close()

    def test_add_pdbs_stages_the_relationship_and_merge_commits_it(self, db_path):
        record = self._linked_record(db_path)
        conn = sqlite3.connect(db_path)

        new_pdbs = add_pdbs([record], conn)
        assert new_pdbs == 1
        conn.close()

        merge_temp_pdb_relationships(db_path)

        conn = sqlite3.connect(db_path)
        row = conn.execute("SELECT pdb_accession, method, resolution FROM Pdbs").fetchone()
        assert row == ("1ABC", "X-ray", 1.5)
        conn.close()

    def test_add_go_terms_stages_the_relationship_and_merge_commits_it(self, db_path):
        record = self._linked_record(db_path)
        conn = sqlite3.connect(db_path)

        new_gos = add_go_terms([record], conn)
        assert new_gos == 1
        conn.close()

        merge_temp_go_relationships(db_path)

        conn = sqlite3.connect(db_path)
        row = conn.execute("SELECT goterm_id, description FROM GoTerms").fetchone()
        assert row == ("GO:0004252", "peptidase activity")
        conn.close()
