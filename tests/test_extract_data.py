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
"""Unit tests for the extract_data subcommand: src.sql.interface.get_data.get_extract_data
(the `--include` assembly logic) and src.export.tabular (the csv/tsv/json/jsonl writers).

Uses the shared `db_path` fixture from tests/conftest.py for a fresh, empty v3-schema db.
"""


import json
import sqlite3

import pytest

from src.export.tabular import write_csv, write_json, write_jsonl, write_tsv
from src.sql.interface.get_data.get_extract_data import iter_protein_records, iter_sequences


@pytest.fixture
def pfam_extract_db(db_path):
    """db_path, populated with one protein carrying two Pfam matches."""
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("INSERT INTO Uniprots (uniprot_accession) VALUES ('P00734')")
    uniprot_id = cur.lastrowid
    cur.execute(
        "INSERT INTO Proteins (protein_accession, source, uniprot_id) VALUES ('TEST_ACC', 'ncbi', ?)",
        (uniprot_id,),
    )
    protein_id = cur.lastrowid

    cur.execute(
        "INSERT INTO Pfams (accession, name, annotation_type, release) VALUES "
        "('PF00051', 'Kringle domain', 'domain', '109.0'), "
        "('PF00089', 'Trypsin', 'domain', '109.0')"
    )
    pfam_ids = [row[0] for row in cur.execute("SELECT pfam_id FROM Pfams ORDER BY pfam_id")]
    cur.executemany(
        "INSERT INTO Proteins_Pfams (protein_id, pfam_id, interpro_accession, match_start, match_end) "
        "VALUES (?, ?, 'IPRXXXXXX', 1, 10)",
        [(protein_id, pfam_id) for pfam_id in pfam_ids],
    )
    conn.commit()
    conn.close()

    return db_path


class TestPfamInclude:
    def test_pfam_accessions_are_included_when_requested(self, pfam_extract_db):
        batches = list(iter_protein_records(["TEST_ACC"], pfam_extract_db, {"pfam"}))

        assert len(batches) == 1
        record = batches[0]["TEST_ACC"]
        assert record["pfam"] == {"PF00051", "PF00089"}

    def test_pfam_is_absent_from_the_record_when_not_requested(self, pfam_extract_db):
        batches = list(iter_protein_records(["TEST_ACC"], pfam_extract_db, {"class"}))

        assert "pfam" not in batches[0]["TEST_ACC"]

    def test_protein_with_no_pfam_matches_gets_an_empty_set_not_a_missing_key(self, db_path):
        conn = sqlite3.connect(db_path)
        conn.execute("INSERT INTO Proteins (protein_accession, source) VALUES ('NO_PFAM', 'ncbi')")
        conn.commit()
        conn.close()

        batches = list(iter_protein_records(["NO_PFAM"], db_path, {"pfam"}))

        assert batches[0]["NO_PFAM"]["pfam"] == set()


@pytest.fixture
def full_extract_db(db_path):
    """db_path, populated with one protein carrying class/family/tax/ec/pdb/uniprot data."""
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute("INSERT INTO Kingdoms (kingdom) VALUES ('Bacteria')")
    kingdom_id = cur.lastrowid
    cur.execute(
        "INSERT INTO Taxs (genus, species, kingdom_id) VALUES ('Escherichia', 'coli', ?)",
        (kingdom_id,),
    )
    taxonomy_id = cur.lastrowid

    cur.execute(
        "INSERT INTO Uniprots (uniprot_accession, name, sequence) VALUES ('P00734', 'Prothrombin', 'MAAA')"
    )
    uniprot_id = cur.lastrowid

    cur.execute(
        "INSERT INTO Proteins (protein_accession, source, sequence, taxonomy_id, uniprot_id) "
        "VALUES ('TEST_ACC', 'ncbi', 'MGGG', ?, ?)",
        (taxonomy_id, uniprot_id),
    )
    protein_id = cur.lastrowid

    cur.execute("INSERT INTO CazyFamilies (family, subfamily) VALUES ('GH13', '_1')")
    family_id = cur.lastrowid
    cur.execute(
        "INSERT INTO Proteins_CazyFamilies (protein_id, family_id) VALUES (?, ?)",
        (protein_id, family_id),
    )

    cur.execute("INSERT INTO Ecs (ec_number) VALUES ('3.2.1.1')")
    ec_id = cur.lastrowid
    cur.execute("INSERT INTO Proteins_Ecs (protein_id, ec_id) VALUES (?, ?)", (protein_id, ec_id))

    cur.execute("INSERT INTO Pdbs (pdb_accession, method, resolution) VALUES ('1ABC', 'X-ray', 1.5)")
    pdb_id = cur.lastrowid
    cur.execute("INSERT INTO Proteins_Pdbs (protein_id, pdb_id) VALUES (?, ?)", (protein_id, pdb_id))

    conn.commit()
    conn.close()
    return db_path


class TestOtherIncludeFields:
    def test_class_family_subfamily(self, full_extract_db):
        record = list(iter_protein_records(["TEST_ACC"], full_extract_db, {"class", "family", "subfamily"}))[0]["TEST_ACC"]

        assert record["class"] == {"GH"}
        assert record["family"] == {"GH13"}
        assert record["subfamily"] == {"_1"}

    def test_kingdom_genus_organism(self, full_extract_db):
        record = list(iter_protein_records(["TEST_ACC"], full_extract_db, {"kingdom", "genus", "organism"}))[0]["TEST_ACC"]

        assert record["kingdom"] == "Bacteria"
        assert record["genus"] == "Escherichia"
        assert record["organism"] == "Escherichia coli"

    def test_ec(self, full_extract_db):
        record = list(iter_protein_records(["TEST_ACC"], full_extract_db, {"ec"}))[0]["TEST_ACC"]
        assert record["ec"] == {"3.2.1.1"}

    def test_pdb(self, full_extract_db):
        record = list(iter_protein_records(["TEST_ACC"], full_extract_db, {"pdb"}))[0]["TEST_ACC"]
        assert record["pdb"] == {"1ABC"}

    def test_uniprot_acc_and_name(self, full_extract_db):
        record = list(iter_protein_records(["TEST_ACC"], full_extract_db, {"uniprot_acc", "uniprot_name"}))[0]["TEST_ACC"]
        assert record["uniprot_acc"] == "P00734"
        assert record["uniprot_name"] == "Prothrombin"

    def test_genbank_and_uniprot_sequences(self, full_extract_db):
        record = list(iter_protein_records(["TEST_ACC"], full_extract_db, {"genbank_seq", "uniprot_seq"}))[0]["TEST_ACC"]
        assert record["genbank_seq"] == "MGGG"
        assert record["uniprot_seq"] == "MAAA"

    def test_only_requested_fields_are_populated(self, full_extract_db):
        record = list(iter_protein_records(["TEST_ACC"], full_extract_db, {"ec"}))[0]["TEST_ACC"]
        assert set(record) == {"ec"}


class TestIterSequences:
    def test_genbank_source_yields_ncbi_sequence(self, full_extract_db):
        results = list(iter_sequences(["TEST_ACC"], full_extract_db, {"genbank"}))
        assert results == [("TEST_ACC", "GenBank", "MGGG")]

    def test_uniprot_source_yields_uniprot_sequence_keyed_by_uniprot_accession(self, full_extract_db):
        results = list(iter_sequences(["TEST_ACC"], full_extract_db, {"uniprot"}))
        assert results == [("P00734", "UniProt", "MAAA")]

    def test_both_sources_together(self, full_extract_db):
        results = list(iter_sequences(["TEST_ACC"], full_extract_db, {"genbank", "uniprot"}))
        assert set(results) == {("TEST_ACC", "GenBank", "MGGG"), ("P00734", "UniProt", "MAAA")}

    def test_protein_with_no_sequence_yields_nothing(self, db_path):
        conn = sqlite3.connect(db_path)
        conn.execute("INSERT INTO Proteins (protein_accession, source) VALUES ('NO_SEQ', 'ncbi')")
        conn.commit()
        conn.close()

        assert list(iter_sequences(["NO_SEQ"], db_path, {"genbank"})) == []


# --- writer tests: pure functions, no db needed ---

RECORD_BATCHES = [
    {
        "B_ACC": {"pfam": {"PF00089"}, "ec": set()},
        "A_ACC": {"pfam": {"PF00051", "PF00089"}, "ec": {"1.1.1.1"}},
    },
]


class TestWriteCsvAndTsv:
    def test_csv_uses_commas_and_sorts_rows_by_accession(self, tmp_path):
        out = tmp_path / "out.csv"

        rows = write_csv(RECORD_BATCHES, out, {"ec", "pfam"})

        assert rows == 2
        lines = out.read_text().splitlines()
        assert lines[0] == "protein_accession,ec,pfam"
        assert lines[1] == 'A_ACC,1.1.1.1,"PF00051,PF00089"'
        assert lines[2] == "B_ACC,,PF00089"

    def test_tsv_uses_tabs(self, tmp_path):
        out = tmp_path / "out.tsv"

        rows = write_tsv(RECORD_BATCHES, out, {"ec", "pfam"})

        assert rows == 2
        lines = out.read_text().splitlines()
        assert lines[0] == "protein_accession\tec\tpfam"
        assert lines[1] == "A_ACC\t1.1.1.1\tPF00051,PF00089"


class TestWriteJsonAndJsonl:
    def test_json_is_one_document_keyed_by_accession(self, tmp_path):
        out = tmp_path / "out.json"

        records = write_json(RECORD_BATCHES, out, {"pfam"})

        assert records == 2
        doc = json.loads(out.read_text())
        assert set(doc) == {"A_ACC", "B_ACC"}
        assert sorted(doc["A_ACC"]["pfam"]) == ["PF00051", "PF00089"]

    def test_jsonl_is_one_self_contained_object_per_line(self, tmp_path):
        out = tmp_path / "out.jsonl"

        records = write_jsonl(RECORD_BATCHES, out, {"pfam"})

        assert records == 2
        lines = out.read_text().splitlines()
        assert len(lines) == 2

        parsed = [json.loads(line) for line in lines]
        by_accession = {p["protein_accession"]: p for p in parsed}
        assert sorted(by_accession["A_ACC"]["pfam"]) == ["PF00051", "PF00089"]
        assert by_accession["B_ACC"]["pfam"] == ["PF00089"]

    def test_jsonl_lines_are_independently_parseable(self, tmp_path):
        # the point of jsonl over json: each line stands alone, not just the whole file
        out = tmp_path / "out.jsonl"
        write_jsonl(RECORD_BATCHES, out, {"pfam"})

        for line in out.read_text().splitlines():
            json.loads(line)  # raises if not valid standalone JSON
