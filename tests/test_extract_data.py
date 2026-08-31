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
"""Unit tests for the extract_data additions: `--include pfam` in
src.sql.interface.get_data.get_extract_data, and the tsv/jsonl writers (plus the refactored
csv/json writers they share code with) in src.export.tabular.

Self-contained: builds its own v3-schema db per test via tmp_path, independent of the
(v2-only, unrelated) fixtures in tests/conftest.py.
"""


import json
import sqlite3

import pytest

from src.export.tabular import write_csv, write_json, write_jsonl, write_tsv
from src.sql import sql_orm
from src.sql.interface.get_data.get_extract_data import iter_protein_records


@pytest.fixture
def db_path(tmp_path):
    path = tmp_path / "test.db"
    sql_orm.get_db_connection(path, False, new=True)

    conn = sqlite3.connect(path)
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

    return path


class TestPfamInclude:
    def test_pfam_accessions_are_included_when_requested(self, db_path):
        batches = list(iter_protein_records(["TEST_ACC"], db_path, {"pfam"}))

        assert len(batches) == 1
        record = batches[0]["TEST_ACC"]
        assert record["pfam"] == {"PF00051", "PF00089"}

    def test_pfam_is_absent_from_the_record_when_not_requested(self, db_path):
        batches = list(iter_protein_records(["TEST_ACC"], db_path, {"class"}))

        assert "pfam" not in batches[0]["TEST_ACC"]

    def test_protein_with_no_pfam_matches_gets_an_empty_set_not_a_missing_key(self, db_path):
        conn = sqlite3.connect(db_path)
        conn.execute("INSERT INTO Proteins (protein_accession, source) VALUES ('NO_PFAM', 'ncbi')")
        conn.commit()
        conn.close()

        batches = list(iter_protein_records(["NO_PFAM"], db_path, {"pfam"}))

        assert batches[0]["NO_PFAM"]["pfam"] == set()


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
