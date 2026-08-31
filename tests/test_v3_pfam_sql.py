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
"""Unit tests for the SQL layer added for Pfam support: src.sql.sql_orm's ProteinPfam model,
src.sql.interface.get_data/add_data.*pfam_data, and the UniProt-accession selector added to
src.sql.interface.get_data.get_proteins.

Self-contained: builds its own v3-schema db per test via tmp_path, independent of the
(v2-only, unrelated) fixtures in tests/conftest.py.
"""


import sqlite3

import pytest

from src.databases.interpro import PfamMatchRecord
from src.sql import sql_orm
from src.sql.interface.add_data.add_pfam_data import add_pfam_matches, merge_temp_pfam_relationships
from src.sql.interface.get_data.get_pfam_data import get_db_pfams
from src.sql.interface.get_data.get_proteins import get_ncbi_acc_for_uniprot_acc, get_uniprot_accessions
from src.sql.interface.temp_tables import create_temp_pfam_protein_table, drop_temp_pfam_protein_table


NO_FILTER = set()
NO_TAXONOMY_FILTER = {"genera": set(), "species": set(), "strains": set()}


@pytest.fixture
def db_path(tmp_path):
    path = tmp_path / "test.db"
    sql_orm.get_db_connection(path, False, new=True)
    return path


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
    def _record(self, protein_id, accession="PF00051", start=1, end=10, interpro="IPR000001"):
        record = PfamMatchRecord(
            "P00734", accession, "Kringle domain", "domain", interpro, start, end, "109.0",
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
