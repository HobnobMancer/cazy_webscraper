#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License

"""Tests the module sql which builds and interacts with an SQL database.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest
import shutil
import sys

from argparse import Namespace
from datetime import datetime
from pathlib import Path

from sqlalchemy.orm.exc import ObjectDeletedError

from scraper.utilities import file_io
from scraper.sql import sql_orm, sql_interface
from scraper.sql.sql_orm import (
    Cazyme,
    CazyFamily,
    Genbank,
    Taxonomy,
    Kingdom,
    EC,
    Uniprot,
    Pdb,
    Log,
    Cazymes_Genbanks,
)


@pytest.fixture
def output_dir(test_dir):
    path = Path("test_outputs")
    path = path / "test_sql_outputs"
    return path


@pytest.fixture
def database_args(output_dir):
    args_dict = {
        "args": Namespace(
            output=output_dir,
        )
    }
    return args_dict


@pytest.fixture
def existing_db_args(output_dir):
    db_path = output_dir / "cazy_scrape_2021-02-07--15-35-13.db"
    args_dict = {
        "args": Namespace(
            output=output_dir,
            database=db_path
        )
    }
    return args_dict


# Unit tests for sql_orm


def test_regex_search(db_session):
    """Test a regular expression search can be performed successfully."""

    cazy_class = "PL"
    class_query = db_session.query(Cazyme.cazyme_id).\
        join(CazyFamily, Cazyme.families).\
        filter(CazyFamily.family.regexp(rf"{cazy_class}\d+")).\
        all()

    assert len(class_query) == 24


def test_calling_class_instances():
    """Test calling Class instances strs and reprs."""
    cazyme = Cazyme(cazyme_name="test_cazyme")
    str(cazyme)
    repr(cazyme)

    tax = Taxonomy(genus="genus", species="species")
    str(tax)
    repr(tax)

    kng = Kingdom(kingdom="kingdom")
    str(kng)
    repr(kng)

    fam = CazyFamily(family="family")
    str(fam)
    repr(fam)

    fam = CazyFamily(family="family", subfamily="subfamily")
    str(fam)
    repr(fam)

    gb = Genbank(genbank_accession="accession")
    str(gb)
    repr(gb)

    cg = Cazymes_Genbanks(cazyme_id=1, genbank_id=2)
    str(cg)
    repr(cg)

    ec = EC(ec_number="EC1.2.3.4")
    str(ec)
    repr(ec)

    up = Uniprot(uniprot_accession="accession")
    str(up)
    repr(up)

    pdb = Pdb(pdb_accession="accession")
    str(pdb)
    repr(pdb)

    log = Log(date="date", time="time", classes="classes", families="families")
    str(log)
    repr(log)


def test_build_db():
    """Test building a database."""
    path_ = Path("tests")
    path_ = path_ / "test_outputs" / "test_outputs_sql" / "temp_dir"
    file_io.make_output_directory(path_, True, False)
    args = {"args": Namespace(output=path_)}
    sql_orm.build_db("time_stamp", args["args"])
    shutil.rmtree(path_)


def test_build_db_stdout():
    """Test building a database when output is stdout"""
    path_ = Path("tests")
    path_ = path_ / "test_outputs" / "test_outputs_sql" / "temp_dir"
    file_io.make_output_directory(path_, True, False)
    args = {"args": Namespace(output=sys.stdout)}
    sql_orm.build_db("time_stamp", args["args"])
    shutil.rmtree(path_)


# Unit tests for sql_interface


def test_custom_error():
    """Test custom error."""
    with pytest.raises(sql_interface.SqlInterfaceException) as pytest_wrapped_err:
        raise sql_interface.SqlInterfaceException("message")
    assert pytest_wrapped_err.type == sql_interface.SqlInterfaceException


def test_adding_a_new_protein(db_session):
    """Test adding a new protein to the local database."""
    time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    try:
        sql_interface.add_protein_to_db(
            "test_cazyme",
            "test_fam",
            "test_genus test_species",
            "kingdom",
            time_stamp,
            db_session,
            ["EC1.2.3.4"],
            ["Genbank1", "Genbank2"],
            ["Uni1P"],
            ["UniNP1"],
            ["PDB1", "PDB2"],
        )
    except ObjectDeletedError as e:
        pass


def test_add_data_to_an_existing_record_in_db(db_session):
    """Test adding data to an existing CAZyme in the local database."""

    sql_interface.add_protein_to_db(
        "JCM19301_1832",
        "PL28",
        "pallidilutea JCM 19301",
        "Bacteria",
        "GAL67220.1",
        db_session,
        ec_numbers=["EC4.2.2.-"],
        gbk_nonprimary=["Genbank1", "Genbank2"],
        uni_nonprimary=["Uni1", "Uni2"],
    )


def test_genbank_no_cazymes(db_session, monkeypatch):
    """test adding a protein to a database, GenBank accession is found with no linked CAZymes."""

    def mock_adding_a_new_protein(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_new_protein_to_db", mock_adding_a_new_protein)

    existing_genbank_with_no_cazyme = "test_genbank_no_cazyme"

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        existing_genbank_with_no_cazyme,
        db_session,
        ec_numbers=["EC4.2.2.-"],
    )


def test_genbank_no_cazymes_raise_error(db_session, monkeypatch):
    """test adding a protein to a database, GenBank accession is found with no linked CAZymes."""

    def mock_adding_a_new_protein(*args, **kwargs):
        return "error_message"

    monkeypatch.setattr(sql_interface, "add_new_protein_to_db", mock_adding_a_new_protein)

    existing_genbank_with_no_cazyme = "test_genbank_no_cazyme"

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        existing_genbank_with_no_cazyme,
        db_session,
        ec_numbers=["EC4.2.2.-"],
        gbk_nonprimary=["Genbank1", "Genbank2"],
        uni_nonprimary=["Uni1", "Uni2"],
    )


def test_one_genbank_multiple_cazymes(db_session, monkeypatch):
    """test adding protein to db when genbank is found linked to muliple CAZymes."""

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_data_to_protein_record", mock_add_protein_to_db)

    genbank_with_multiple_cazymes = "one_genbank_multi_cazymes"

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        genbank_with_multiple_cazymes,
        db_session,
    )


def test_multiple_genbanks_multiple_cazymes(db_session, monkeypatch):
    """test adding protein to db when finding multiple identical CAZymes and GenBank accesisons."""

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_data_to_protein_record", mock_add_protein_to_db)

    identical_genbank_accession = "identical_accession"

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        identical_genbank_accession,
        db_session,
    )


def test_multiple_genbanks_one_cazyme(db_session, monkeypatch):
    """test adding protien to db when identical GenBank accessions with one CAZyme link."""

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_data_to_protein_record", mock_add_protein_to_db)

    multiple_accession = "multiple_accession"

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        multiple_accession,
        db_session,
    )


def test_multiple_genbanks_no_cazymes(db_session, monkeypatch):
    """test adding protien to db when identical GenBank accessions are linked to no CAZymes."""

    identical_genbanks_no_cazymes = "multi_accession_no_caz"

    def mock_add_protein_to_db(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_new_protein_to_db", mock_add_protein_to_db)

    sql_interface.add_protein_to_db(
        "test_cazyme_name",
        "cazy_family",
        "source_genus organism",
        "kingdom",
        identical_genbanks_no_cazymes,
        db_session,
    )


# Unit tests for add_new_protein_to_db()


def test_adding_new_protein_and_new_species(db_session, monkeypatch):
    """Test add_new_protein_to_db and a new species to the database."""

    def mock_return_none(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_ec_numbers", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_nonprimary_gbk_accessions", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_uniprot_accessions", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_return_none)

    time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    new_species = f"genus-{time_stamp} species-{time_stamp}"

    sql_interface.add_new_protein_to_db(
        "cazyme_name",
        "GH5_1",
        new_species,
        "Bacteria",
        Genbank(genbank_accession=f"primary_genbank{time_stamp}"),
        db_session,
        ec_numbers=["EC number", "ec number"],
        gbk_nonprimary=["gen1", "gen2"],
        uni_primary=["primary_uni"],
        uni_nonprimary=["uni1", "uni2"],
        pdb_accessions=["pdb1", "pdb2"],
    )


def test_adding_new_protein_and_new_species_and_fam(db_session, monkeypatch, time_stamp):
    """Test add_new_protein_to_db and a new species to the database."""

    def mock_return_none(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_ec_numbers", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_nonprimary_gbk_accessions", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_uniprot_accessions", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_return_none)

    new_species = f"genus-{time_stamp} species-{time_stamp}"

    try:
        sql_interface.add_new_protein_to_db(
            "cazyme_name_test",
            "GH5",
            new_species,
            f"kingdom:{time_stamp}",
            Genbank(genbank_accession=f"primary_genbank{time_stamp}"),
            db_session,
            ec_numbers=["EC number", "ec number"],
            gbk_nonprimary=["gen1", "gen2"],
            uni_nonprimary=["uni1", "uni2"],
            pdb_accessions=["pdb1", "pdb2"],
        )
    except ObjectDeletedError as e:
        pass


def test_addding_new_protein_with_existing_species(db_session, monkeypatch):
    """Test add_new_protein_to_db when the species exists in the database."""

    def mock_return_none(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_ec_numbers", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_nonprimary_gbk_accessions", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_uniprot_accessions", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_return_none)

    existing_species = "test_existing_genus test_existing_species"

    gb = db_session.query(Genbank).filter(Genbank.genbank_accession == 'one_genbank_multi_cazymes')\
        .all()[0]

    sql_interface.add_new_protein_to_db(
        cazyme_name="cazyme_name",
        family="GH5_1",
        source_organism=existing_species,
        tax_kingdom='Bacteria',
        primary_genbank_object=gb,
        session=db_session,
        ec_numbers=["EC number", "ec number"],
        gbk_nonprimary=["gen1", "gen2"],
        uni_primary=["primary_uni"],
        uni_nonprimary=["uni1", "uni2"],
        pdb_accessions=["pdb1", "pdb2"],
    )


def test_adding_new_protein_with_multiple_species(db_session, monkeypatch):
    """Test add_new_protein_to_db when there are multiple records for the same species."""

    def mock_return_none(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_ec_numbers", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_nonprimary_gbk_accessions", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_uniprot_accessions", mock_return_none)
    monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_return_none)

    duplicate_species = "duplicate_genus duplicate_species"
    time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    sql_interface.add_new_protein_to_db(
        "cazyme_name",
        "GH5_1",
        duplicate_species,
        "Eukaryota",
        Genbank(genbank_accession=f"primary_genbank{time_stamp}"),
        db_session,
        ec_numbers=["EC number", "ec number"],
        gbk_nonprimary=["gen1", "gen2"],
        uni_primary=["primary_uni"],
        uni_nonprimary=["uni1", "uni2"],
        pdb_accessions=["pdb1", "pdb2"],
    )


# test add_data_to_protein_record()


def test_add_data_to_protein_record(db_session, monkeypatch):
    """Test add_data_to_protein_record()."""
    def mock_no_return(*args, **kwargs):
        return

    monkeypatch.setattr(sql_interface, "add_cazy_family", mock_no_return)
    monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_no_return)

    sql_interface.add_data_to_protein_record(
        Cazyme(cazyme_name="test_cazyme"),
        "GH3",
        db_session,
        ec_numbers=["EC number", "ec number"],
        gbk_nonprimary=["gen1", "gen2"],
        uni_primary=["primary_uni"],
        uni_nonprimary=["uni1", "uni2"],
        pdb_accessions=["pdb1", "pdb2"],
    )


# Unit tests for add_cazy_family


def test_adding_new_family(db_session):
    """Test adding a new CAZy family to the local database."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    new_fam = f"GH{time_stamp}"

    sql_interface.add_cazy_family(
        new_fam,
        cazyme,
        db_session,
    )


def test_add_existing_family(db_session):
    """Test adding an existing family to a CAZyme record."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    existing_fam = "FamOnly"

    sql_interface.add_cazy_family(
        existing_fam,
        cazyme,
        db_session,
    )


def test_add_new_fam_cos_old_fam_has_subfam(db_session):
    """Test adding a new family because the existing family record is a subfamily."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    existing_fam = "FamWithSubFam"

    sql_interface.add_cazy_family(
        existing_fam,
        cazyme,
        db_session,
    )


def test_duplicate_families_no_nonsubfams(db_session):
    """Test when multiple families are found with none without subfamilies."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    duplicate_families = "dupFam"

    sql_interface.add_cazy_family(
        duplicate_families,
        cazyme,
        db_session,
    )


def test_duplicate_families_one_with_no_subfams(db_session):
    """test when multiple families are found and one has no subfamiy submission."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    identical_fam = "identFam"

    sql_interface.add_cazy_family(
        identical_fam,
        cazyme,
        db_session,
    )


def test_duplicate_families_multiple_with_no_subfams(db_session):
    """test when multiple families have no subfamily assoication."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    identical_fam = "IdenticalFamily"

    sql_interface.add_cazy_family(
        identical_fam,
        cazyme,
        db_session,
    )


# Unit tests for adding subfamilies


def test_adding_new_subfamily(db_session):
    """Test adding a new subfamily to the local database."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    time_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    new_fam = f"GH{time_stamp}_1"

    sql_interface.add_cazy_subfamily(
        new_fam,
        cazyme,
        db_session,
    )


def test_adding_cazyme_to_existing_db(db_session):
    """Test adding a CAZyme to an existing subfamily in the local database."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    existing_fam = "testGH1_53"

    sql_interface.add_cazy_subfamily(
        existing_fam,
        cazyme,
        db_session,
    )


def test_multiple_subfamilies_found(db_session):
    """Test when multiple identical subfamilies are found in the local database."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    identical_subfam = "testident_ident123"

    sql_interface.add_cazy_subfamily(
        identical_subfam,
        cazyme,
        db_session,
    )


# Unit tests for adding non-primary GenBank accessions


def test_adding_new_non_prim_gb_acc(db_session):
    """Test adding a new non-primary GenBank accession."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    time_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    new_acc = [f"gb_acc_{time_stamp}"]

    sql_interface.add_nonprimary_gbk_accessions(
        new_acc,
        cazyme,
        db_session,
    )


def test_duplicate_non_prim_db_acc(db_session):
    """Test when finding multiple identical non-primary duplicate GenBank accessions."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    duplicate_acc = ["dupNonPrimGBAcc"]

    sql_interface.add_nonprimary_gbk_accessions(
        duplicate_acc,
        cazyme,
        db_session,
    )


# Unit tests for adding EC numbers


def test_adding_new_ec_num(db_session):
    """Testing adding a new EC# to the local database."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    time_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    ec = [f"EC{time_stamp}"]

    sql_interface.add_ec_numbers(
        ec,
        cazyme,
        db_session,
    )


def test_adding_existing_ec_num(db_session):
    """Testing adding an existing EC# to a CAZyme."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    ec = ["existing_ec"]

    sql_interface.add_ec_numbers(
        ec,
        cazyme,
        db_session,
    )


def test_finding_multiple_ecs(db_session):
    """Testing handling when multiple duplicate EC#s are found."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    ec = ["dupEC"]

    sql_interface.add_ec_numbers(
        ec,
        cazyme,
        db_session,
    )


# Unit tests for adding primary UniProt accessions


def test_new_primary_uniprot(db_session):
    """Test adding a new primary UniProt accession."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    time_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    new_acc = f"uni_acc_test{time_stamp}"

    sql_interface.add_uniprot_accessions(
        [new_acc],
        cazyme,
        True,
        db_session,
    )


def test_existing_primary_uniprot(db_session):
    """Test adding an existing primary UniProt accession to a CAZyme."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    existing_uni = "existing_acc_test"

    sql_interface.add_uniprot_accessions(
        [existing_uni],
        cazyme,
        True,
        db_session,
    )


def test_add_nonprimery_uniprot(db_session, time_stamp):
    """Test adding adding non-primary UniProt accessions."""

    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

    sql_interface.add_uniprot_accessions(
        [f'acc{time_stamp}', 'dupp_np'],
        cazyme,
        True,
        db_session,
    )


# Unit tests for adding PDB accessions to the local database


def test_adding_pdb_accessions(time_stamp, db_session):
    """Test adding PDB accessions to the local database."""
    cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]
    accessions = [f'new_accession-{time_stamp}', 'existing_pdb', 'dupEC']

    sql_interface.add_pdb_accessions(
        accessions,
        cazyme,
        db_session,
    )
