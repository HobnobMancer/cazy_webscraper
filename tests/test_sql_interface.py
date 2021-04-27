#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Tests script that add CAZyme data to the SQL database.

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace
from datetime import datetime

from sqlalchemy.orm.exc import ObjectDeletedError

from scraper.sql import sql_interface


# # Unit tests for add_new_protein_to_db()


# def test_adding_new_protein_and_new_species(db_session, monkeypatch):
#     """Test add_new_protein_to_db and a new species to the database."""

#     def mock_return_none(*args, **kwargs):
#         return

#     monkeypatch.setattr(sql_interface, "add_ec_numbers", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_nonprimary_gbk_accessions", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_uniprot_accessions", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_return_none)

#     time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#     new_species = f"genus-{time_stamp} species-{time_stamp}"

#     sql_interface.add_new_protein_to_db(
#         "cazyme_name",
#         "GH5_1",
#         new_species,
#         "Bacteria",
#         Genbank(genbank_accession=f"primary_genbank{time_stamp}"),
#         db_session,
#         ec_numbers=["EC number", "ec number"],
#         gbk_nonprimary=["gen1", "gen2"],
#         uni_primary=["primary_uni"],
#         uni_nonprimary=["uni1", "uni2"],
#         pdb_accessions=["pdb1", "pdb2"],
#     )


# def test_adding_new_protein_and_new_species_and_fam(db_session, monkeypatch, time_stamp):
#     """Test add_new_protein_to_db and a new species to the database."""

#     def mock_return_none(*args, **kwargs):
#         return

#     monkeypatch.setattr(sql_interface, "add_ec_numbers", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_nonprimary_gbk_accessions", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_uniprot_accessions", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_return_none)

#     new_species = f"genus-{time_stamp} species-{time_stamp}"

#     try:
#         sql_interface.add_new_protein_to_db(
#             "cazyme_name_test",
#             "GH5",
#             new_species,
#             f"kingdom:{time_stamp}",
#             Genbank(genbank_accession=f"primary_genbank{time_stamp}"),
#             db_session,
#             ec_numbers=["EC number", "ec number"],
#             gbk_nonprimary=["gen1", "gen2"],
#             uni_nonprimary=["uni1", "uni2"],
#             pdb_accessions=["pdb1", "pdb2"],
#         )
#     except ObjectDeletedError as e:
#         pass


# def test_addding_new_protein_with_existing_species(db_session, monkeypatch):
#     """Test add_new_protein_to_db when the species exists in the database."""

#     def mock_return_none(*args, **kwargs):
#         return

#     monkeypatch.setattr(sql_interface, "add_ec_numbers", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_nonprimary_gbk_accessions", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_uniprot_accessions", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_return_none)

#     existing_species = "test_existing_genus test_existing_species"

#     gb = db_session.query(Genbank).filter(Genbank.genbank_accession == 'one_genbank_multi_cazymes')\
#         .all()[0]

#     sql_interface.add_new_protein_to_db(
#         cazyme_name="cazyme_name",
#         family="GH5_1",
#         source_organism=existing_species,
#         tax_kingdom='Bacteria',
#         primary_genbank_object=gb,
#         session=db_session,
#         ec_numbers=["EC number", "ec number"],
#         gbk_nonprimary=["gen1", "gen2"],
#         uni_primary=["primary_uni"],
#         uni_nonprimary=["uni1", "uni2"],
#         pdb_accessions=["pdb1", "pdb2"],
#     )


# def test_adding_new_protein_with_multiple_species(db_session, monkeypatch):
#     """Test add_new_protein_to_db when there are multiple records for the same species."""

#     def mock_return_none(*args, **kwargs):
#         return

#     monkeypatch.setattr(sql_interface, "add_ec_numbers", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_nonprimary_gbk_accessions", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_uniprot_accessions", mock_return_none)
#     monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_return_none)

#     duplicate_species = "duplicate_genus duplicate_species"
#     time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

#     sql_interface.add_new_protein_to_db(
#         "cazyme_name",
#         "GH5_1",
#         duplicate_species,
#         "Eukaryota",
#         Genbank(genbank_accession=f"primary_genbank{time_stamp}"),
#         db_session,
#         ec_numbers=["EC number", "ec number"],
#         gbk_nonprimary=["gen1", "gen2"],
#         uni_primary=["primary_uni"],
#         uni_nonprimary=["uni1", "uni2"],
#         pdb_accessions=["pdb1", "pdb2"],
#     )


# # test add_data_to_protein_record()


# def test_add_data_to_protein_record(db_session, monkeypatch):
#     """Test add_data_to_protein_record()."""
#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(sql_interface, "add_cazy_family", mock_no_return)
#     monkeypatch.setattr(sql_interface, "add_pdb_accessions", mock_no_return)

#     sql_interface.add_data_to_protein_record(
#         Cazyme(cazyme_name="test_cazyme"),
#         "GH3",
#         db_session,
#         ec_numbers=["EC number", "ec number"],
#         gbk_nonprimary=["gen1", "gen2"],
#         uni_primary=["primary_uni"],
#         uni_nonprimary=["uni1", "uni2"],
#         pdb_accessions=["pdb1", "pdb2"],
#     )


# # Unit tests for add_cazy_family


# def test_adding_new_family(db_session):
#     """Test adding a new CAZy family to the local database."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#     new_fam = f"GH{time_stamp}"

#     sql_interface.add_cazy_family(
#         new_fam,
#         cazyme,
#         db_session,
#     )


# def test_add_existing_family(db_session):
#     """Test adding an existing family to a CAZyme record."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     existing_fam = "FamOnly"

#     sql_interface.add_cazy_family(
#         existing_fam,
#         cazyme,
#         db_session,
#     )


# def test_add_new_fam_cos_old_fam_has_subfam(db_session):
#     """Test adding a new family because the existing family record is a subfamily."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     existing_fam = "FamWithSubFam"

#     sql_interface.add_cazy_family(
#         existing_fam,
#         cazyme,
#         db_session,
#     )


# def test_duplicate_families_no_nonsubfams(db_session):
#     """Test when multiple families are found with none without subfamilies."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     duplicate_families = "dupFam"

#     sql_interface.add_cazy_family(
#         duplicate_families,
#         cazyme,
#         db_session,
#     )


# def test_duplicate_families_one_with_no_subfams(db_session):
#     """test when multiple families are found and one has no subfamiy submission."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     identical_fam = "identFam"

#     sql_interface.add_cazy_family(
#         identical_fam,
#         cazyme,
#         db_session,
#     )


# def test_duplicate_families_multiple_with_no_subfams(db_session):
#     """test when multiple families have no subfamily assoication."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     identical_fam = "IdenticalFamily"

#     sql_interface.add_cazy_family(
#         identical_fam,
#         cazyme,
#         db_session,
#     )


# # Unit tests for adding subfamilies


# def test_adding_new_subfamily(db_session):
#     """Test adding a new subfamily to the local database."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     time_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
#     new_fam = f"GH{time_stamp}_1"

#     sql_interface.add_cazy_subfamily(
#         new_fam,
#         cazyme,
#         db_session,
#     )


# def test_adding_cazyme_to_existing_db(db_session):
#     """Test adding a CAZyme to an existing subfamily in the local database."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     existing_fam = "testGH1_53"

#     sql_interface.add_cazy_subfamily(
#         existing_fam,
#         cazyme,
#         db_session,
#     )


# def test_multiple_subfamilies_found(db_session):
#     """Test when multiple identical subfamilies are found in the local database."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     identical_subfam = "testident_ident123"

#     sql_interface.add_cazy_subfamily(
#         identical_subfam,
#         cazyme,
#         db_session,
#     )


# # Unit tests for adding non-primary GenBank accessions


# def test_adding_new_non_prim_gb_acc(db_session):
#     """Test adding a new non-primary GenBank accession."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     time_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
#     new_acc = [f"gb_acc_{time_stamp}"]

#     sql_interface.add_nonprimary_gbk_accessions(
#         new_acc,
#         cazyme,
#         db_session,
#     )


# def test_duplicate_non_prim_db_acc(db_session):
#     """Test when finding multiple identical non-primary duplicate GenBank accessions."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     duplicate_acc = ["dupNonPrimGBAcc"]

#     sql_interface.add_nonprimary_gbk_accessions(
#         duplicate_acc,
#         cazyme,
#         db_session,
#     )


# # Unit tests for adding EC numbers


# def test_adding_new_ec_num(db_session):
#     """Testing adding a new EC# to the local database."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     time_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
#     ec = [f"EC{time_stamp}"]

#     sql_interface.add_ec_numbers(
#         ec,
#         cazyme,
#         db_session,
#     )


# def test_adding_existing_ec_num(db_session):
#     """Testing adding an existing EC# to a CAZyme."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     ec = ["existing_ec"]

#     sql_interface.add_ec_numbers(
#         ec,
#         cazyme,
#         db_session,
#     )


# def test_finding_multiple_ecs(db_session):
#     """Testing handling when multiple duplicate EC#s are found."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     ec = ["dupEC"]

#     sql_interface.add_ec_numbers(
#         ec,
#         cazyme,
#         db_session,
#     )


# # Unit tests for adding primary UniProt accessions


# def test_new_primary_uniprot(db_session):
#     """Test adding a new primary UniProt accession."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     time_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
#     new_acc = f"uni_acc_test{time_stamp}"

#     sql_interface.add_uniprot_accessions(
#         [new_acc],
#         cazyme,
#         True,
#         db_session,
#     )


# def test_existing_primary_uniprot(db_session):
#     """Test adding an existing primary UniProt accession to a CAZyme."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     existing_uni = "existing_acc_test"

#     sql_interface.add_uniprot_accessions(
#         [existing_uni],
#         cazyme,
#         True,
#         db_session,
#     )


# def test_add_nonprimery_uniprot(db_session, time_stamp):
#     """Test adding adding non-primary UniProt accessions."""

#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]

#     sql_interface.add_uniprot_accessions(
#         [f'acc{time_stamp}', 'dupp_np'],
#         cazyme,
#         True,
#         db_session,
#     )


# # Unit tests for adding PDB accessions to the local database


# def test_adding_pdb_accessions(time_stamp, db_session):
#     """Test adding PDB accessions to the local database."""
#     cazyme = db_session.query(Cazyme).filter(Cazyme.cazyme_id == 50).all()[0]
#     accessions = [f'new_accession-{time_stamp}', 'existing_pdb', 'dupEC']

#     sql_interface.add_pdb_accessions(
#         accessions,
#         cazyme,
#         db_session,
#     )
