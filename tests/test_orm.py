#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Tests sql_orm.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


from argparse import Namespace
from datetime import datetime

import pytest

from sqlalchemy.exc import IntegrityError

from cazy_webscraper.sql import sql_orm


def test_genbank_table():
    assert "<class 'cazy_webscraper.sql.sql_orm.Genbank'>" == str(sql_orm.Genbank)


def test_genbank_instance():
    assert "-Genbank accession=genbank_acc-" == str(sql_orm.Genbank(genbank_accession='genbank_acc'))


def test_genbank_repr():
    assert "<Class GenBank acc=genbank_acc>" == str(repr(sql_orm.Genbank(genbank_accession='genbank_acc')))


def test_ncbi_tax_table():
    assert "<class 'cazy_webscraper.sql.sql_orm.NcbiTax'>" == str(sql_orm.NcbiTax)


def test_ncbi_tax_instance():
    assert "-Ncbi Tax, Kingdom=Bacteria, genus=Trichoderma, species=None-" == str(sql_orm.NcbiTax(ncbi_tax_id=1, kingdom='Bacteria', genus='Trichoderma'))


def test_ncbi_repr():
    assert "<Class NcbiTax, Kingdom=Bacteria, genus=Trichoderma, species=None>" == str(repr(sql_orm.NcbiTax(ncbi_tax_id=1, kingdom='Bacteria', genus='Trichoderma')))


def test_tax_table():
    assert "<class 'cazy_webscraper.sql.sql_orm.Taxonomy'>" == str(sql_orm.Taxonomy)


def test_tax_instance():
    assert "-Source organism, Genus=genus, Species=species-" == str(sql_orm.Taxonomy(genus='genus', species='species'))


def test_tax_repr():
    assert "<Class Taxonomy: genus=genus, species=species, id=None, kndgm=None>" == str(repr(sql_orm.Taxonomy(genus='genus', species='species')))


def test_kingdom_table():
    assert "<class 'cazy_webscraper.sql.sql_orm.Kingdom'>" == str(sql_orm.Kingdom)


def test_kingdom_instance():
    assert "-Kingdom, kingdom=Bacteria-" == str(sql_orm.Kingdom(kingdom='Bacteria'))


def test_kingdom_repr():
    assert "<Class Kingdom, kingdom=Bacteria, kingdom_id=None>" == str(repr(sql_orm.Kingdom(kingdom='Bacteria')))


def test_ec_table():
    assert "<class 'cazy_webscraper.sql.sql_orm.Ec'>" == str(sql_orm.Ec)


def test_ec_instance():
    assert "-EC1.2.3.4-ec_id=1.2.3.4-" == str(sql_orm.Ec(ec_number='1.2.3.4'))


def test_ec_repr():
    assert "<Class EC, EC1.2.3.4, ec_id=1.2.3.4>" == str(repr(sql_orm.Ec(ec_number='1.2.3.4')))


def test_fam_table():
    assert "<class 'cazy_webscraper.sql.sql_orm.CazyFamily'>" == str(sql_orm.CazyFamily)


def test_fam_instance():
    assert "-CAZy Family, Family=GH1, Subfamily=None, id=None-" == str(sql_orm.CazyFamily(family='GH1'))


def test_fam_repr():
    assert "<Class Family, family=GH1, subfamily=None, id=None>" == str(repr(sql_orm.CazyFamily(family='GH1')))


def test_log_table():
    assert "<class 'cazy_webscraper.sql.sql_orm.Log'>" == str(sql_orm.Log)


def test_log_instance():
    assert "Log: date=date, scraped classes=GH, scraped families=GH1, cmd line commands=cmd_line" == str(sql_orm.Log(date='date', classes='GH', families='GH1', cmd_line='cmd_line'))


def test_log_repr():
    assert "Class Log: date=date, scraped classes=GH, scraped families=GH1, cmd line commands=cmd_line>" == str(repr(sql_orm.Log(date='date', classes='GH', families='GH1', cmd_line='cmd_line')))


def test_uniprot_table():
    assert "<class 'cazy_webscraper.sql.sql_orm.Uniprot'>" == str(sql_orm.Uniprot)


def test_uniprot_instance():
    assert "-Uniprot, accession=acc, name=name, id=1-" == str(sql_orm.Uniprot(uniprot_id=1, uniprot_name='name', uniprot_accession='acc'))


def test_uniprot_repr():
    assert "<Uniprot, accession=acc, name=name, id=1>" == str(repr(sql_orm.Uniprot(uniprot_id=1, uniprot_name='name', uniprot_accession='acc')))


def test_pdb_table():
    assert "<class 'cazy_webscraper.sql.sql_orm.Pdb'>" == str(sql_orm.Pdb)


def test_pdb_instance():
    assert "-PDB accession=pdb_acc, id=None-" == str(sql_orm.Pdb(pdb_accession='pdb_acc'))


def test_pdb_repr():
    assert "<Class Pdb accession=pdb_acc, id=None>" == str(repr(sql_orm.Pdb(pdb_accession='pdb_acc')))
