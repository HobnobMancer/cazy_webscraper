#!/usr/bin/env python3
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
"""Module for retrieving protein structures from PDB."""

from bioservices import PDB

from pymol import cmd


accession = "test"

# user will specify the file format at at the cmd when enabling pdb

# put in try execpt in case error is raised because cant connect
structure_file = PDB.get_file(identifier=accession, frmt=args.pdb)
cmd.load(accession, object=structure_file)
cmd.save(accession)

# need to find out how to use pymol to save the structure file
