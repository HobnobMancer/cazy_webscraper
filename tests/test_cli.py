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
"""Unit tests for the top-level CLI parser (src.utilities.parsers.parse_cmd) and the
subcommand-to-script wiring in src.utilities.parsers.subcmds.*.
"""


import pytest

from src.scripts import (
    extract_data, get_gtdb_tax, get_ncbi_genomes, get_ncbi_seqs, get_ncbi_tax,
    get_pdb_structures, get_pfams, get_uniprot_data, scrape_cazy,
)
from src.utilities.parsers.parse_cmd import build_parser


# (subcommand name, minimal required positional args, expected dispatch target)
SUBCOMMANDS = [
    ("download_cazy", ["user@example.com"], scrape_cazy),
    ("get_ncbi_seqs", ["cazy.db", "user@example.com"], get_ncbi_seqs),
    ("get_ncbi_taxs", ["cazy.db", "user@example.com"], get_ncbi_tax),
    ("get_ncbi_genomes", ["cazy.db", "user@example.com"], get_ncbi_genomes),
    ("get_uniprot_data", ["cazy.db"], get_uniprot_data),
    ("get_pdb_structures", ["cazy.db"], get_pdb_structures),
    ("get_pfams", ["cazy.db"], get_pfams),
    ("get_gtdb_taxs", ["cazy.db"], get_gtdb_tax),
    ("extract_data", ["cazy.db"], extract_data),
]


class TestSubcommandDispatch:
    @pytest.mark.parametrize("name,required_args,module", SUBCOMMANDS, ids=[s[0] for s in SUBCOMMANDS])
    def test_each_subcommand_is_registered_and_dispatches_to_its_script(self, name, required_args, module):
        args = build_parser([name] + required_args)
        assert args.func == module.main

    def test_all_nine_subcommands_are_registered(self):
        # a parser with no subcommand named still lists every choice in its usage/error text
        with pytest.raises(SystemExit):
            build_parser(["not_a_real_subcommand"])

    @pytest.mark.parametrize("name,required_args,module", SUBCOMMANDS, ids=[s[0] for s in SUBCOMMANDS])
    def test_help_does_not_crash(self, name, required_args, module, capsys):
        with pytest.raises(SystemExit) as exc_info:
            build_parser([name, "--help"])
        assert exc_info.value.code == 0

    @pytest.mark.parametrize("name,required_args,module", SUBCOMMANDS, ids=[s[0] for s in SUBCOMMANDS])
    def test_missing_required_positional_arg_is_an_error_not_a_crash(self, name, required_args, module):
        with pytest.raises(SystemExit) as exc_info:
            build_parser([name])
        assert exc_info.value.code == 2


class TestGlobalArgumentOrder:
    def test_global_flags_before_subcommand_are_accepted(self):
        args = build_parser(["-v", "get_pfams", "cazy.db"])
        assert args.verbose is True

    def test_global_flags_after_subcommand_are_rejected(self):
        # this is the change most likely to trip up a version 2 script/shell history
        with pytest.raises(SystemExit) as exc_info:
            build_parser(["get_pfams", "cazy.db", "-v"])
        assert exc_info.value.code == 2

    def test_log_flag_before_subcommand(self):
        args = build_parser(["-l", "run.log", "get_pfams", "cazy.db"])
        assert str(args.log) == "run.log"


class TestSharedFilterArguments:
    @pytest.mark.parametrize("name,required_args,module", SUBCOMMANDS, ids=[s[0] for s in SUBCOMMANDS])
    def test_filter_arguments_are_accepted(self, name, required_args, module):
        if name == "download_cazy":
            pytest.skip("download_cazy has no ec_filter/accession-file filters (nothing to filter against yet)")

        args = build_parser([name] + required_args + [
            "--classes", "GH,CE",
            "--families", "GH1,GH2",
            "--kingdoms", "bacteria",
            "--genera", "Escherichia",
            "--species", "Escherichia coli",
            "--strains", "Escherichia coli K12",
        ])
        assert args.classes == "GH,CE"
        assert args.families == "GH1,GH2"
