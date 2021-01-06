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

"""Tests the module utilities which builds the logger and args parser.

These test are intened to be run from the root of the repository using:
pytest -v

"""

import pytest

from argparse import Namespace

from scraper import utilities


@pytest.fixture
def args_v_false():
    args_dict = {
        "args": Namespace(
            verbose=False,
            log=None,
        )
    }
    return args_dict


@pytest.fixture
def args_v_true():
    args_dict = {
        "args": Namespace(
            verbose=True,
            log="tests/test_outputs/test_log",
        )
    }
    return args_dict


# test building the logger

def test_verbose_false(args_v_false):
    """Test build_logger when args.verbose and args.log are false"""
    logger = utilities.build_logger("test", args_v_false["args"])
    assert logger.level == 30


def test_verbose_true(args_v_true):
    """Test build_logger when args.verbose and args.log are true"""
    logger = utilities.build_logger("test", args_v_true["args"])
    assert logger.level == 20


# test building the parser


def test_parser():
    """Test building the parser when argsv is None"""
    utilities.build_parser()


def test_parser_arsv():
    """Test building the parser when argsv is not None"""
    utilities.build_parser(["-f"])
