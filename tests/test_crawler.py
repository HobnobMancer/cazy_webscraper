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
"""Tests crawler.__init__.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


from urllib import request
from argparse import Namespace
from cazy_webscraper import crawler


class MockedResponse:
    """Mock response to CAZy"""
    def __init__(self, name):
        self.name = name

    def info(self):
        return MockContent("name")

    def read(self, content):
        return ""

    def __str__(self):
        return "Mock CAZy response"

    def __repr__(self):
        return "<Mock CAZy response>"


class MockContent:
    """Mock getting content length"""
    def __init__(self, name):
        self.name = name

    def get(self, func):
        return 0

    def __str__(self):
        return "Mock MockContent"

    def __repr__(self):
        return "<Mock MockContent>"


def test_failed_file_download(monkeypatch):
    """Test failed page download using get_cazy_file()"""
    def mock_failed_connection(*args, **kwards):
        return MockedResponse(2)

    def mock_get(*args, **kwards):
        return 0

    monkeypatch.setattr(request, "urlopen", mock_failed_connection)
    monkeypatch.setattr(crawler, "urlopen", mock_failed_connection)

    argsdict = {
        "args": Namespace(
            force=True,
            nodelete_cache=True,
            retries=1,
            timeout=1,
        )
    }

    a = crawler.get_cazy_file("tests/test_outputs/test_crawler/test_crawler", argsdict["args"], max_tries=2)
    assert a is None
