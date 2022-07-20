#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
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
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""A module for expanding the the local CAZY database beyond what is provided within CAZy.

Submodules are organised by databased from which data is sourced"""


def get_chunks_gen(lst, chunk_length):
    """Separate the long list into separate chunks.

    :param lst: list to be separated into smaller lists (or chunks)
    :param chunk_length: int, the length of the lists the longer list is to be split up into

    Return a generator object containing lists.
    """
    for i in range(0, len(lst), chunk_length):
        yield lst[i:i + chunk_length]


def get_chunks_list(lst, chunk_length):
    """Separate the long list into separate chunks.

    :param lst: list to be separated into smaller lists (or chunks)
    :param chunk_length: int, the length of the lists the longer list is to be split up into

    Return a list of nested lists.
    """
    chunks = []
    for i in range(0, len(lst), chunk_length):
        chunks.append(lst[i:i + chunk_length])
    return chunks
