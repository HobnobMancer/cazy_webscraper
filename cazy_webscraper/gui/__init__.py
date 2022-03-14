#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
"""Module for creating GUIs for cazy_webscraper."""


from pathlib import Path


def build_and_covert_to_paths(gooey_args):
    """Build path to log file and convert strings from file/dir choosers to paths.
    
    :param gooey_args: arguments from gui
    
    Return arguments from gui"""
    # compile path for the log file 
    if gooey_args.log is not None and gooey_args.log_dir is not None:
        gooey_args.log = Path(gooey_args.log_dir) / gooey_args.log

    # Convert strings to Paths
    gooey_args.database = Path(gooey_args.database)

    if gooey_args.cache_dir is not None:
        gooey_args.cache_dir = Path(gooey_args.cache_dir)
    
    if gooey_args.cazy_synonyms is not None:
        gooey_args.cazy_synonyms = Path(gooey_args.cazy_synonyms)

    if gooey_args.config is not None:
        gooey_args.config = Path(gooey_args.config)
    
    return gooey_args
