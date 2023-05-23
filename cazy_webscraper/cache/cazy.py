#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Cache data retrieved from the remove CAZy database"""


import logging
import json


def cache_cazy_data(cazy_data, cache_dir):
    """Cache dict of data retrieved from CAZy.

    :param cazy_data: dict of data from CAZy
    :param cache_dir: Path

    Return nothing
    """
    logger = logging.getLogger(__name__)

    # convert all sets to lists
    for gbk_acc in cazy_data:
        for key in cazy_data[gbk_acc]:
            if key == 'taxonomy':
                taxs = []
                for tax_tuple in cazy_data[gbk_acc][key]:
                    taxs.append(list(tax_tuple))
                cazy_data[gbk_acc][key] = taxs

            elif key == 'families':
                for fam in cazy_data[gbk_acc][key]:
                    cazy_data[gbk_acc][key][fam] = list(cazy_data[gbk_acc][key][fam])

            elif type(cazy_data[gbk_acc][key]) is set:
                cazy_data[gbk_acc][key] = list(cazy_data[gbk_acc][key])

    cache_path = cache_dir / "cached_cazy_data.json"

    logger.warning(f"Caching CAZy data to {cache_dir}")

    with open(cache_path, "w") as fh:
        json.dump(cazy_data, fh)
