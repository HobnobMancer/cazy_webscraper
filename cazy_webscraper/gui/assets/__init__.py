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
"""Store assets and additional data for the GUIs"""

from cazy_webscraper import (
    COPYRIGHT,
    DEVELOPERS,
    DOCUMENTATION,
    ISSUES_PAGE,
    LICENCE,
    REPOSITORY,
    WEBSITE,
    GUI_VERSION,
)


def build_menus(subcommand, desc):
    """Build the menu for the GUI.
    
    :param subcommand: str, name of the subcommand the gui is generated for
    :param desc: str, description of the subcommand
    
    Return list of dicts.
    """
    menu = [
        { # Define the FILE menus
            'name': 'File',
            'items': [
                {
                    'type': 'AboutDialog',
                    'menuTitle': 'About',
                    'name': subcommand,
                    'description': desc,
                    'version': GUI_VERSION,
                    'copyright': COPYRIGHT,
                    'website': WEBSITE,
                    'developer': DEVELOPERS,
                    'license': LICENCE,
                },
                {
                    'type': 'Link',
                    'menuTitle': 'Documentation',
                    'url': DOCUMENTATION,
                },
                {
                    'type': 'Link',
                    'menuTitle': 'Getting started poster',
                    'url': 'https://hobnobmancer.github.io/cazy_webscraper/getting_started_cazy_webscraper_v2.pdf',
                },
                {
                    'type': 'MessageDialog',
                    'menuTitle': 'Citation',
                    'message': 'If use use cazy_webscraper in your work, please cite our work: \n    Hobbs, Emma E. M.; Pritchard, Leighton; Chapman, Sean; Gloster, Tracey M. (2021): cazy_webscraper Microbiology Society Annual Conference 2021 poster. figshare. Poster. https://doi.org/10.6084/m9.figshare.14370860.v7',
                    'caption': 'Plesae cite our work',
                },
                {
                    'type': 'Link',
                    'menuTitle': 'Visit Our Site',
                    'url': WEBSITE,
                },
                {
                    'type': 'Link',
                    'menuTitle': 'Visit the GitHub repository',
                    'url': REPOSITORY,
                }
            ],
        },
        { # Define the HELP menu
            'name': 'Help',
            'items': [
                {
                    'type': 'Link',
                    'menuTitle': 'Documentation',
                    'url': DOCUMENTATION,
                },
                {
                    'type': 'AboutDialog',
                    'name': 'getting help',
                    'menuTitle': 'Help',
                    'description': 'If you need help, think an additional feature could be included and/or have found a bug, please raise an issue at GitHub.',
                    'website': ISSUES_PAGE,
                },
            ],
        },
        { # Define the CONTRIBUTE menu
            'name': 'Contributing',
            'items': [
                {
                    'type': 'AboutDialog',
                    'menuTitle': 'reporting bugs and errors',
                    'name': 'Reporting',
                    'description': 'If you find a bug, or an error in the code or documentation, please report this by raising an issue at the GitHub issues page for cazy_webscraper.',
                    'website': ISSUES_PAGE,
                },
                {
                    'type': 'AboutDialog',
                    'name': 'Contributing',
                    'menuTitle': 'contributing code or documentation',
                    'description': "We gratefully accept code contributions, if you would like to fix a bug, improve documentation, or extend cazy_webscraper and its subcommands please see the guidelines:",
                    'website': 'https://cazy-webscraper.readthedocs.io/en/latest/contributing.html',
                },
                {
                    'type': 'AboutDialog',
                    'name': 'Suggesting improvments',
                    'menuTitle': 'suggestions for improvement',
                    'description': "If you would like to make a suggestion for how we could improve cazy_webscraper, we welcome contributions at the GitHub issues page.",
                    'website': ISSUES_PAGE,
                },
            ],
        },
    ]

    return menu
