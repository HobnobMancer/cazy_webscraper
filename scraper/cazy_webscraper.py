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
"""
Web scraper to scrape CAZy website and retrieve all protein data.

:cmd_args:...

:func ...:...

"""

from collections import defaultdict

import mechanicalsoup

from tqdm.notebook import tqdm


class Protein:
    """A single protein.

    Each protein has a name, source organism (source), and links to external databases.
    The links to external databases are stored in a dictionary, keyed by the external database
    name ('str') with 'list' values becuase there may be multiple links per database.
    """

    def __init__(self, name, source, links=None):
        self.name = name
        self.source = source
        if links is None:
            self.links = defaultdict(list)
        else:
            self.links = links
    
    def __str__(self):
        """Create representative string of class object"""
        return f"{self.name} ({self.source}): links to {self.links.keys()}"

    def __repr__(self):
        """Create representative object"""
        return(
            (
                f"<Protein: {id(self)}: {self.name}, "
                f"({self.source}), {len(self.links)} to external databases>"
            )
        )


class Family:
    """A single CAZy family.
    
    Each family has a name and a set of proteins that are members of
    the family."""
    
    members = set()  # holds Protein instances

    def __init__(self, name):
        self.name = name

    def __str__(self):
        return f"CAZy family {self.name}: {len(self.members)} protein members"
    
    def __repr__(self):
        return f"<Family: {id(self)}: {self.name}, {len(self.members)} protein members"


if __name__ == "__main__":
    main()
