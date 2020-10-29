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


def browser_decorator(func):
    """Decorator to retry the wrapped function up to 'retries' times."""

    def wrapper(*args, retries=10, **kwargs):
        tries, success = 0, False
        while not success and (tries < retries):
            response = func(*args, **kwargs)
            if str(response) == "<Response [200]>":  # response was successful
                success = True
            # if response from webpage was not successful
            tries += 1
        if not success:
            print("Ran out of connection retries, and was unable to connect")
            return None
        else:
            return response

    return wrapper


@browser_decorator
def get_page(url):
    """Create browser and use browser to retrieve page for given URL.

    :param url: str, url to webpage

    Return browser response object (the page).
    """
    # create browser object
    browser = mechanicalsoup.Browser()

    # create response object
    page = browser.get(url)

    return page


if __name__ == "__main__":
    main()
