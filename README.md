# cazy_webscraper

-------------------------------

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4300858.svg)](https://doi.org/10.5281/zenodo.4300858)
[![licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/cazy_webscraper/blob/master/LICENSE)
[![CircleCI](https://circleci.com/gh/HobnobMancer/cazy_webscraper.svg?style=shield)](https://circleci.com/gh/HobnobMancer/cazy_webscraper)
[![codecov](https://codecov.io/gh/HobnobMancer/cazy_webscraper/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/cazy_webscraper)
[![Documentation Status](https://readthedocs.org/projects/cazy-webscraper/badge/?version=latest)](https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cazy_webscraper/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda) 
[![pyani PyPi version](https://img.shields.io/pypi/v/cazy_webscraper "PyPI version")](https://pypi.python.org/pypi/cazy_webscraper)

-------------------------------

# Attention! With the new CAZy update, cazy_webscraper is currently defunct for the scraping of CAZy. We are addressing this, and will update shortly.

`cazy_webscraper` is an application and Python3 package for the automated retrieval of protein data from the [CAZy](http://wwww.cazy.org/) database. The code is distributed under the MIT license.

`cazy_webscraper` retrieves protein data from the [CAZy database](https://www.cazy.org) into a local SQLite3 database. This enables users to integrate the dataset into analytical pipelines, and interrogate the data in a manner unachievable through the CAZy website.

Using the `expand` subcommand, a user can retrieve CAZyme protein sequence data from [GenBank](https://www.ncbi.nlm.nih.gov/genbank/), and protein structure files from the Research Collaboratory for Structural Bioinformatics (RCSB) Protein Data Bank [(PDB)](https://www.rcsb.org/).

`cazy_webscraper` can recover specified CAZy Classes and/or CAZy families, and these queries can be filtered by taxonomy at Kingdoms, genus, species or strain level. Successive CAZy queries can be collated into a single local database. A log of each query is recorded in the database for transparency, reproducibility and shareablity.

If you use `cazy_webscraper, please cite the following publication:

> Hobbs, Emma E. M.; Pritchard, Leighton; Chapman, Sean; Gloster, Tracey M. (2021): cazy_webscraper Microbiology Society Annual Conference 2021 poster. FigShare. Poster. [https://doi.org/10.6084/m9.figshare.14370860.v7](https://doi.org/10.6084/m9.figshare.14370860.v7)

## Best practice

**Please do not perform a complete scrape of the CAZy database unless you specifically require to reproduce the entire CAZy dataset. A complete scrape will take several hours and may unintentionally deny the service to others.**

When performing a series of many automated calls to a server it is best to do this when traffic is lowest, such as at weekends or overnight *at the server*.

## Documentation

Please see the [full documentation at ReadTheDocs](https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest).

### Installation

`cazy_webscraper` can be installed *via* `conda` or `pip`:

```bash
conda install -c bioconda cazy_webscraper
```

Please see the [`conda` documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) and [`bioconda` documentation](https://bioconda.github.io/) for further details.

```bash
pip install cazy_webscraper
```

Please see the [`pip` documentation](https://pypi.org/project/pip/) for further details.

### Quickstart

We have produced a "Getting Started With `cazy_webscraper`" [poster](https://hobnobmancer.github.io/cazy_webscraper/getting_started_poster.pdf).

### Creating a local CAZyme database

`--config`, `-c` - path to a configuration YAML file. Default: None.

`--classes` - list of classes from which all families are to be scrape.

`--cazy_synonyms` - path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

`--citation`, `-C` - Print the `cazy_webscraper` citation. When called, the program terminates after printng the citation and CAZy is **not** scraped.

`--db_output`, `-d` - path to write out a **new** local CAZyme database to. Include the name of the new database, including the `.db` extension. Default: None.

`--database`, `-D` - path to an **existings** local CAZyme database to add newly scraped too. Default: None.

_Do not use `--db_output` **and** `--database` at the same time._

`--dict_output` - path to write out a JSON, keyed by GenBank accession and valued by list of CAZy (sub)family annotations. Default: None.

_`--dict_output` can be called instead of `--db_output` and `--database`, but can also be called along side these two flags, producing an output database **and** JSON file._

_If `--db_output`, `--database` and `--dict_output` are **not** called, no target output is provided to `cazy_webscraper`. `cazy_webscraper` will write out a local CAZyme database to the cwd with the standardised name `cazy_webscraper_<date>_<time>.db`_

`--force`, `-f` - force overwriting existing output file. Default: False.

`--families` - List of CAZy (sub)families to scrape.

`--genera` - List of genera to restrict the scrape to. Default: None, filter not applied to scrape.

`--log`, `-l` - Target path to write out a log file. If not called, no log file is written. Default: None.

`--no_db` - When called, no local CAZyme database will be written out. Use when only a output JSON file is desired. Default: False.

_If `--db_output`, `--database` and `--dict_output` are **not** called, and `--no_db` is called no target output is provided to `cazy_webscraper`. `cazy_webscraper` will terminate and **not** scrape CAZy._

`--nodelete_cache` - When called, content in the existing cache dir will **not** be deleted. Default: False.

`--retries`, `-r` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

`--subfamilies`, `-s` - Enable retrival of CAZy subfamilies, otherwise **only** CAZy family annotations will be retrieved. Default: False.

`--species` - List of species written as Genus Species) to restrict the scraping of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

`--strains` - List of specific species strains to restrict the scraping of CAZymes to.

`--timeout`, `-t` - Connection timout limit (seconds). Default: 45.

`--verbose`, `-v` - Enable verbose logging. Default: False.

`--version`, `-V` - Print `cazy_webscraper` version number. When called and the version number is printed, `cazy_webscraper` is immediately terminated.

### Default CAZy class synonyms

Below are the defauls accepted syonyms for CAZy class names.



### Retrieving protein sequences and structure files

The `expand` subcommand is used to update a local CAZy database. It manages retrieval of CAZyme protein sequences from GenBank and protein structure files from RCSB/PDB.

### Roadmap

Our roadmap for development and improvement is shared on the repository wiki

- [`cazy_webscraper` roadmap](https://github.com/HobnobMancer/cazy_webscraper/wiki/Roadmap)

## Contributions

We welcome contributions and suggestions. You can raise issues at this repository, or fork the repository and submit pull requests, at the links below:

- [Issues](https://github.com/HobnobMancer/cazy_webscraper/issues)
- [Pull Requests](https://github.com/HobnobMancer/cazy_webscraper/pulls)

## License and copyright

MIT License

Copyright (c) 2020-2021 University of St Andrews
Copyright (c) 2020-2021 University of Strathclyde
Copyright (c) 2020-2021 UJames Hutton Institute
