# cazy_webscraper

-------------------------------

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4300858.svg)](https://doi.org/10.5281/zenodo.4300858)
[![licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/cazy_webscraper/blob/master/LICENSE)
[![CircleCI](https://circleci.com/gh/HobnobMancer/cazy_webscraper.svg?style=shield)](https://circleci.com/gh/HobnobMancer/cazy_webscraper)
[![codecov](https://codecov.io/gh/HobnobMancer/cazy_webscraper/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/cazy_webscraper)
[![Documentation Status](https://readthedocs.org/projects/cazy-webscraper/badge/?version=latest)](https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cazy_webscraper/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Conda-downloads](https://img.shields.io/conda/dn/bioconda/cazy_webscraper?label=Bioconda%20downloads)](https://anaconda.org/bioconda/cazy_webscraper)  
[![pyani PyPi version](https://img.shields.io/pypi/v/cazy_webscraper "PyPI version")](https://pypi.python.org/pypi/cazy_webscraper)
[![PypiDownload](https://img.shields.io/pypi/dm/cazy_webscraper?label=Pypi%20downloads)](https://pypi.org/project/cazy-webscraper/)

-------------------------------

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

### Retrieving protein sequences and structure files

The `expand` subcommand is used to update a local CAZy database. It manages retrieval of CAZyme protein sequences from GenBank and protein structure files from RCSB/PDB.

### Roadmap

Our roadmap for development and improvement is shared on the repository wiki

- [`cazy_webscraper` roadmap](https://github.com/HobnobMancer/cazy_webscraper/wiki/Roadmap)

We welcome contributions and suggestions. You can raise issues at this repository, or fork the repository and submit pull requests, at the links below:

- [Issues](https://github.com/HobnobMancer/cazy_webscraper/issues)
- [Pull Requests](https://github.com/HobnobMancer/cazy_webscraper/pulls)

## LICENSE

MIT License

Copyright (c) 2020-2021 University of St Andrews

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE
