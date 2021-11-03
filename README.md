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

> `cazy_webscraper` version 1 is depracted. Please ensure you are using `cazy_webscraper` version 2 or newer.

`cazy_webscraper` is an application and Python3 package for the automated retrieval of protein data from the [CAZy](http://wwww.cazy.org/) database. The code is distributed under the MIT license.

`cazy_webscraper` retrieves protein data from the [CAZy database](https://www.cazy.org) into a local SQLite3 database. This enables users to integrate the dataset into analytical pipelines, and interrogate the data in a manner unachievable through the CAZy website.

Using the `expand` subcommand, a user can retrieve:
- CAZyme protein sequence data from [GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
- Protein structure files from the Research Collaboratory for Structural Bioinformatics (RCSB) Protein Data Bank [(PDB)](https://www.rcsb.org/)
- EC number and Uniprot protein IDs from the [UniProtKB database](https://www.uniprot.org/)

`cazy_webscraper` can recover specified CAZy Classes and/or CAZy families. These queries can be filtered by taxonomy at Kingdoms, genus, species or strain level. Successive CAZy queries can be collated into a single local database. A log of each query is recorded in the database for transparency, reproducibility and shareablity.

## Citation

If you use `cazy_webscraper`, please cite the following publication:

> Hobbs, Emma E. M.; Pritchard, Leighton; Chapman, Sean; Gloster, Tracey M. (2021): cazy_webscraper Microbiology Society Annual Conference 2021 poster. FigShare. Poster. [https://doi.org/10.6084/m9.figshare.14370860.v7](https://doi.org/10.6084/m9.figshare.14370860.v7)

## Table of Contents
<!-- TOC -->
- [`cazy_webscraper`](#cazy_webscraper)
- [Citation](#citation)
- [Best practice](#best-practice)
- [Documentation](#documentation)
    - [Installation](#installation)
    - [Quick start](#quick-start)
- [Creating a local CAZyme database](#creating-a-local-cazyme-database)
    - [Combining configuration filters](#combining-configuration-filters)
    - [Default CAZy class synonyms](#default-cazy-class-synonyms)
- [Retrieve data from UniProt](#retrieve-data-from-uniprot)
    - [Configuring UniProt data retrieval](#configuring-uniprot-data-retrieval)
- [Retrieving protein sequences from GenBank](#retrieving-protein-sequences-from-genbank)
    - [Configuring GenBank protein sequence data retrieval](#configuring-genbank-protein-sequence-retrieval)
- [Extracing protein sequences from the local CAZyme database and building a BLAST database](#extracting-protein-sequences-from-the-local-cazyme-database-and-building-a-blast-database)
- [Retrieving protein structure files from PDB](#retrieving-protein-structure-files-from-pdb)
    - [Configuring PDB protein structure file retrieval](#configuring-pdb-protein-structure-file-retrieval)
- [Configuring `cazy_webscraper` using a YAML file](#configuring-using-a-yaml-file)
- [Contributions](#contributions)
- [License and copyright](#license-and-copyright)
<!-- /TOC -->

## Best practice

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

To download all of CAZy and save the database in the default location (the cwd) with the default name (`cazy_webscraper_<date>_<time>.db`) use the following command:  
```bash
cazy_webscraper <user_email>
```

## Creating a local CAZyme database
Command line options for `cazy_webscraper`, which is used to scrape CAZy and compile a local SQLite database. 
Options are written in alphabetical order.

`email` - \[REQUIRED\] User email address. This is required by NCBI Entrez for querying the Entrez server.

`--cache_dir` - Path to cache dir to be used instead of default cache dir path.

`--cazy_synonyms` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

`--config`, `-c` - Path to a configuration YAML file. Default: None.

`--citation`, `-C` - Print the `cazy_webscraper` citation. When called, the program terminates after printng the citation and CAZy is **not** scraped.

`--classes` - list of classes from which all families are to be scrape.

`--database`, `-D` - Path to an **existings** local CAZyme database to add newly scraped too. Default: None.

`--db_output`, `-d` - Path to write out a **new** local CAZyme database to. Include the name of the new database, including the `.db` extension. Default: None.

_Do not use `--db_output` **and** `--database` at the same time._

_If `--db_output` **and** `--database` are **not** called, `cazy_webscraper` write out a local CAZyme database to the cwd with the standardised name `cazy_webscraper_<date>_<time>.db`_

`--families` - List of CAZy (sub)families to scrape.

`--force`, `-f` - force overwriting existing output file. Default: False.

`--genera` - List of genera to restrict the scrape to. Default: None, filter not applied to scrape.

`--log`, `-l` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

`--nodelete_cache` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

`--nodelete_log` - When called, content in the existing log dir will **not** be deleted. Default: False (existing content is deleted).

`--retries`, `-r` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

`--sql_echo` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

`--subfamilies`, `-s` - Enable retrival of CAZy subfamilies, otherwise **only** CAZy family annotations will be retrieved. Default: False.

`--species` - List of species written as Genus Species) to restrict the scraping of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

`--strains` - List of specific species strains to restrict the scraping of CAZymes to.

`--timeout`, `-t` - Connection timout limit (seconds). Default: 45.

`--validate`, - Retrieve CAZy family population sizes from the CAZy website and check against the number of family members added to the local CAZyme database, as a method for validating the complete retrieval of CAZy data.

`--verbose`, `-v` - Enable verbose logging. This does not set the SQLite engine `echo` parameter to True. Default: False.

`--version`, `-V` - Print `cazy_webscraper` version number. When called and the version number is printed, `cazy_webscraper` is immediately terminated.

### Combining configuration filters

`cazy_webscraper` applies filters in a successive and layered structure.

CAZy class and family filters are applied first.

Kingdom filters are applied second.

Lastly, taxonomy (genus, species and strain) filters are applied.

### Default CAZy class synonyms

CAZy classes are accepted in the written long form (such as Glycoside Hydrolases) and in their abbreviated form (e.g. GH).

Both the plural and singular abbreviated form of a CAZy class name is accepted, e.g. GH and GHs.

Spaces, hythens, underscores and no space or extract character can be used in the CAZy class names. Therefore, Glycoside Hydrolases, Glycoside-Hydrolases, Glycoside_Hydrolases and GlycosideHydrolases are all accepted.

Class names can be written in all upper case, all lower case, or mixed case, such as GLYCOSIDE-HYDROLASES, glycoside hydrolases and Glycoside Hydrolases. All lower or all upper case CAZy class name abbreviations (such as GH and gh) are accepted.

## Retrieve data from UniProt

[UniProtKB] is one of the largest protein database, incorporating data from the [PDB] structure database and other protein annotation databases.

`cazy_webscraper` can retrieve protein data from UniProt for proteins catalogued in a local CAZyme database created using `cazy_webscraper`. Specifically, for each protein, `cazy_webscraper` can retrieve:
- The UniProt accession
- PDB accessions of associated structure files from the PDB database
- EC number annnotations
- Protein sequence from the UniProt

`cazy_webscraper` always retrieves the UniProt accession, but the retrieval of PDB accession, EC numbers and protein sequences is optional.

Data can be retrieived for all proteins in the local CAZyme database, or a specific subset. CAZy class, CAZy family, genus, species, strains, kingdom and EC number filters can be defined in order to define a dataset to retrieve protein data from UniProt for.

To retrieve all UniProt data for all proteins in a local CAZyme datbase, using the following command:
```bash
cw_get_uniprot_data <path_to_local_CAZyme_db> --ec --pdb --seq
```

### Configuring UniProt data retrieval

Below are listed the command-line flags for configuring the retrieval of UniProt data.

The first positional argument is the path to the local CAZyme database, which is **required**.

`--bioservices_batch_size` - Change the query batch size submitted via [`bioservices`]() to UniProt to retrieve protein data. Default is 150. `bioservices` recommands queries not larger than 200 objects.

`--cache_dir` - Path to cache dir to be used instead of default cache dir path.

`--cazy_synonyms` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

`--config`, `-c` - Path to a configuration YAML file. Default: None.

`--classes` - List of classes from which all families are to be scrape.

`--ec`, `-e` - Enable retrieval of EC number annotations from UniProt

`--ec_filter` - Limist retrieval of protein data to proteins annotated with a provided list of EC numbers. Separate the EC numbers bu single commas without spaces. Recommend to wrap the entire str in quotation marks, for example:
```bash
cq_get_uniprot_data my_cazyme_db/cazyme_db.db --ec_filter 'EC1.2.3.4,EC2.3.1.-'
```

`--families` - List of CAZy (sub)families to scrape.

`--genera` - List of genera to restrict the scrape to. Default: None, filter not applied to scrape.

`--log`, `-l` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

`--nodelete_cache` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

`--retries`, `-r` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

`--sql_echo` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

`--species` - List of species written as Genus Species) to restrict the scraping of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

`--strains` - List of specific species strains to restrict the scraping of CAZymes to.

`--timeout`, `-t` - Connection timout limit (seconds). Default: 45.

`--uniprot_batch_size` - Size of an individual batch query submitted to the [UniProt REST API]() to retrieve the UniProt accessions of proteins identified by the GenBank accession. Default is 150. The UniProt API documentation recommands batch sizes of less than 20,000 but batch sizes of 1,000 often result in HTTP 400 errors. It is recommend to keep batch sizes less than 1,000, and ideally less than 200.

`--update_seq` - If a newer version of the protein sequence is available, overwrite the existing sequence for the protein in the database. Default is false, the protein sequence is **not** overwritten and updated.

`--verbose`, `-v` - Enable verbose logging. This does not set the SQLite engine `echo` parameter to True. Default: False.

## Retrieveing protein seqences from GenBank

Protein amino acid sequences can be retrieved for proteins in a local CAZyme database using `cazy_webscraper`. Protein sequences can be retrieved for a specific subset of proteins, identified through the use of CAZy class, CAZy family, taxonomy (kingdom, genus, species and strain) filters, and EC number filters. The retrieved protein sequences are written to the local CAZyme database.

_Extracting protein sequences from the local CAZyme database and writing them to a BLAST database and/or FASTA file(s) is covered in the next section._

To retrieve all GenBank protein seuqneces for all proteins in a local CAZyme datbase, using the following command:
```bash
cw_get_genbank_seq <path_to_local_CAZyme_db>
```

### Configuring GenBank protein sequence retrieval

Below are listed the command-line flags for configuring the retrieval of protein sequences from GenBank.

The first positional argument is the path to the local CAZyme database, which is **required**.

`--cache_dir` - Path to cache dir to be used instead of default cache dir path.

`--cazy_synonyms` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

`--config`, `-c` - Path to a configuration YAML file. Default: None.

`--classes` - List of classes from which all families are to be scrape.

`--ec_filter` - Limist retrieval of protein data to proteins annotated with a provided list of EC numbers. Separate the EC numbers bu single commas without spaces. Recommend to wrap the entire str in quotation marks, for example:
```bash
cq_get_uniprot_data my_cazyme_db/cazyme_db.db --ec_filter 'EC1.2.3.4,EC2.3.1.-'
```

`--entrez_batch_size` - Change the query batch size submitted via [`Entrez`]() to retrieve protein sequences from GenBank data. Default is 150. `Entrez` recommands queries not larger than XXX objects in length.

`--families` - List of CAZy (sub)families to scrape.

`--genera` - List of genera to restrict the scrape to. Default: None, filter not applied to scrape.

`--log`, `-l` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

`--nodelete_cache` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

`--retries`, `-r` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

`--sql_echo` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

`--species` - List of species written as Genus Species) to restrict the scraping of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

`--strains` - List of specific species strains to restrict the scraping of CAZymes to.

`--timeout`, `-t` - Connection timout limit (seconds). Default: 45.

`--update_seq` - If a newer version of the protein sequence is available, overwrite the existing sequence for the protein in the database. Default is false, the protein sequence is **not** overwritten and updated.

`--verbose`, `-v` - Enable verbose logging. This does not set the SQLite engine `echo` parameter to True. Default: False.

## Extract protein sequences from the local CAZyme database and building a BLAST database

Protein sequences from GenBank and UniProt that are stored in the local CAZyme database can be extracted using `cazy_webscraper`, and written to:
- 1 FASTA file per unique protein
- A single FASTA file containing all extracted seqences
- A BLAST database

**FASTA file format:** The protein ID line in the FASTA files compiled by `cazy_webscraper`

## Retrieving protein structure files from PDB

`cazy_webscraper` can retrieve protein structure files for proteins catalogued in a local CAZyme database. Structure files can be retrieved for all proteins in the database or a subset of proteins, chosen by defining CAZy class, CAZy family, taxonomy (kingdom, genus, species and strain) filters, and EC number filters.

Retrieval of structure files from PDB is performed by the `BioPython` module `PDB` [Cock _et al._, 2009], which writes the downloaded structure files to the local disk. Therefore, the downloaded structure files are **not** stored in the local CAZyme database at the present.

> Cock, P. J. A, Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A. _et al._ (2009) 'Biopython: freely available Python tools for computaitonal molecular biology and bioinformatics', _Bioinformatics_, 25(11), pp. 1422-3.

To retrieve structure files for all proteins in a local CAZyme database in `mmCif` and `pdb` format, use the following command:
```bash
cw_get_pdb_structures <path_to_local_CAZyme_db> mmcif,pdb
```

### Configuring PDB protein structure file retrieval

Below are listed the command-line flags for configuring the retrieval of protein structure files from PDB.

The first positional argument is the path to the local CAZyme database, which is **required**.

The second positional argument (which is also **required**) is the file types to be retrieved from PDB. The following file types are supported:  
- `mmCif`
- `pdb`
- `xml`
- `mmft`
- `bundle`
To chose multiple file types, list all desired file types, separting the files using a single comma. For example:
```bash
cw_get_genbank_seq my_cazyme_db/cazyme_db.db mmcif,pdb,xml
```
Providing the file types is **not** case sensitive, and the order the file types are listed does **not** matter.

Optional flags are listed below.

`--cache_dir` - Path to cache dir to be used instead of default cache dir path.

`--cazy_synonyms` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

`--config`, `-c` - Path to a configuration YAML file. Default: None.

`--classes` - List of classes from which all families are to be scrape.

`--ec_filter` - Limist retrieval of protein data to proteins annotated with a provided list of EC numbers. Separate the EC numbers bu single commas without spaces. Recommend to wrap the entire str in quotation marks, for example:
```bash
cq_get_uniprot_data my_cazyme_db/cazyme_db.db --ec_filter 'EC1.2.3.4,EC2.3.1.-'
```

`--entrez_batch_size` - Change the query batch size submitted via [`Entrez`]() to retrieve protein sequences from GenBank data. Default is 150. `Entrez` recommands queries not larger than XXX objects in length.

`--families` - List of CAZy (sub)families to scrape.

`--genera` - List of genera to restrict the scrape to. Default: None, filter not applied to scrape.

`--log`, `-l` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

`--nodelete_cache` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

`--outdir`, `-o` - Output directory to write out downloaded protein structure files to. Default is to write out the downloaded structure files to the current working directory.

`--retries`, `-r` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

`--sql_echo` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

`--species` - List of species written as Genus Species) to restrict the scraping of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

`--strains` - List of specific species strains to restrict the scraping of CAZymes to.

`--timeout`, `-t` - Connection timout limit (seconds). Default: 45.

`--update_seq` - If a newer version of the protein sequence is available, overwrite the existing sequence for the protein in the database. Default is false, the protein sequence is **not** overwritten and updated.

`--verbose`, `-v` - Enable verbose logging. This does not set the SQLite engine `echo` parameter to True. Default: False.

## Configuring using a YAML file

The retrieval of data from CAZy, UniProt, GenBank and PDB can be configured at the command-line **and** via a YAML file.

The YAML file must have the following structure, specifically the YAML file must have the exact keys presented below and the values can be customised to configure the behaviour of `cazy_webscraper`:
```yaml
classes:  # classes from which all proteins will be retrieved
  - "GH"
  - "CE"
Glycoside Hydrolases (GHs):
GlycosylTransferases (GTs):
Polysaccharide Lyases (PLs):
  - "GT1"
  - "GT5"
  - "GT6"
Carbohydrate Esterases (CEs):
Auxiliary Activities (AAs):
Carbohydrate-Binding Modules (CBMs):
genera:  # list genera to be scraped
 - "Trichoderma"
 - "Aspergillus"
species:  # list species, this will scrape all strains under the species
- "Pythium ultimum"
strains:  # list specific strains to be scraped
kingdoms:  # Archaea, Bacteria, Eukaryota, Viruses, Unclassified
```

For configuring the retrieval of data from UniProt, GenBank and PDB (_but not CAZy) the additional `ec` tag can be included to limit the retrieval of data to proteins annotated with specific EC numbers.

When listing EC numbers, the 'EC' prefix can be included or excluded. For example, 'EC1.2.3.4' and '1.2.3.4' are accepted. Additionally, both dashes ('-') and astrixes ('*') can be used to represent missing digits, both '1.2.3.-' and '1.2.3.\*' are accepted.

`cazy_webscraper` performs a direct EC number comparison. Therefore, supplying `cazy_webscraper` with the EC number EC1.2.3.- will only retrieve protein specifically annotated with EC1.2.3.-. `cazy_webscraper` will **not** retrieve proteins will all completed EC numbers under EC1.2.3.-, thus, `cazy_webscraper` will **not** retrieve data for proteins annotated with EC1.2.3.1, EC1.2.3.2, EC1.2.3.3, etc.

Example configuration files, and an empty configuraiton file template are located in the [`config_files`]() directory of this repo.

## Contributions

We welcome contributions and suggestions. You can raise issues at this repository, or fork the repository and submit pull requests, at the links below:

- [Issues](https://github.com/HobnobMancer/cazy_webscraper/issues)
- [Pull Requests](https://github.com/HobnobMancer/cazy_webscraper/pulls)

## License and copyright

MIT License

Copyright (c) 2020-2021 University of St Andrews  
Copyright (c) 2020-2021 University of Strathclyde  
Copyright (c) 2020-2021 UJames Hutton Institute  
