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

# Beta Version 2
This release of `cazy_webscraper` is a beta release. The main component of `cazy_webscraper` (for the downloading of data from CAZy) is stable. Expanding the data set to incorproate data from UniProt, GenBank and PDB still require some work to workout all the remaining bugs.

The documentation is being actively updated to match the new `cazy_webscraper` version 2.

**The `bioconda` installation method is not currently support, but we aim to get this fixed soon**. For now please install via pypi or from source.

## New features in version 2:
- **Faster scraping:** The entirtity of CAZy can be scraped in 15 minutes
- **Retrieval of UniProt data:** UniProt accessions, EC numbers, protein sequences, and PDB accessions can be retrieved from UniProt and added to the local CAZyme database
- **Addition of an API:** As well as retrieving data from the local CAZyme database via an SQL interface, `cazy_webscraper` can retrieved user-specfied data (e.g. the GenBank protein accession and EC number annotations) for proteins matching user-specified critieria. The extracted data can be written to a `JSON` and/or `CSV` file, to facilitate inclusion in downstream analyses.
- **Caching:** Data downloaded from CAZy is not only parsed and written to a local CAZyme database. The raw data files are written to cache. Data can be scraped directly from a cache (ideal if CAZy updates during the retrieval of multiple datasets from the CAZy database).

## Future work for version 2:
- Fix any remaining bugs we can find (if you find a bug, please report it!)
- Update the unit tests to work with the new `cazy_webscraper` architecture
- Update the documentation
- Create video tutorials

## cazy_webscraper

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
- [Extracting protein sequences from the local CAZyme database and building a BLAST database](#extracting-protein-sequences-from-the-local-cazyme-database-and-building-a-blast-database)
    - [Configuring extracting sequences from a local CAZyme db](#configuring-extracting-sequences-from-a-local-cazyme-db)
- [Retrieving protein structure files from PDB](#retrieving-protein-structure-files-from-pdb)
    - [Configuring PDB protein structure file retrieval](#configuring-pdb-protein-structure-file-retrieval)
- [The `cazy_webscraper` API or Interrogating the local CAZyme database](#the_cazy_webscraper_api_or_interrogating_the_local_cazyme_database)
- [Configuring `cazy_webscraper` using a YAML file](#configuring-using-a-yaml-file)
- [CAZy coverage of GenBank](#cazy-coverage-of-genbank)
    - [Configure calculating CAZy coverage of GenBank](#configure-calculating-cazy-coverage-of-genbank)
- [Contributions](#contributions)
- [License and copyright](#license-and-copyright)
<!-- /TOC -->

## Best practice

When performing a series of many automated calls to a server it is best to do this when traffic is lowest, such as at weekends or overnight *at the server*.

## Documentation

Please see the [full documentation at ReadTheDocs](https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest).

### Installation

`cazy_webscraper` can be installed via `conda` or `pip`:

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

`--cazy_data` - Path to a txt file downloaded from CAZy containing a CAZy db data dump.

`--cazy_synonyms` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

`--classes` - list of classes from which all families are to be scrape.

`--config`, `-c` - Path to a configuration YAML file. Default: None.

`--citation`, `-C` - Print the `cazy_webscraper` citation. When called, the program terminates after printng the citation and CAZy is **not** scraped.

`--database`, `-d` - Path to an **existings** local CAZyme database to add newly scraped too. Default: None.


_Do not use `--db_output` **and** `--database` at the same time._

_If `--db_output` **and** `--database` are **not** called, `cazy_webscraper` write out a local CAZyme database to the cwd with the standardised name `cazy_webscraper_<date>_<time>.db`_

`--delete_old_relationships` - Detele old CAZy family annotations of GenBank accessions. These are CAZy family annotations of a given GenBank accession are in the local database but the accession is not longer associated with those CAZy families, so delete old accession-family relationships.

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

`database` - \[REQUIRED\] Path to a local CAZyme database to add UniProt data to.

`--bioservices_batch_size` - Change the query batch size submitted via [`bioservices`](https://bioservices.readthedocs.io/en/master/) to UniProt to retrieve protein data. Default is 150. `bioservices` recommands queries not larger than 200 objects.

`--cache_dir` - Path to cache dir to be used instead of default cache dir path.

`--cazy_synonyms` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

`--classes` - List of classes to retrieve UniProt data for.

`--config`, `-c` - Path to a configuration YAML file. Default: None.

`--ec`, `-e` - Enable retrieval of EC number annotations from UniProt

`--ec_filter` - Limist retrieval of protein data to proteins annotated with a provided list of EC numbers. Separate the EC numbers bu single commas without spaces. Recommend to wrap the entire str in quotation marks, for example:
```bash
cw_get_uniprot_data my_cazyme_db/cazyme_db.db --ec_filter 'EC1.2.3.4,EC2.3.1.-'
```

`--families` - List of CAZy (sub)families to scrape.

`--genera` - List of genera to restrict the scrape to. Default: None, filter not applied to scrape.

`--log`, `-l` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

`--nodelete_cache` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

`--nodelete_log` - When called, content in the existing log dir will **not** be deleted. Default: False (existing content is deleted).

`--retries`, `-r` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

`--sequence`, `-s` - Retrieve protein amino acid sequences from UniProt

`--sql_echo` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

`--species` - List of species written as Genus Species) to restrict the scraping of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

`--strains` - List of specific species strains to restrict the scraping of CAZymes to.

`--timeout`, `-t` - Connection timout limit (seconds). Default: 45.

`--uniprot_batch_size` - Size of an individual batch query submitted to the [UniProt REST API](https://www.uniprot.org/help/programmatic_access) to retrieve the UniProt accessions of proteins identified by the GenBank accession. Default is 150. The UniProt API documentation recommands batch sizes of less than 20,000 but batch sizes of 1,000 often result in HTTP 400 errors. It is recommend to keep batch sizes less than 1,000, and ideally less than 200.

`--seq_update` - If a newer version of the protein sequence is available, overwrite the existing sequence for the protein in the database. Default is false, the protein sequence is **not** overwritten and updated.

`--verbose`, `-v` - Enable verbose logging. This does not set the SQLite engine `echo` parameter to True. Default: False.

## Retrieveing protein seqences from GenBank

Protein amino acid sequences can be retrieved for proteins in a local CAZyme database using `cazy_webscraper`. Protein sequences can be retrieved for a specific subset of proteins, identified through the use of CAZy class, CAZy family, taxonomy (kingdom, genus, species and strain) filters, and EC number filters. The retrieved protein sequences are written to the local CAZyme database.

_Extracting protein sequences from the local CAZyme database and writing them to a BLAST database and/or FASTA file(s) is covered in the next section._

To retrieve all GenBank protein seuqneces for all proteins in a local CAZyme datbase, using the following command:
```bash
cw_get_genbank_seqs <path_to_local_CAZyme_db>
```

`cazy_webscraper` produces to cache files, which are written to the cache dir:
1. `no_seq_retrieved.txt` which lists the GenBank accessions for which no sequence could be retrieved from GenBank
2. `seq_retrieved.txt` which list GenBank accessiosn for which a sequence was retrieved from GenBank

### Configuring GenBank protein sequence retrieval

Below are listed the command-line flags for configuring the retrieval of protein sequences from GenBank.

`database` - \[REQUIRED\] Path to a local CAZyme database to add UniProt data to.

`email` - \[REQUIRED\] User email address, required by NCBI Entrez.

`--cache_dir` - Path to cache dir to be used instead of default cache dir path.

`--cazy_synonyms` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

`--config`, `-c` - Path to a configuration YAML file. Default: None.

`--classes` - List of classes from which all families are to be scrape.

`--ec_filter` - Limist retrieval of protein data to proteins annotated with a provided list of EC numbers. Separate the EC numbers bu single commas without spaces. Recommend to wrap the entire str in quotation marks, for example:
```bash
cw_get_uniprot_data my_cazyme_db/cazyme_db.db --ec_filter 'EC1.2.3.4,EC2.3.1.-'
```

`--entrez_batch_size` - Change the query batch size submitted via [`Entrez`]() to retrieve protein sequences from GenBank data. Default is 150. `Entrez` recommands queries not larger than XXX objects in length.

`--families` - List of CAZy (sub)families to scrape.#

`--kingdoms` - List of taxonomy kingdoms to retrieve UniProt data for.

`--genera` - List of genera to restrict the scrape to. Default: None, filter not applied to scrape.

`--log`, `-l` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

`--nodelete_cache` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

`--retries`, `-r` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

`--seq_update` - If a newer version of the protein sequence is available, overwrite the existing sequence for the protein in the database. Default is false, the protein sequence is **not** overwritten and updated.

`--sql_echo` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

`--species` - List of species written as Genus Species) to restrict the scraping of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

`--strains` - List of specific species strains to restrict the scraping of CAZymes to.

`--verbose`, `-v` - Enable verbose logging. This does not set the SQLite engine `echo` parameter to True. Default: False.

## Extracting protein sequences from the local CAZyme database and building a BLAST database

Protein sequences from GenBank and UniProt that are stored in the local CAZyme database can be extracted using `cazy_webscraper`, and written to any combination of:
- 1 FASTA file per unique protein
- A single FASTA file containing all extracted seqences
- A BLAST database

**FASTA file format:** Protein sequences extracted from a local CAZyme database are written out with the GenBank/UniProt accession as the protein ID, and the name of the source database ('GenBank' or 'UniProt') as the description.

To extract all protein seqeunces from the local CAZyme database using the following command structure:
```bash
cw_extract_sequences <path_to_local_CAZyme_db> --genbank --uniprot
```

### Configuring extracting sequences from a local CAZyme db

Below are listed the command-line flags for configuring the extraction of protein sequences from the local CAZyme db.

`database` - \[REQUIRED\] Path to a local CAZyme database to add UniProt data to.

`-g`, `--genbank` - Extract protein sequences retrieved from GenBank

`-u`, `--uniprot` - Extract protein sequences retrieved from UniProt.

_Note: At least one of `--genbank` and `--uniprot` must be called, otherwise not protein sequences will be extracted._

`-b`, `--blastdb` - Create BLAST database of extracted protein sequences. Provide the path to the directory to store the BLAST database in.

`--fasta_dir` - Write out each extracted sequence to a separate FASTA file in the provided dir. Provide a path to a directory to write out the FASTA files.

`--fasta_file` - Write out all extracted sequences to a single FASTA file. Provide a path to write out the FASTA file.

_Note: at least one of `--blastdb`, `--fasta_dir`, and `--fasta_file` must be called to inform `cazy_webscraper` where to write the output to. If none are called sequences will be extracted._

`--cache_dir` - Path to cache dir to be used instead of default cache dir path.

`--cazy_synonyms` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

`--config`, `-c` - Path to a configuration YAML file. Default: None.

`--classes` - List of classes from which all families are to be scrape.

`--ec_filter` - Limist retrieval of protein data to proteins annotated with a provided list of EC numbers. Separate the EC numbers bu single commas without spaces. Recommend to wrap the entire str in quotation marks, for example:
```bash
cw_get_uniprot_data my_cazyme_db/cazyme_db.db --ec_filter 'EC1.2.3.4,EC2.3.1.-'
```

`--families` - List of CAZy (sub)families to scrape.#

`--kingdoms` - List of taxonomy kingdoms to retrieve UniProt data for.

`--genera` - List of genera to restrict the scrape to. Default: None, filter not applied to scrape.

`--log`, `-l` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

`--nodelete` - When called, content in the existing output dir will **not** be deleted. Default: False (existing content is deleted).

`--nodelete_cache` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

`--sql_echo` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

`--species` - List of species written as Genus Species) to restrict the scraping of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

`--strains` - List of specific species strains to restrict the scraping of CAZymes to.

`--verbose`, `-v` - Enable verbose logging. This does not set the SQLite engine `echo` parameter to True. Default: False.


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

`database` - \[REQUIRED\] Path to a local CAZyme database to add UniProt data to.

`pdb` \[REQUIRED\] The file types to be retrieved from PDB. The following file types are supported:  
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

`--cache_dir` - Path to cache dir to be used instead of default cache dir path.

`--cazy_synonyms` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

`--config`, `-c` - Path to a configuration YAML file. Default: None.

`--classes` - List of classes from which all families are to be scrape.

`--ec_filter` - Limist retrieval of protein data to proteins annotated with a provided list of EC numbers. Separate the EC numbers bu single commas without spaces. Recommend to wrap the entire str in quotation marks, for example:
```bash
cw_get_uniprot_data my_cazyme_db/cazyme_db.db --ec_filter 'EC1.2.3.4,EC2.3.1.-'
```

`--families` - List of CAZy (sub)families to scrape.

`--genera` - List of genera to restrict the scrape to. Default: None, filter not applied to scrape.

`--log`, `-l` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

`--nodelete` - When called, content in the existing output dir will **not** be deleted. Default: False (existing content is deleted).

`--nodelete_cache` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

`--outdir`, `-o` - Output directory to write out downloaded protein structure files to. Default is to write out the downloaded structure files to the current working directory.

`--retries`, `-r` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

`--sql_echo` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

`--species` - List of species written as Genus Species) to restrict the scraping of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

`--strains` - List of specific species strains to restrict the scraping of CAZymes to.

`--timeout`, `-t` - Connection timout limit (seconds). Default: 45.

`--verbose`, `-v` - Enable verbose logging. This does not set the SQLite engine `echo` parameter to True. Default: False.


## The `cazy_webscraper` API or Interrogating the local CAZyme database

The SQLite3 database compiled by `cazy_webscraper` can be interrogated in the native interface (i.e. queries written in SQL can be used to interrogate the database). This can be achieved via the command-line or via an SQL database browser (such as [DB Browser for SQLite](https://sqlitebrowser.org/)).

`cazy_webscraper` also provides its own API (Application Programming Interface) for interrogating the local CAZyme database: `cw_query_database`. The API faciliates the intergration of the dataset in the local CAZyme database into downstream bioinformatic pipelines, _and_ provides a method of interrograting the dataset for those who do not use SQL.

`cw_query_database` is the command that can be used to interrogate the dataset in the local CAZyme database, and extract protein data of interest for the proteins matching the user's cirteria of 
interest.

By default `cw_query_database` retrieves only the GenBank accessions of proteins matching the user's criteria of interest. If not criteria of interest are provided, all GenBank accessions are retrieved. Optional flags can be applied to retrieve additional data about CAZymes that match the user's criteria of interest.

`cw_query_database` currently supports writing the output in two file formats:
- `json`, by default this is written to the current working directory and with name `<database_name>_<retrieved_data>_YYYY-MM-DD_HH-MM-SS.json`
- `csv`, by default this is written to the current working directory and with name `<database_name>_<retrieved_data>_YYYY-MM-DD_HH-MM-SS.csv`

`cw_query_database` takes two positional arguments:
1. The path to the local CAZyme database
2. The file formats of the output files, presented as a list with each file type separated by as single comma (e.g. `json,csv`). This is **not** case sensitive and the order does **not** matter.
For example, to retrieve all GenBank accessions for all proteins in the local CAZyme database, and 
write them to a `json` file, the following command could be used for a database called `cazy.db`:
```bash
cw_query_database cazy.db json
```

### Configuring interrogating the local CAZyme database

Below are listed the command-line flags for configuring the interrogation of the local CAZyme database.

`database` - \[REQUIRED\] Path to a local CAZyme database to add UniProt data to.

`file_types` - \[REQUIRED\] file types to write the interrogation output to. Accepted file types are
`JSON` and `CSV`. These are *not* case sensitive, and the order does not matter.

`--cazy_synonyms` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

`--config`, `-c` - Path to a configuration YAML file. Default: None.

`--class` - Include a 'Class' column in the output `csv` file, listing the CAZy class of all retrieved CAZymes

`--classes` - List of classes from which all families are to be retrieval.

`--ec` - Included EC number annotations in the output

`--ec_filter` - Limist retrieval of protein data to proteins annotated with a provided list of EC numbers. Separate the EC numbers bu single commas without spaces. Recommend to wrap the entire str in quotation marks, for example:
```bash
cw_get_uniprot_data my_cazyme_db/cazyme_db.db --ec_filter 'EC1.2.3.4,EC2.3.1.-'
```

`--family` - Include CAZy family annotations in the output.

`--families` - List of CAZy (sub)families to retrieve CAZymes from. This includes families and SUBfamilies.

`--genus` - Include a 'Genus' column in the output `csv` file, listing the genus for each retrieved CAZymes

`--genera` - List of genera to restrict the retrieval to. Default: None, filter not applied to scrape.

`--log`, `-l` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

`--nodelete` - When called, content in the existing output dir will **not** be deleted. Default: False (existing content is deleted).

`--organism` - Include the scientific (organism) name of the source organis

`--output`, `-o` - Output path to write the compiled `csv` file. Default is to write out the `csv` file to the current working directory.

`--seq_genbank` - Add the protein sequences retrieved from GenBank to the query output

`--seq_uniprot` - Add the protein sequences retrieved from GenBank to the query output

`--sql_echo` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

`--species` - List of species written as (Genus Species) to restrict the retrieval of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

`--strains` - List of specific species strains to restrict the retrieval of CAZymes to.

`--subfamily` - Include the CAZy subfamily annotations for proteins matching the criteria of interest.

`--uniprot` - Include the UniProt protein accession and protein name in the output.

`--verbose`, `-v` - Enable verbose logging. This does not set the SQLite engine `echo` parameter to True. Default: False.


## Configuring `cazy_webscraper` using a YAML file

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


## CAZy coverage of GenBank

The number of genomes represented by the local CAZyme database per taxonomy Kingdom can be compared to the number 
GenBank genomic assemblies in NCBI. This done using the command `cw_cazy_genbank_coverage`, which requires two 
positional arguments:
1. Path to the local CAZyme database
2. The users email address
For example:
```bash
cw_cazy_genbank_coverage cazymes.db my_email@domain.com
```

This produces several output files (where time stamp is the date and time the command was invoked):
1. `cazy_genomic_accessions_<time_stamp>.json`, a JSON file contained a multi-layer dictionary: {kingdom: {genus: {species: {genomic_accession: {proteins: {protein_accessions}, count=int(# of proteins)}}}}}
2. `genomic_accessions_<time_stamp>.csv`, containing the columns 'Kingdom','Genus','Species','Genomic_accession','#ofProteins', where the number of proteins represents the number of proteins from the genome which are catalogued in CAZy
3. `protein_genomic_accessions_<time_stamp>.csv`, containing the 'Kingdom', 'Genus', 'Species', 'Genomic_accession', 'Protein_accession'. A unique protein accession is listed on each row, and lists which protein accesison is derived from each genomic assembly.
4. `cazy_genbank_genome_coverage_<time_stamp>.csv`, containing the columns 'Kingdom', 'NCBI_genomes', 'CAZy_genomes', 'Coverage_percent'. The dataframe lists the number of genomes catalogued in GenBank (NCBI) and the local CAZyme database per taxonomy Kingdom.
5. `gbk_cazy_genomes_plot_<time_stamp>.png`, a plot of a stacked bar chart with the number of genomes in GenBank (NCBI) and CAZy (the local CAZyme database) per taxonomy Kingdom.

### Configure calculating CAZy coverage of GenBank

Optional cmd-line arguments for `cw_cazy_genbank_coverage` are listed below:

`--batch_size` - The number of accessions posted to NCBI per epost, advice to be max 200. Default=150.

``--force``, ``-f`` - Force writing in existing output directory.

``--force_cache`` - Force writing in existing cache directory.

`--log`, `-l` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

`--nodelete` - When called, content in the existing output dir will **not** be deleted. Default: False (existing content is deleted).

`--nodelete_cache` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

`--output_dir`, `-o` - Path to output directory.

`--retries` - Number of times to retry connection to NCBI Entrez if connection fails.

`--sql_echo` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

`--verbose`, `-v` - Enable verbose logging. This does not set the SQLite engine `echo` parameter to True. Default: False.

## Contributions

We welcome contributions and suggestions. You can raise issues at this repository, or fork the repository and submit pull requests, at the links below:

- [Issues](https://github.com/HobnobMancer/cazy_webscraper/issues)
- [Pull Requests](https://github.com/HobnobMancer/cazy_webscraper/pulls)

## License and copyright

MIT License

Copyright (c) 2020-2021 University of St Andrews  
Copyright (c) 2020-2021 University of Strathclyde  
Copyright (c) 2020-2021 James Hutton Institute  
