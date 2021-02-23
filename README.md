[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4300858.svg)](https://doi.org/10.5281/zenodo.4300858)
[![licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/cazy_webscraper/blob/master/LICENSE)
[![CircleCI](https://circleci.com/gh/HobnobMancer/cazy_webscraper.svg?style=shield)](https://circleci.com/gh/HobnobMancer/cazy_webscraper)
[![codecov](https://codecov.io/gh/HobnobMancer/cazy_webscraper/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/cazy_webscraper)
[![Documentation Status](https://readthedocs.org/projects/cazy-webscraper/badge/?version=latest)](https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest)
[![Python](https://img.shields.io/badge/Python-v3.8.---orange)](https://www.python.org/about/)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)

# `cazy_webscraper`

The `cazy_webscraper` is a Python3 package for the automated retrieval of protein data from the [CAZy](http://wwww.cazy.org/) database. This program is free to use under the MIT license when proper recognition is given.

The cazy_webscraper retrieves protein data from CAZy, producing a local SQL database which enables uses to throughly interrogate the data in a manner unachievable through the CAZy website.

The program additionally, includes an `expand` module, which can retrieve the protein sequences from GenBank and protein structure files from the Research Collaboratory for Structural Bioinformatics (RCSB) Protein Data Bank [(PDB)](https://www.rcsb.org/).

The cazy_webscraper can be configured to scrape the entire database, selected Classes and/or selected CAZy families. An additional taxonomy filter can be applied to restrict the scraping of CAZymes to specific genera, species and/or strains.

_For detailed documentation see the [full documentation](https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest)._

_An ER model can be found in the root of the `cazy_webscraper` GitHub repo, demonstrating the structure of the SQL database create by `cazy_webscraper`._

## Requirements

POISx or Mac OS, or a Linux emulator  
Python version 3.8+  
Internet access  
The python libraries listed within `requirements.txt`  

## Installation

### Quick and Easy 

`cazy_webscraper` is available in the [`bioconda`](https://bioconda.github.io/user/install.html) channel of [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/), and can be installed from the command-line with:

```bash
conda install cazy_webscraper
```

if you have the `bioconda` channel available, or

```bash
conda install -c bioconda cazy_webscraper
```

if you do not. This installs the full cazy_webscraper and all dependencies.

### From Source

First clone the GitHub repository. This can be done at the command-line with the command:

```bash
git clone https://github.com/HobnobMancer/cazy_webscraper.git
```

Then change directory to the repository root, and use Python's setup tools to install:

```bash
cd cazy_webscraper
python3 setup.py install
```

## Best practise

When performing a series of many, automated, repeated calls to a server it is best practise to do this during the period of the day when traffic is lowest. This typically includes over the weekend and overnight. Therefore, when scraping entire CAZy database, entire class(es), and/or multiple CAZy families it is advised perform these scrapes over night and/or over the weekend.

The webscraper can appear to run slowly but this is due to the limited access speed of the CAZy server. When a CAZy family is being parsed by the scraper, and protein records are being retrieved for the CAZy family a progress bar is produced in the terminal to provide an indicator that webscraper is working. **However, expect an entire scrape of the CAZy database to take several hours.**

## Output

**Dataframes:**

The basic function of the `cazy_webscraper` is to retrieve the protein data stored within and presented in the CAZy database. The data is written out to a dataframe with the same headings as found in CAZy, to reflect the way CAZy presented data in its webpages. The resulting dataframe includes the additional column “CAZy family”, which includes the CAZy family/subfamily under which the respective CAZy is catalogued. Therefore, the resulting dataframe of the webscraper contains the following columns:

- Protein_name
- CAZy_family
- EC#
- Source_organism
- GenBank
- UniProt
- PDB/3D

The scraped CAZymes can be written to a single dataframe, or separated out into different dataframes depending on how the data has been set to be split:

- Per family: a single dataframe is created per scraped CAZy family
- Per class: a singel datafarme is created per scraped CAZy family
- Not split: a single dataframe containing all data scrapend from CAZy is created

The dataframes can be written out to a specified directory or written to STDOUT to facilitate piping to a subsequent program. If the dataframes are written to the disk they are saved as .csv files.

**Primary and non-primary GenBank accessions:**  
Often multiple GenBank accession numbers are listed for a given CAZyme within CAZy. However, only the first listed accession number is hyperlinked to the GenBank database. Examination of the other listed synonyms (referred to as genbank synonyms in the webscraper) shows that these GenBank synonyms are the result of submission of identical protein sequences, splice site and protein isoforms. It has been interpreted that it is the first GenBank accession that is listed and hyperlinked to GenBank is the accession number of protein sequence which was used by CAZy to catalogue the CAZyme, and is referred to as the **primary GenBank accession**. The other GenBank accessions are **non-primary GenBank accessions**, and represent protein records in GenBank that have high sequence identity to the catalogued CAZyme.

CAZy retrieves protein records from GenBank in order to make its own (Lombard _et al_., 2013), therefore, in the local CAZy database created by `cazy_webscraper` each unique CAZyme is identified by its unique **primary GenBank accession**.


## Configuration

The operation of the `cazy_webscraper` is configured by command-line arguments and a Yaml configuration file.

For the basic invoking of the `cazy_webscraper` use:  
`python3 cazy_webscraper.py`


### Command line arguments and operation

- `-c` `--config` Path to configuration yaml file
- `--classes` Define classes to be scraped. Separate classes with a single comma and **no spaces**
- `-d` `--database` Path to an existing database to add data to
- `-f` `--force` Force writing out in output directory that already exists
- `-families` Define families to be scraped. Separate families with a single comma (e.g. GH1,GH2) and **no spaces**
- `genera` Genera to restrict the scrape to. Pass a string of generas separated by a comman and **no spaces**, e.g. `--genera "Trichoderma,Acidianus"`
- `-g` `--genbank` Entable retrieval of FASTA files, also pass user email
- `-h` `--help` Print option descriptions
- `l` `--log` Write a log file, pass file to output log file
- `n` `--nodelete` Do not delete content present in already existent output directory
- `o` `--output` Specify output directory for CAZy dataframes
- `r` `--retries` Number of times to attempt rescraping CAZy families if an error is raised. Default is 0, meaning no reattempted scrapes will be performed. Pass an interger
- `-s` `--subfamilies` Enable retrieval of subfamilies
- `--species` Species to restrict the scrape to, this will scrape CAZymes from all strains under the specified species. Pass as a string of species separated by a comma and **no spaces**
- `--strains` Specific strains of species to limit the scrape to.Pass as a string of species separated by a comma and **no spaces**
- `-v` `--verbose` Enable verbose logging


### Configuration files

For shareable documentation of the scraping CAZy a configuration file can by used to specify CAZy classes and families to be scraped.

Under the **classes** heading list any classes to be scrapped. For classes listed under 'classes', all proteins catalogued under that class will be retrieved, unless specific families have been listed under the respective classes heading in the configuration file. Then scraping only the specific classes takes precident and the entire class is not scraped. _If you believe this should be changed please raise an issue. It is invisioned that very few users would want to scrape an entire class an also scrape only specific families from that class simultanious._

A `cazy_dictionary.json` has been created and packaged within the `cazy_webscraper`. This allows users to use a variety of synonoms for the CAZy classes, for example both "GH" and "Glycoside-Hydrolases" are accepted as synonoms for "Glycoside Hydrolases (GHs)". This dictionary is packaged within the `scraper/file_io` directory. If you having issues with the scraper retrieving the list of CAZy classes that are written under 'classes' in the configuration file please check the dictionary first to see the full list of accepted synonoms. If you are comfortable modifying json files then feel free to add your own synonoms to the dictionary.

Under the each of the specific class names listed in the configuration file list the names of specific **families** to be scraped from that class. You do not have to list the class of the families to be scraped under 'classes' as well, this is handled by the web scraper.

Write the true name of the family not only it's number, for example **GH1** is excepted by **1** is not. Additionally, use the standard CAZy notation for subfamilies (**GH3_1**). If subfamilies are specified in the configuration file `--subfamilies` **must be enabled when invoking the cazy_webscraper.**

If the parent family, e.g GH3, is listed in the configuration file and `--subfamilies` is enabled, all proteins catalogued under GH3 and its subfamilies will be retrieved. This is to
save time having to write out all the subfamilies for a given CAZy family.

Each family must be listed on a separate line and the name surrounded by double or single quotation marks. For example:

```
Glycoside Hydrolases (GHs):
    - "GH1"
    - "GH2"
```

Please find more information on writing lists in Yaml files [here](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html).

A blank configuration file is packaged within `cazy_webscraper`, within the `scraper` directory, called `scraper_config.yaml`. This configuration file contains comments to assit filling in the file correctly. If a new configuration file is created it **must** be a Yaml file and it **must** use the same headings as used in the configuration file `scraper_config.yaml`.

## Retrieving protein sequences and structure files

The `expand` module in `cazy_webscraper` manages the retrieval of CAZyme protein sequences from GenBank and protein structure files from PDB.
The `expand` module is invoked to act upon a local CAZy database, not invoked while scraping CAZy.

**Protein Sequences**
`get_genbank_sequences.py` retrieves protein sequences of CAZymes from GenBank.

To invoke the script use the following command:  
`python3 <path to get_genbank_sequences.py> <path to local CAZy database> <user email address>`  
A user email address is not stored by `cazy_webscraper`, it is a requirement of NCBI to use the Entrez server for automating the retrieval of protein sequences from GenBank. For more information about Entrez and NCBI please refer to their documentation...

Optional arguments are: _The format to pass arguments to the script is the same format use to pass the same argument to `cazy_webscraper`_
- `c` `--config` Define path to a configuration yaml file to define CAZy classes and/or families to retrieve the protein sequences for
- `--classes` Define CAZy classes to retrieve protein sequences for
- `-e` `--epost` (Interget) Define the number of accessions to post to NCBI in each batch call. Default is 150. Adviced maximum is 200.
- `--families` Define CAZy families to retrieve protein sequences for
- `--genera` Genera of species to retrieve sequences for
- `l` `--log` Path to write out a log file to
- `p` `--primary` Default False. If enabled only retrieve protein sequences for **primary** GenBank accessions
- `--species` Species to retrieve sequences for
- `--strains` Specific species strains to retrieve sequences for
- `u` `--update` Default False. Enables overwriting sequences in the database only if the sequence in GenBank has been updated since the sequence was added to the local database. It still retrieves protein sequences for proteins that do not have a sequence stored in the local database.
- `w` `--write` Write out the protein sequences to FASTA file, one protein per FASTA file
- `v` `--verbose` Default False. Enables verbose logging, changing the logger level from `WARNING` to `INFO`.

**Protein Structures**
`get_pdb_structures.py` retrieves protein structure files from PDB. The downloading of the structure files is handled by the `BioPython` `PDB` module, further information on the module can be found [here](https://biopython.org/docs/1.75/api/Bio.PDB.html).

to invoke the script, use the following command:  
`python3 <path to get_pdb_structures.py> <path to local CAZy database> <file format of downloaded files>`  
The file formats and descriptions are taken straight from the `BioPython` documentation [found here](https://biopython.org/docs/1.75/api/Bio.PDB.PDBList.html?highlight=pdblist#Bio.PDB.PDBList.PDBList.retrieve_pdb_file).
- "mmCif" (default, PDBx/mmCif file),
- "pdb" (format PDB),
- "xml" (PDBML/XML format),
- "mmtf" (highly compressed),
- "bundle" (PDB formatted archive for large structure}

Optional arguments: _The format to pass arguments to the script is the same format use to pass the same argument to `cazy_webscraper`_
- `c` `--config` Define path to a configuration yaml file to define CAZy classes and/or families to retrieve the protein sequences for
- `--classes` Define CAZy classes to retrieve protein structures for
- `--families` Define CAZy families to retrieve protein structures for
- `--genera` Genera of species to retrieve sequences for
- `f` `--force` Force writing to directory that already exists
- `l` `--log` Path to write out a log file to
- `n` `--nodelete` Default False. If enables it stops content in the existing output directory being deleted, else the contents in the existing output directory are deleted first before proceeding
- `o` `--outdir` Default is the current working directory. Define path to output directory. If the output directory does not exist the program will create it.
- `p` `--primary` Default False. If enabled only retrieve protein structures for **primary** GenBank accessions
- `--species` Species to retrieve sequences for
- `--strains` Specific species strains to retrieve sequences for
- `v` `--verbose` Default False. Enables verbose logging, changing the logger level from `WARNING` to `INFO`.


## Directories

Below is a directory plan of this repository, followed by a brief overview of each directories role , to facilitate navigation through the repository. This structure is the same structre found within the `cazy_webscraper` program.

### root directory

This contains the README and files necessary for installing `cazy_webscraper`.

### scraper

This directory houses all modules and Python scripts of the webscraper.

`cazy_webscraper.py` is the entry point of the program.

`crawler` is the module that coordinates crawling through the CAZy website and parsing retrived data.

`utilties` contains functions for building cmd-line args parsers and loggers.

`file_io` contains all functions related to handingling input/ouput files and directories, this includes creating the output directory, parsing configuraiton files and writing the output dataframes.

`parse` contains functions related to parsing the protein data retrieved from CAZy into Pandas dataframes.

`expand` is the module for expanding the retrieved CAZyme data from CAZy beyond what can be retrieved directly from CAZy.
Specifically, the module enables retrieving protein sequences from GenBank and protein structure files from RSCB/PDB.

### tests

This directory contains all the unit tests for the webscraper. These tests are designed to be invoked from the root directory of the repository.

## Development and issues

If you have any issues or wish to see an additional feature added to `cazy_webscraper` please raise an _issue_ on GitHub.
