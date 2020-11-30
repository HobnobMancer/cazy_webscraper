![bioconda](assets/bioconda-badge-wide.png)

[![licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/cazy_webscraper/blob/master/LICENSE)
[![CircleCI](https://circleci.com/gh/HobnobMancer/cazy_webscraper.svg?style=shield)](https://circleci.com/gh/HobnobMancer/cazy_webscraper)
[![codecov](https://codecov.io/gh/HobnobMancer/cazy_webscraper/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/cazy_webscraper)
[![Documentation Status](https://readthedocs.org/projects/cazy-webscraper/badge/?version=latest)](https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest)
[![Python](https://img.shields.io/badge/Python-v3.8.---orange)](https://www.python.org/about/)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)

# cazy_webscrapper

The `cazy_webscraper` is a Python3 package for the automated retrieval of protein data from the [CAZy](http://wwww.cazy.org/) database. This program is free to use under the MIT license when proper recognition is given.

The cazy_webscraper retrieves protein data from CAZy, writing out the data to a dataframe in the same manner as the data is presented in the CAZy website. If enabled, the webscraper will also retrieve the protein sequence of the scraped CAZymes from GenBank, writing out the sequences in FASTA format. Additionally, if enabled the webscraper will retrieve all protein structures for each scraped CAZyme from the Research Collaboratory for Structural Bioinformatics (RCSB) Protein Data Bank [(PDB)](https://www.rcsb.org/).

The cazy_webscraper can be configured to scrape the entire database, selected Classes and/or selected CAZy families. Additionally, the retrieved data can be separated/split by CAZy class, CAZy family or not at all, thus gathering all data into a single dataframe.

_For detailed documentation see the [full documentation](https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest)._

## Requirements

POISx or Mac OS, or a Linux emulator  
Python version 3.8+  
Internet access  
The python libraries listed within `requirements.txt`  

## Installation

**Quick and Easy:** The easiest method for installing the cazy_webscraper is to use Conda, using the following command at the command-line in the terminal: `conda asdasd cazy_webscraper` This method installs the full cazy_webscraper and all dependencies.

If Conda is not installed, please see the Conda website for installation [instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

**Alternative, script access:** Alternatively, for easier access to Python scripts that make up the cazy_webscraper first clone the webscrapers GitHub repository, then install the webscraper using pips: `git clone https://github.com/HobnobMancer/cazy_webscraper pip3 install -e <path_to_dir_containing_the_setup.py_file>`

Then install remaining requirements: `pip3 install -r <path_to_requirements.txt_file>`

In both commands to do not forget the additional pips options (-e and -r)!

To write the webscraper repository to a specific directory use the following command: `git clone https://github.com/HobnobMancer/cazy_webscraper> <path_to_dir>`  The directory to which the path points will form the root of the local copy of the repository.

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

**GenBank synonyms:**  
Often multiple GenBank accession numbers are listed for a given CAZyme within CAZy. However, only the first listed accession number is hyperlinked to the GenBank database. Examination of the other listed synonyms (referred to as genbank synonyms in the webscraper) shows that these GenBank synonyms are the result of submission of identical protein sequences, splice site and protein isoforms. It has been interpreted that it is the first GenBank accession that is listed and hyperlinked to GenBank is the accession number of protein sequence which was used by CAZy to catalogue the CAZyme and the GenBank synonyms were identified and listed by having extremely high sequence identity to the catalogued CAZyme.

Therefore, the webscraper writes only the first GenBank accession listed for each CAZyme in the resulting dataframe. The remaining GenBank synonyms are written out to a JSON file, keyed by the first GenBank accession given for each CAZyme, and valued by a list of GenBank synonyms. If no GenBank synonyms are retrieved for a CAZyme then the CAZyme’s GenBank accession is not written out to the JSON GenBank synonyms file.

The GenBank synonyms file is written out to the same directory as specificed for the dataframes. Additionally, the data is split as is specified for the dataframes.

**Protein sequences:**  

If enabled, the protein sequence of the scraped CAZymes are retrieved from GenBank are retrieved in the FASTA format and can be written to STDOUT to facilitate piping to a subsequent program or written out to disk, within a specified directory.

**Protein structures:**  

If enabled, the protein structures will be written out to the disk, to a specified directory. The protein structures cannot be written to STDOUT due to using the BioPython module PDB, which currently does not facilitate writing out the protein structures to STDOUT. The format of the the structure file is specified at the command line.


## Configuration

The operation of the `cazy_webscraper` is configured by command-line arguments and a Yaml configuration file.

For the basic invoking of the `cazy_webscraper` use:  
`python3 cazy_webscraper.py`

### Command line arguments and operation

- `-c` `--config` Path to configuration yaml file.
- `-d` `--data_split` Split data by CAZy class, CAZy family or not at all
- `-f` `--force` Force writing out in output directory that already exists
- `-g` `--genbank` Entable retrieval of FASTA files, also pass user email
- `-genbank_output` Specify output directory for FASTA files
- `-h` `--help` Print option help descriptions
- `l` `--log` Write a log file, pass file to output log file
- `n` `--nodelete` Do not delete content present in already existent output directory
- `o` `--output` Specify output directory for CAZy dataframes
- `p` `--pdb` Enable retrieval of protein structures from PDB, also pass desired file format of structure files
- `-pdb_output` Specify output directory for protein structure files
- `-s` `--subfamilies` Enable retrieval of subfamilies
- `-v` `--verbose` Enable verbose logging


### Configuration files

The configuration file is for specifying specific CAZy classes and families to be scraped.

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

## Directories

Below is a directory plan of this repository, followed by a brief overview of each directories role , to facilitate navigation through the repository. This structure is the same structre found within the `cazy_webscraper` program.

### root directory

This contains the README and files necessary for installing `cazy_webscraper`.

### scraper

This directory houses all modules and Python scripts of the webscraper.

`cazy_webscraper.py` is the entry point of the program.

`utilties` contains functions for building cmd-line args parsers and loggers.

`file_io` contains all functions related to handingling input/ouput files and directories, this includes creating the output directory, parsing configuraiton files and writing the output dataframes.

`parse` contains functions related to parsing the protein data retrieved from CAZy into Pandas dataframes.

### tests

This directory contains all the unit tests for the webscraper. These tests are designed to be invoked from the root directory of the repository.

## Development and issues

This webscraper is still in early development and is not fully opperational yet. Furthermore, additional features may be added later, such as ability to configure the webscraper to only scrape the data for a single CAZy class instead of scraping the entire database.

For more information on planning/development please see the Wiki.

If you have additional features wish to be added or any other issues please raise and issue in this repoistory.
