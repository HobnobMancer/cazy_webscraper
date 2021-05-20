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

The program additionally, includes an `expand` module, which can retrieve the protein sequences from [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) and protein structure files from the Research Collaboratory for Structural Bioinformatics (RCSB) Protein Data Bank [(PDB)](https://www.rcsb.org/).

The cazy_webscraper can be configured to scrape the entire database, selected Classes and/or selected CAZy families. An additional taxonomy filter can be applied to restrict the scraping of CAZymes to specific Kingdoms, genera, species and/or strains.

_For detailed documentation see the [full documentation](https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest)._

_An ER model can be found in the root of the `cazy_webscraper` GitHub repo, demonstrating the structure of the SQL database create by `cazy_webscraper`._

## Referencing

If you use `cazy_webscraper` in your work *please* do cite our work (including the provided DOI), as well as citing the specific version you use. This is not only helpful for us as developers to get out work out into the world, but it is also **essential for the reproducibility and integrity of scientific research**.  

**Citation:** Hobbs, Emma E. M.; Pritchard, Leighton; Chapman, Sean; Gloster, Tracey M. (2021): cazy_webscraper Microbiology Society Annual Conference 2021 poster. figshare. Poster. https://doi.org/10.6084/m9.figshare.14370860.v7 

## Our last developments before releasing version 1!

We're almost at the stage of releasing the first version of `cazy_webscraper`, but there are a few things we want to do first.

- Unit Tests: The coverage of the unit tests will be increased to cover at least 95% of the entire package

- The Expand module: This module is retrieving protein sequences from GenBank and PDB structures from PDB. The code is workable but inefficient. The code will be factorised out into multiple submodules, and include allow applying an additional EC number filter for selecting CAZyme to retrieve structural and sequence data for.

- Progress on the Expand module: Retrieving sequences for specific classes and families from a dictionary created using `cazy_webscraper` is complete, although no unit tests have been written yet.

- Documentation: The ReadTheDocs documentation has been updated to include all cmd-line flag. However, I am working on developing written and video tutorials to help users will less experience using the cmd-line to get started with this tool!

All this work is being done on the currently active `update_expand` branch.

## Future releases and features

To make `cazy_webscraper` accessible to a wider audience, I will (after the release of version 1) add on a GUI for `cazy_webscraper`. This will likely be slow in development but it will come in the next 6-9 months.

If there any additional features you would like to see added to `cazy_webscraper` or the documentation, please do raise an issue in the GitHub repo. All feedback is greatly appreciated and the aim is provide a tool users will find extremely useful!



## Requirements

POISx or Mac OS, or a Linux emulator  
Python version 3.8+  
Internet access  
The python libraries listed within `requirements.txt`  

## Installation


First clone the GitHub repository. This can be done at the command-line with the command:

```bash
git clone https://github.com/HobnobMancer/cazy_webscraper.git
```

Then change directory to the repository root, and use Python's setup tools to install:

```bash
cd cazy_webscraper
pip3 install -e .
```

## Getting started

For a summary for getting started, have a look at the [getting started poster](https://hobnobmancer.github.io/cazy_webscraper/getting_started_poster.pdf).

## Best practise

When performing a series of many, automated, repeated calls to a server it is best practise to do this during the period of the day when traffic is lowest. This typically includes over the weekend and overnight. Therefore, when scraping entire CAZy database, entire class(es), and/or multiple CAZy families it is advised perform these scrapes over night and/or over the weekend.

The webscraper can appear to run slowly but this is due to the limited access speed of the CAZy server. When a CAZy family is being parsed by the scraper, and protein records are being retrieved for the CAZy family a progress bar is produced in the terminal to provide an indicator that webscraper is working. **Expect an entire scrape of the CAZy database to take several hours.** Scraping individual families can be extremely rapid; scraping the entirity of the CAZy family GH1 (containing 42,647 proteins) takes approximately 42 minutes when using a 6-core AMD fx 63000 processor and 16GB RAM.

## Output


**Database:**

To facilitate thorough interrogation of data retrieve from CAZy, minimise storing duplicate and redundant data, data retrieved from CAZy is stored in a local SQL database. Every CAZyme scraped from CAZy has the following data:
- Protein name
- CAZy (sub)family
- GenBank accession(s)

Each CAZyme may or may not have the following data:
- EC number(s)
- UniProt acession(s)
- PDB accession(s)


**Primary and non-primary GenBank accessions:**  
Often multiple GenBank accession numbers are listed for a given CAZyme within CAZy. CAZy writes the 'best model' in bold, see [here for details](http://www.cazy.org/Help.html), this is interpretted by `cazy_webscraper` as the **primary GenBank accession**. Each unique CAZyme scraped from CAZy, is identified by it's **primary GenBank accession**. This enables associating all CAZy family annotations for a single CAZyme together. To be explicit, a single CAZyme record is for a CAZyme containing, for example, in GH32 and GH64, with both GH32 and GH64 annotations instead of creating a CAZyme record for the instance of the CAZyme in GH32 and another instance for the CAZyme in GH64. Another advantage of this is handling when there are duplicate records in CAZy. Identical entries are identified by two rows in the same HTML table containing identical data. If only a single GenBank accession is listed for a CAZyme (even if it is not written in bold), the lone GenBank accession is defined as the **primary GenBank accession**.

When multiple GenBank accessions are listed, the accession written in bold is listed as the **primary GenBank accession** and all other listed GenBank accessions are listed as **non-primary GenBank accessions**. For the instances when there are multiple GenBank accessions written in bold, only the first listed bold GenBank accession is listed as the **primary GenBank accession**. The remaining GenBank accessions written in bold are flagged up to the user and are listed as **non-primary GenBank accessions**. This method is to enable identifying each unique CAZyme from CAZy by a single unique **primary GenBank accession**.


**Primary and non-primary UniProt accessions**
IF multiple UniProt accessions are listed, all those written in bold are identified by CAZy and `cazy_webscraper` as the 'best' model. Consequently, all UniProt accessions listed in bold are defined as **primary UniProt accessions**, and all UniProt accessions not written in bold are listed as **non-primary UniProt accessions**. In cases when only a single UniProt accession, this lone accession is defined as the **primary UniProt accession** for the respective CAZyme.


**PDB accessions**
*All* PDB accessions stored in CAZy are retrieved, and are *not* differentitated between primary and non-primary accessions.

It is important to note that not all PDB accessions listed in CAZy are also represent in PDB. Therefore, not all PDB accessions retrieved from CAZy can be used to retrieve a protein structure file form PDB. For example, the PDB accession 'BBB[nd]' was previously listed in CAZy but did not represent a record in PDB.


## Configuration

The operation of the `cazy_webscraper` is configured by command-line arguments and a Yaml configuration file.

For the basic invoking of the `cazy_webscraper` use:  
`python3 cazy_webscraper.py`
This will result in `cazy_webscraper` scraping CAZymes from all CAZy families and species, and write the database to the memory.

### Command line arguments and operation

- `-c` `--config` Path to configuration yaml file
- `--classes` Define classes to be scraped. Separate classes with a single comma and **no spaces**
- `-d` `--database` Path to an existing database to add data to
- `-f` `--force` Force writing out in output directory that already exists
- `--families` Define families to be scraped. Separate families with a single comma (e.g. GH1,GH2) and **no spaces**
- `--genera` Genera to restrict the scrape to. Pass a string of generas separated by a comman and **no spaces**, e.g. `--genera "Trichoderma,Acidianus"`
- `--kingdoms` Taxonomy Kingdoms to be scraped. Separate Kingdoms with a single comma and **no spaces**
- `-h` `--help` Print option descriptions
- `-l` `--log` Write a log file, pass file to output log file
- `-n` `--nodelete` Do not delete content present in already existent output directory
- `-o` `--output` Specify output directory for CAZy dataframes
- `-r` `--retries` Number of times to attempt rescraping CAZy families if an error is raised. Default is 10.
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

To specify taxonomy filters (Kingdoms, genera, species and strains), list the respective filter under the appropriate keys (kingdoms, genera, species and strains) in the yaml file.


## Retrieving protein sequences and structure files

The `expand` module in `cazy_webscraper` manages the retrieval of CAZyme protein sequences from GenBank and protein structure files from PDB.
The `expand` module is invoked to act upon a local CAZy database, not invoked while scraping CAZy.

**Protein Sequences**
`get_genbank_sequences.py` retrieves protein sequences of CAZymes from GenBank.

To invoke the script use the following command:  
`python3 <path to get_genbank_sequences.py> <path to local CAZy database> <user email address>`  
A user email address is not stored by `cazy_webscraper`, it is a requirement of NCBI to use the Entrez server for automating the retrieval of protein sequences from GenBank. For more information about Entrez and NCBI please refer to their documentation...

Optional arguments are: _The format to pass arguments to the script is the same format use to pass the same argument to `cazy_webscraper`_
- `-c` `--config` Define path to a configuration yaml file to define CAZy classes and/or families to retrieve the protein sequences for
- `--classes` Define CAZy classes to retrieve protein sequences for
- `-e` `--epost` (Interget) Define the number of accessions to post to NCBI in each batch call. Default is 150. Adviced maximum is 200.
- `--families` Define CAZy families to retrieve protein sequences for
- `--genera` Genera of species to retrieve sequences for
- `--kingdoms` Taxonomy Kingdoms to retrieve sequences for
- `-l` `--log` Path to write out a log file to
- `-p` `--primary` Default False. If enabled only retrieve protein sequences for **primary** GenBank accessions
- `--species` Species to retrieve sequences for
- `--strains` Specific species strains to retrieve sequences for
- `-u` `--update` Default False. Enables overwriting sequences in the database only if the sequence in GenBank has been updated since the sequence was added to the local database. It still retrieves protein sequences for proteins that do not have a sequence stored in the local database.
- `-w` `--write` Write out the protein sequences to FASTA file, one protein per FASTA file
- `-v` `--verbose` Default False. Enables verbose logging, changing the logger level from `WARNING` to `INFO`.

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
- `-c` `--config` Define path to a configuration yaml file to define CAZy classes and/or families to retrieve the protein structures for
- `--classes` Define CAZy classes to retrieve protein structures for
- `--families` Define CAZy families to retrieve protein structures for
- `-f` `--force` Force writing to directory that already exists
- `--genera` Genera of species to retrieve structures for
- `--kingdoms` Taxonomy Kingdoms to retrieve structures for
- `-l` `--log` Path to write out a log file to
- `-n` `--nodelete` Default False. If enables it stops content in the existing output directory being deleted, else the contents in the existing output directory are deleted first before proceeding
- `-o` `--outdir` Default is the current working directory. Define path to output directory. If the output directory does not exist the program will create it.
- `-p` `--primary` Default False. If enabled only retrieve protein structures for **primary** GenBank accessions
- `--species` Species to retrieve structures for
- `--strains` Specific species strains to retrieve structures for
- `v` `--verbose` Default False. Enables verbose logging, changing the logger level from `WARNING` to `INFO`.


## Directories

Below is a directory plan of this repository, followed by a brief overview of each directories role , to facilitate navigation through the repository. This structure is the same structre found within the `cazy_webscraper` program.

### root directory

This contains the README and files necessary for installing `cazy_webscraper`.

### scraper

This directory houses all modules and Python scripts of the webscraper.

`cazy_webscraper.py` is the entry point of the program.

`crawler` is the module that coordinates crawling through the CAZy website and parsing retrived data.

`utilties` contains submodules for handling input/output files and directories (`file_io`), parsing configuration input (`parse_configuration`), building cmd-line argument parsers (`parsers`) and configuring the logger.

`expand` is the module for expanding the retrieved CAZyme data from CAZy beyond what can be retrieved directly from CAZy.
Specifically, the module enables retrieving protein sequences from GenBank and protein structure files from RSCB/PDB.

### tests

This directory contains all the unit tests for the webscraper. These tests are designed to be invoked from the root directory of the repository.

## Development and issues

If you have any issues or wish to see an additional feature added to `cazy_webscraper` please raise an _issue_ on GitHub.
