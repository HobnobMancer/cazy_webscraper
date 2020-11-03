# cazy_webscrapper
The `cazy_webscraper` is a Python3 package for the automated retrieval of protein data from the [CAZy](http://wwww.cazy.org/) database. This program is free to use under the MIT license when proper recognition is given.

The CAZy database is one of the largest databases for cataloging Carbohydrate Active enZymes, and is one of the most frequently visited protein data bases for researchers in this field. However, CAZy provides not method of automated retrieval of data, especially for large dataset requests. The `cazy_webscraper` provides a way for the automated retrieval of large datasets from the CAZy database. The retrieved protein data is written out into dataframes, and protein sequences written out to FASTA files.

The `cazy_webscraper` is configurable to scrape specific CAZy classes and/or families, and determine how the data split when written out. The data can be split:

- Per family: a single dataframe is created per scraped CAZy family
- Per class: a singel datafarme is created per scraped CAZy family
- Not split: a single dataframe containing all data scrapend from CAZy is created

## Requirements

POISx or Mac OS, or a Linux emulator
Python version 3.8+
Internet access
The python libraries listed within `requirements.txt`, packaged within the `cazy_webscraper`

## Installation

The easiest method for install the `cazy_webscraper` is to use pips.  
`pip3 install -e <path_to_dir_containing_the_setup.py_file>`  

Then install remaining requirements:  
`pip3 install -r <path_to_requirements.txt_file>`  

In both commands to not forget the additional pips options (`-e` and `-r`)!

## Configuration

The operation of the `cazy_webscraper` is configured by command-line arguments and a Yaml configuration file.

For the basic invoking of the `cazy_webscraper` use:  
`python3 cazy_webscraper.py`

### Command line arguments and operation

- `-c`, `--config`: Path to the configuration file. Default: None, scrapes entire CAZy database
- `-d`, `--data_split`: Choices: None, class, family. Default: None, not to split data when written out
- `-f`, `--force`: (True/False). Force over writing in output directory if specified output directory already exists. Default: False, does not over write in already exising output dictory. 
- `-l`, `--log`: Path to write out a logger file. Default: None. Logger messages will be written out to the terminal and out to the specified file.
- `-n`, `--nodelete`: (True or False). Do not delete content in already exisiting output directory. Default: False, will delete content in already existing output directory. If set to true then the content in the output directory will not be deleted first before writing out output from the scrape.
- `-o`, `--output`: Path to output DIRECTORY for all output to be written to. Default: STDOUT. The dataframe names are pre-formated by the scraper so only pass the path to the directory into which the output data is to be written. If the directory does not already exist the `cazy_webscraper` will create the output dataframe.
- `-v`, `--verbose`: (True or False) Change the logger level from Warning to Info, resulting in logging of the scrapers progress. Default: false.

### Configuration files

The configuration file is for specifying specific CAZy classes and families to be scraped.

Under the **classes** heading list any classes to be scrapped. For classes listed under 'classes', all proteins catalgoued under that class will be retrieved, unless specific families have been listed under the respective classes heading in the configuration file. Then scraping only the specific classes takes precident and the entire class is not scraped. _If you believe this should be changed please raise an issue. It is invisioned that very few users would want to scrape an entire class an also scrape only specific families from that class simultanious._

A `cazy_dictionary.json` has been created and packaged within the `cazy_webscraper`. This allows users to use a variety of synonoms for the CAZy classes, for example both "GH" and "Glycoside-Hydrolases" are accepted as synonoms for "Glycoside Hydrolases (GHs)". This dictionary is packaged within the `scraper/file_io` directory. If you having issues with the scraper retrieving the list of CAZy classes that are written under 'classes' in the configuration file please check the dictionary first to see the full list of accepted synonoms. If you are comfortable modifying json files then feel free to add your own synonoms to the dictionary.

Under the each of the specific class names listed in the configuration file list the names of specific **families** to be scraped from that class. You do not have to list the class of the families to be scraped under 'classes' as well, this is handled by the web scraper.

Write the true name of the family not only it's number, for example **GH1** is excepted by **1** is not. Each family must be listed on a separate line and the name surrounded by double or single quotation marks. For example:

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
