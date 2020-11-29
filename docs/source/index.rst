.. cazy_webscraper documentation master file, created by
   sphinx-quickstart on Fri Nov 20 15:33:10 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root ``toctree`` directive.

Welcome to cazy_webscraper's documentation!
===========================================

| Version v1.0.0 2020/11/30
| DOI: ---
| `GitHub repository <https://github.com/HobnobMancer/cazy_webscraper>`_

The ``cazy_webscraper`` is a Python3 package for the automated retrieval of protein data from the 
[CAZy](http://wwww.cazy.org/) database. This program is free to use under the MIT license, 
preferably with proper recognition is given.

The ``cazy_webscraper`` retrieves protein data from CAZy, writing out the data to a dataframe in 
the same manner as the data is presented in the CAZy website. If enables, the webscraper will also 
retrieve the protein sequence of the scraped CAZymes from GenBank, writing out the sequences in 
FASTA format. Additionally, if enabled the webscraper will retrieve all protein structures for each 
scraped CAZyme from the Research Collaboratory for Structural Bioinformatics (RCSB) Protein Data 
Bank [(PDB)](https://www.rcsb.org/).

The ``cazy_webscraper`` can be configured to scrape the entire database, selected Classes and/or 
selected CAZy families. Additionally, the retrieved data can be separated/split by CAZy class, CAZy 
family or not at all, thus gathering all data into a single dataframe.


.. toctree::
   :maxdepth: 4
   
   configuration_scraper
   genbank
   pdb
   notebooks
   license


Requirements
------------------

POISx or Mac OS, or a Linux emulator
Python version 3.8+
Internet access
The python libraries listed within ``requirements.txt``, packaged within the ``cazy_webscraper``


Installation
-------------


**Quick and Easy:**
The easiest method for installing the ``cazy_webscraper`` is to use Conda, using the following 
command at the command-line in the terminal:  
``conda asdasd cazy_webscraper``
This method installs the full ``cazy_webscraper`` and all dependencies.

If Conda is not installed, please see the Conda website for installation 
[instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).


**Alternative, script access**
Alternatively, for easier access to Python scripts that make up the ``cazy_webscraper`` first clone 
the webscrapers GitHub repository, then install the webscraper using pips:
``git clone https://github.com/HobnobMancer/cazy_webscraper``  
``pip3 install -e <path_to_dir_containing_the_setup.py_file>``  

Then install remaining requirements:  
``pip3 install -r <path_to_requirements.txt_file>``  

In both commands to **do not** forget the additional pips options (``-e`` and ``-r``)!

To write the webscraper repository to a specific directory use the following command:    
``git clone https://github.com/HobnobMancer/cazy_webscraper> <path_to_dir>``  
The directory to which the path points will form the root of the local copy of the repository.  


Best practice
--------------

When performing a series of many, automated, repeated calls to a server it is best practise to do 
this during the period of the day when traffic is lowest. This typically includes over the weekend 
and overnight. Therefore, when scraping entire CAZy database, entire class(es), and/or multiple 
CAZy families it is advised perform these scrapes over night and/or over the weekend.

The webscraper can appear to run slowly but this is due to the limited access speed of the CAZy 
server. When a CAZy family is being parsed by the scraper, and protein records are being retrieved 
for the CAZy family a progress bar is produced in the terminal to provide an indicator that 
webscraper is working. **However, expect an entire scrape of the CAZy database to take several hours.**


Output
-------

The basic function of the ``cazy_webscraper`` is to retrieve the protein data stored within and 
presented in the CAZy database. The data is written out to a dataframe with the same headings as 
found in CAZy, to reflect the way CAZy presented data in its webpages. The resulting dataframe 
includes the additional column "CAZy family", which includes the CAZy family/subfamily under which 
the respective CAZy is catalogued. Therefore, the resulting dataframe of the webscraper contains 
the following columns:

* Protein_name
* CAZy_family
* EC#
* Source_organism
* GenBank
* UniProt
* PDB/3D

The scraped CAZymes can be written to a single dataframe, or separated out into different 
dataframes depending on how the data has been set to be split:

* Per family: a single dataframe is created per scraped CAZy family
* Per class: a singel datafarme is created per scraped CAZy family
* Not split: a single dataframe containing all data scrapend from CAZy is created

The dataframes can be written out to a specified directory or written to STDOUT to facilitate 
piping to a subsequent program. If the dataframes are written to the disk they are saved as 
**.csv** files.

If enabled, the protein sequence of the scraped CAZymes are retrieved from GenBank are retrieved in 
the FASTA format and can be written to STDOUT to facilitate piping to a subsequent program or 
written out to disk, within a specified directory.

If enabled, the protein structures will be written out to the disk, to a specified directory. The 
protein structures cannot be written to STDOUT due to using the ``BioPython`` module ``PDB``, 
which currently does not facilitate writing out the protein structures to STDOUT. The format of the 
the structure file is specified at the command line.


Quick Start
--------------

To invoke the webscraper with its default functionality simply call the webscraper at the command 
line:  
``cazy_webscraper``  

The default behaviour of the scraper is:

* Scrape the entire CAZy databases
* Not to split/separate the data, producing a single dataframe
* Write the resulting dataframe to standard out (STDOUT)
* Not to retrieve subfamilies (members of subfamilies will be retrieved but their parent family be listed)
* Not to retrieve FASTA files from GenBank
* Not to retrieve protein sequences from PDB


Configuration
--------------

The scraping of CAZy is entirely configurable to suit your purpose. Entire classes and/or specific 
families can be scraped, and the ability to scrape subfamilies can be enabled or disabled.

For more detailed see the configuration section of the documetation.


Development and issues
-----------------------

If there are additional features you wish to be added, or having consistent problemms with the 
scraper, please raise and issue at the GitHub repository.
