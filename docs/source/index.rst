.. cazy_webscraper documentation master file, created by
   sphinx-quickstart on Fri Nov 20 15:33:10 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root ``toctree`` directive.

Welcome to cazy_webscraper's documentation!
===========================================

| For all the latest updates, and development progress make sure to check the `GitHub repository <https://github.com/HobnobMancer/cazy_webscraper>`_

.. image:: https://img.shields.io/badge/Version-v1.0.2-yellowgreen
   :target: https://github.com/HobnobMancer/cazy_webscraper
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4300858.svg
   :target: https://doi.org/10.5281/zenodo.4300858
.. image:: https://img.shields.io/badge/Licence-MIT-brightgreen
   :target: https://img.shields.io/badge/Licence-MIT-brightgreen
.. image:: https://circleci.com/gh/HobnobMancer/cazy_webscraper.svg?style=shield
   :target: https://circleci.com/gh/HobnobMancer/cazy_webscraper
.. image:: https://codecov.io/gh/HobnobMancer/cazy_webscraper/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/HobnobMancer/cazy_webscraper
.. image:: https://img.shields.io/badge/Python-v3.8.---orange
   :target: https://www.python.org/about/
.. image:: https://img.shields.io/badge/Bioinformatics-EASTBio-ff69b4
   :target: http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019
  
The ``cazy_webscraper`` is a Python3 package for the automated retrieval of protein data from the 
`CAZy <http://wwww.cazy.org/>`_ database. This program is free to use under the MIT license, 
preferably with proper recognition is given.

The ``cazy_webscraper`` retrieves protein data from CAZy, writing out the data to a dataframe in 
the same manner as the data is presented in the CAZy website. If enables, the webscraper will also 
retrieve the protein sequence of the scraped CAZymes from GenBank, writing out the sequences in 
FASTA format. Additionally, if enabled the webscraper will retrieve all protein structures for each 
scraped CAZyme from the Research Collaboratory for Structural Bioinformatics (RCSB) Protein Data 
Bank (`PDB <https://www.rcsb.org/>_).

The ``cazy_webscraper`` can be configured to scrape the entire database, selected Classes and/or 
selected CAZy families. Additionally, the retrieved data can be separated/split by CAZy class, CAZy 
family or not at all, thus gathering all data into a single dataframe.

For details and updates on development, please checkout the GitHub repository.


.. toctree::
   :maxdepth: 4
   
   configuration_scraper
   genbank
   pdb
   license


Requirements
------------------

POISx or Mac OS, or a Linux emulator
Python version 3.8+
Internet access
The python libraries listed within ``requirements.txt``, packaged within the ``cazy_webscraper``


Installation
-------------


At the present moment ``cazy_webscraper`` only supports installation from source. In time we aim 
to include installation via bioconda.

[1] First clone the GitHub repository. This can be done at the command-line with the command:

.. code-block:: bash

   git clone https://github.com/HobnobMancer/cazy_webscraper  

[2] Use the Python package manager ``pip`` to install ``cazy_webscraper``.

.. code-block:: bash

   pip3 install -e <path_to_dir_containing_the_setup.py_file>

Do **not** forget the `-e-` flag when invoking pip. The `-e-` flag installs the package as an 
'executable'. If you do not use the `-e-` flag you will run into *constant* issues when trying to 
run ``cazy_webscraper``.

For those new to using the command-line and installing packages using `pip` here are the lines of 
code you can copy and invoke one at a time into your command-line (please note this is for `bash` 
command-line).

.. code-block:: bash

   git clone https://github.com/HobnobMancer/cazy_webscraper
   cd cazy_webscraper  # this line changes the directory
   pip3 install -e .  # the dot represents look for the setup.py file in the current directory



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


**Dataframes:**

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

**GenBank synonyms:**
Often multiple GenBank accession numbers are listed for a given CAZyme within CAZy. However, only the 
first listed accession number is hyperlinked to the GenBank database. Examination of the other listed 
synonyms (referred to as **genbank synonyms** in the webscraper) shows that these GenBank synonyms are 
the result of submission of identical protein sequences, splice site and protein isoforms. It has been 
interpreted that it is the first GenBank accession that is listed and hyperlinked to GenBank is the accession 
number of protein sequence which was used by CAZy to catalogue the CAZyme and the GenBank synonyms were 
identified and listed by having extremely high sequence identity to the catalogued CAZyme.

Therefore, the webscraper writes only the first GenBank accession listed for each CAZyme in the resulting dataframe. 
The remaining GenBank synonyms are written out to a JSON file, keyed by the first GenBank accession given for each CAZyme, and valued 
by a list of GenBank synonyms. If no GenBank synonyms are retrieved for a CAZyme then the CAZyme's GenBank accession is **not** 
written out to the JSON GenBank synonyms file.

The GenBank synonyms file is written out to the same directory as specificed for the dataframes. 
Additionally, the data is split as is specified for the dataframes.


**Protein sequences:**

If enabled, the protein sequence of the scraped CAZymes are retrieved from GenBank are retrieved in 
the FASTA format and can be written to STDOUT to facilitate piping to a subsequent program or 
written out to disk, within a specified directory.

**Protein structures:**

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

If there are additional features you wish to be added, having consistent problemms with the 
scraper, or want to contribute please raise and issue at the GitHub repository.
