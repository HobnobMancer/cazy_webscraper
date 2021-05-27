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
Bank, `PDB <https://www.rcsb.org/>`_.

The ``cazy_webscraper`` can be configured to scrape the entire database, selected Classes and/or 
selected CAZy families. Additionally, the retrieved data can be separated/split by CAZy class, CAZy 
family or not at all, thus gathering all data into a single dataframe.

For details and updates on development, please checkout the GitHub repository.


.. toctree::
   :maxdepth: 4
   
   configuration_scraper
   configuration_tutorials
   genbank
   pdb
   license


Reference
----------------
 
If you use `cazy_webscraper` in your work *please* do cite our work (including the provided DOI), as well as citing the specific version you use. This is not only helpful for us as developers to get out work out into the world, but it is also **essential for the reproducibility and integrity of scientific research**.  

**Citation:** Hobbs, Emma E. M.; Pritchard, Leighton; Chapman, Sean; Gloster, Tracey M. (2021): cazy_webscraper Microbiology Society Annual Conference 2021 poster. figshare. Poster. https://doi.org/10.6084/m9.figshare.14370860.v7 


Requirements
------------------

POISx or Mac OS, or a Linux emulator
Python version 3.8+
Internet access
The python libraries listed within ``requirements.txt``, packaged within the ``cazy_webscraper``


Installation
-------------

**Quick and easy**

You can install ``cazy_webscraper`` via `pip  <https://pypi.org/project/cazy-webscraper/>`_:  

.. code-block:: bash

   pip3 install cazy-webscraper

**Install from source**

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
code you can copy and invoke one at a time into your command-line (please note these comamnds are written in `bash`).

.. code-block:: bash

   git clone https://github.com/HobnobMancer/cazy_webscraper
   cd cazy_webscraper  # this line changes the directory
   pip3 install -e .  # the dot tells the computer to look for the setup.py file in the current directory


Getting started
-------------------

For a quick summary for getting started, checkout the poster: Hobbs, Emma E. M.; Pritchard, Leighton; Gloster, Tracey M.; Chapman, Sean (2021): cazy_webscraper - getting started. figshare. Poster. https://doi.org/10.6084/m9.figshare.14370869.v3 


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
