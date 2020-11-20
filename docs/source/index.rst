.. cazy_webscraper documentation master file, created by
   sphinx-quickstart on Fri Nov 20 15:33:10 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root ``toctree`` directive.

Welcome to cazy_webscraper's documentation!
===========================================

| Version v0.1.2 2020/11/20
| DOI: ---
| `GitHub repository <https://github.com/HobnobMancer/cazy_webscraper>`_

The ``cazy_webscraper`` is a Python3 package for the automated retrieval of protein data from the [CAZy](http://wwww.cazy.org/) database. This program is free to use under the MIT license when proper recognition is given.

The CAZy database is one of the largest databases for cataloging Carbohydrate Active enZymes, and is one of the most frequently visited protein data bases for researchers in this field. However, CAZy provides not method of automated retrieval of data, especially for large dataset requests. The ``cazy_webscraper`` provides a way for the automated retrieval of large datasets from the CAZy database. The retrieved protein data is written out into dataframes, and protein sequences written out to FASTA files.

The ``cazy_webscraper`` is configurable to scrape specific CAZy classes and/or families, and determine how the data split when written out. The data can be split:

- Per family: a single dataframe is created per scraped CAZy family
- Per class: a singel datafarme is created per scraped CAZy family
- Not split: a single dataframe containing all data scrapend from CAZy is created


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

The easiest method for install the ``cazy_webscraper`` is to use pips.  
``pip3 install -e <path_to_dir_containing_the_setup.py_file>``  

Then install remaining requirements:  
``pip3 install -r <path_to_requirements.txt_file>``  

In both commands to **do not** forget the additional pips options (``-e`` and ``-r``)!


Best practice
--------------

When performing a series of many, automated, repeated calls to a server, such as is performed by the ``cazy_webscraper``, it is best practise to do this during the period of the day when traffic is lowest. This typically includes over the weekend and overnight.

The webscraper can appear to run slowly but this is due to the limited access speed of the CAZy server. When a CAZy family is being parsed by the scraper, and protein records are being retrieved for the CAZy family a progress bar is produced in the terminal to provide an indicator the webscraper is working. However, expect an entire scrape of the CAZy database to take severak hours.


Output
-------

At the moment the webscraper retrieves the protein data as presented in the CAZy database, in the table formate as viewed in a webrowser.
This data is then written out a dataframe with the same headings as present in CAZy, with the exception of the additional column 'CAZy family' which lists the proteins CAZy family or subfamily as appropriate. This is in case multiple families are scraped and the proteins are stored in a single dataframe together.


Configuration
--------------

The scraping of CAZy is entirely configurable. An entire class and/or specific families can be scraped, and 
the ability to scrape subfamilies can be enabled or disabled.

For more detailed see the configuration section of the documetation.


Development and issues
-----------------------

If there are additional features you wish to be added, or having consistent problemms with the scraper, please raise and issue at the GitHub repository.
