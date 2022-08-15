.. cazy_webscraper documentation master file, created by
   sphinx-quickstart on Fri Nov 20 15:33:10 2020.

===========================================
Welcome to cazy_webscraper's documentation!
===========================================

.. image:: cazy_web_logo.svg
   :scale: 50 %
   :alt: cazy_webscraper logo, host organisations and funding
   :align: center

| For latest updates and development progress, please check the `GitHub repository <https://github.com/HobnobMancer/cazy_webscraper>`_

-----------------
Build Information
-----------------

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


--------
``PyPI``
--------

.. image:: https://img.shields.io/pypi/v/cazy_webscraper.svg?style=flat-square
    :target: https://pypi.python.org/pypi/cazy_webscraper
.. image:: https://img.shields.io/pypi/dm/cazy_webscraper?label=Pypi%20downloads
   :target: https://pypi.org/project/cazy-webscraper/

------------
``bioconda``
------------

.. image:: https://anaconda.org/bioconda/cazy_webscraper/badges/installer/conda.svg?style=flat-square
     :target: https://conda.anaconda.org/bioconda
.. image:: https://anaconda.org/bioconda/cazy_webscraper/badges/version.svg?style=flat-square
    :target: https://anaconda.org/bioconda/cazy_webscraper
.. image:: https://anaconda.org/bioconda/cazy_webscraper/badges/latest_release_date.svg?style=flat-square
     :target: https://anaconda.org/bioconda/cazy_webscraper
.. image:: https://img.shields.io/conda/dn/bioconda/cazy_webscraper?label=Bioconda%20downloads
   :target: https://bioconda.github.io/user/install.html

-------------------
``cazy_webscraper``
-------------------

``cazy_webscraper`` is a Python3 package for the automated retrieval of Carbohydrate-Active enZyme (CAZyme) data from the `CAZy <http://wwww.cazy.org/>`_ database. This program is free to use under the MIT license, and we kindly request that, if you use this program or Python package, you cite it as indicated below.

   Hobbs, Emma E. M.; Pritchard, Leighton; Chapman, Sean; Gloster, Tracey M. (2021): cazy_webscraper Microbiology Society Annual Conference 2021 poster. figshare. Poster. https://doi.org/10.6084/m9.figshare.14370860.v7 

``cazy_webscraper`` retrieves data from CAZy, writing it to a local SQLite3 file (typically taking 10-15 minutes to scrape the entirety of CAZy). 

**Additionally, ``cazy_webscraper`` can:**

* Retrieve the protein data from `UniProt <https://www.uniprot.org/>`_ for CAZymes in the local database. This data includes:

   * UniProt accession
   * Protien name
   * Protein amino acid sequence
   * EC numbers
   * PDB accessions

* Retrieve protein sequences from NCBI GenBank for CAZymes in the local database.
* Write out protein sequences retrieved from UniProt and NCBI in FASTA format, and build a local BLAST database.
* Retrieve protein structures from the Research Collaboratory for Structural Bioinformatics (RCSB) Protein Data Bank, `PDB <https://www.rcsb.org/>`_, for CAZymes in the local database.
* Be configured to scrape the entire CAZy database, or recover only CAZymes filtered by user-supplied criteria, such as CAZy classes, CAZy (sub)family, or taxonomy. 
* Retrieve the latest taxonomic classifications (including the complete lineage) from the NCBI Taxonomy database

----------
Quickstart
----------

We have produced a "Getting Started With ``cazy_webscraper``" `poster <https://hobnobmancer.github.io/cazy_webscraper/getting_started_poster.pdf>`_.

To download the entire CAZy dataset, and save the data set to the current working directory with the file name 
``cazy_webscraper_<date>_<time>.db``, use the following command structure:  

.. code-block:: bash

   cazy_webscraper <user_email>

.. NOTE::
   The user email address is a requirement of NCBI. NCBI is queried to identify the currect source organism 
   for a given protein, when multiple source organisms are retrieved from CAZy for a single protein. 
   For more information please see the `NCBI Entrez <https://www.ncbi.nlm.nih.gov/books/NBK25497/>`_ documentation.

---------------
Command summary
---------------

Below are the list of commands (excluding required and optional arguments) included in ``cazy_webscraper``.

**CAZy**

To retrieve data from CAZy and compile and SQLite database using ``cazy_webscraper`` command.

**UniProt**

To retrieve protein data from UniProt, use the ``cw_get_uniprot_data`` command.

The following data can be retrieved:
- UniProt accession
- Protein name
- EC numbers
- PDB accession
- Protein sequences

**GenBank**

- To retrieve protein sequences from GenBank use the ``cw_get_genbank_seqs`` command.
- To retrieve the latest taxonomic classifications from NCBI Taxonomy using the ``cw_get_ncbi_taxs`` command.

**Extract sequences**

To extract GenBank and/or UniProt protein sequences from a local CAZyme database, use the ``cw_extract_db_seqs`` command.

**PDB**

To protein structure files from PDB use the ``cw_get_pdb_structures`` command.

**NCBI taxonomies**

Retrieve the latest taxonomic classifications (including the complete lineage from kingdom to strain) using the ``cw_get_ncbi_taxs`` command.

**GTDB taxonomies**

Retrieve the latest taxonomic classifications (incluidng the complete lineage from kingdom to strain) from the GTDB database using the ``cw_get_gtdb_taxs`` command.

**Interrogate the database**

To interrogate the database, use the ``cw_query_database`` command.

-------------
Best practice
-------------

When performing a series of many automated, repeated calls to a server it is polite to do this when internet traffic is lowest *at the server*. This is typically at the weekend and overnight.

When using ``cazy_webscraper`` to retrieve data from UniProt, NCBI or PDB, the webscraper can appear 
to run slowly but this may be due to bandwidth at the database server, or server speed. 
``cazy_webscraper`` provides a progress bar to reassure the user that the webscraper is working.

.. WARNING::
   Please **do not** perform a retrieval of UniProt, NCBI and/or PDB data for the entire CAZy dataset, unless 
   absolutely unavoidable. Retrieving the data from any of these exteranl databases for the entire CAZy 
   dataset will take several hours and may unintentionally deny the service to others.

-------------
Documentation
-------------

For details and updates on development, please consult the `GitHub repository <https://github.com/HobnobMancer/cazy_webscraper>`_.

.. toctree::
   :maxdepth: 2
   
   installation
   quickstart
   usage
   tutorial
   database
   uniprot
   uniprottutorial
   genbank
   genbanktutorial
   sequence
   sequencetutorial
   pdb
   pdbtutorial
   ncbitax
   ncbitaxtutorial
   genomes
   genomestutorial
   gtdb
   gtdbtutorial
   api
   apitutorial
   cache
   integrate
   contributing
   license


--------------------------
Citing ``cazy_webscraper``
--------------------------
 
If you use ``cazy_webscraper`` in your work *please* do cite our work (including the provided DOI), as well as the specific version of the tool you use. This is not only helpful for us as developers to get our work out into the world, but it is also **essential for the reproducibility and integrity of scientific research**. Citation:
   
   Hobbs, Emma E. M.; Pritchard, Leighton; Chapman, Sean; Gloster, Tracey M. (2021): cazy_webscraper Microbiology Society Annual Conference 2021 poster. FigShare. Poster. https://doi.org/10.6084/m9.figshare.14370860.v7 


----------------------
Development and issues
----------------------

If there are additional features you wish to be added, you have problems with the scraper, or would like to contribute please raise an issue at the `GitHub repository <https://github.com/HobnobMancer/cazy_webscraper>`_.

* Issues Page: `https://github.com/HobnobMancer/cazy_webscraper/issues <https://github.com/HobnobMancer/cazy_webscraper/issues>`_
