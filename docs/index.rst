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

.. image:: https://img.shields.io/badge/Version-v3.0.0-yellowgreen
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

``cazy_webscraper`` retrieves data from CAZy, writing it to a local SQLite3 file, typically taking 5-15 minutes to scrape the entirety of CAZy. 

**Additionally, ``cazy_webscraper`` can also supplement the local CAZyme database with data from the following public databases:**


* `NCBI <https://www.ncbi.nlm.nih.gov/>`_:
   * Protein sequence
   * Taxonomic classification, including the NCBI Taxonomy ID
   * Genomic assembly accession
* `UniProt <https://www.uniprot.org/>`_:
   * UniProt accession
   * Protein name
   * Protein amino acid sequence
   * EC numbers
   * PDB accessions
* `RCSB PDB <https://www.rcsb.org/>`_:
   * Protein structure files in variety of file formats
* `GTDB <https://gtdb.ecogenomic.org/>`_:
   * Taxonomic classification, including the GTDB Taxonomy ID

----------
Quickstart
----------

To download the entire CAZy dataset, and save the data set to the current working directory with the file name 
``cazy_webscraper_<date>_<time>.db``, use the following command structure:  

.. code-block:: bash

   cazy_webscraper download_cazy <user_email>

.. NOTE::
   The user email address is a requirement of NCBI. NCBI is queried to identify the currect source organism 
   for a given protein, when multiple source organisms are retrieved from CAZy for a single protein. 
   For more information please see the `NCBI Entrez <https://www.ncbi.nlm.nih.gov/books/NBK25497/>`_ documentation.

---------------
Command summary
---------------

Below are the subcommands (excluding required and optional arguments) included in ``cazy_webscraper``.

.. list-table::
   :header-rows: 1

   * - Command
     - Description
   * - ``download_cazy``
     - Download CAZy data and create a local CAZyme database.
   * - ``get_ncbi_seqs``
     - Supplement an existing local CAZyme database with protein sequences from NCBI.
   * - ``get_ncbi_taxs``
     - Update an existing local CAZyme database with new taxonomic classifications from NCBI.
   * - ``get_ncbi_genomes``
     - Supplement an existing local CAZyme database with genomic assembly data from NCBI.
   * - ``get_uniprot_data``
     - Supplement an existing local CAZyme database with data from UniProt.
   * - ``get_pdb_structures``
     - Download protein structure files from RCSB PDB, and add the experimental method and
       resolution of each structure to an existing local CAZyme database.
   * - ``get_pfams``
     - Supplement an existing local CAZyme database with Pfam domain annotations from InterPro.
   * - ``get_gtdb_taxs``
     - Add GTDB taxonomic classifications for the genomes in an existing local CAZyme database.
   * - ``extract_data``
     - Write data held in an existing local CAZyme database out to CSV, TSV, JSON, JSON Lines,
       FASTA or a BLAST database.

``download_cazy`` is the one subcommand that *builds* a local CAZyme database; every ``get_*``
subcommand adds data to a database that already exists.

.. note::
   ``--version``, ``--citation``, ``-l``/``--log``, ``--sql_echo`` and ``-v``/``--verbose`` belong
   to ``cazy_webscraper`` itself and must be given **before** the subcommand, for example
   ``cazy_webscraper -v get_ncbi_seqs ...``.

.. warning::
   The version 2 ``cw_query_database`` and ``cw_extract_db_seqs`` commands are now the single
   ``extract_data`` subcommand; see :doc:`migration` for the mapping. ``cw_get_db_schema`` has
   **not yet been migrated**.

--------------------------
Upgrading from version 2
--------------------------

Version 3 replaces the ten separate ``cw_*`` commands of version 2 with subcommands of a single
``cazy_webscraper`` command, and changes the schema of the local CAZyme database. A CAZyme database
built with version 2 **cannot be used with version 3** and must be rebuilt.

See :doc:`migration` for the full command mapping, the renamed arguments, and a migration
checklist.

-------------
Documentation
-------------

For details and updates on development, please consult the `GitHub repository <https://github.com/HobnobMancer/cazy_webscraper>`_.

.. toctree::
   :maxdepth: 2
   
   installation
   quickstart
   migration
   01scrapeCazy
   02genbankSeqs
   03ncbiTaxs
   uniprot
   sequence
   sequencetutorial
   pdb
   pdbtutorial
   pfam
   genomes
   genomestutorial
   gtdbtax
   gtdbtaxtutorial
   api
   apitutorial
   database
   cache
   integrate
   contributing
   citation
   license
   citation
   license


--------------------------
Citing ``cazy_webscraper``
--------------------------
 
If you use ``cazy_webscraper`` in your work *please* do cite our work (including the provided DOI), as well as the specific version of the tool you use. This is not only helpful for us as developers to get our work out into the world, but it is also **essential for the reproducibility and integrity of scientific research**.

**Citation:**

   Hobbs, E. E. M., Gloster, T. M., and Pritchard, L. (2023) 'cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets', Microbial Genomics, 9(8), https://doi.org/10.1099/mgen.0.001086

This paper includes a full description of the operation and examples of use.

``cazy_webscraper`` depends on a number of tools. To recognise the contributions that the 
authors and developers have made, please also cite the following:

When making an SQLite database:
   Hipp, R. D. (2020) SQLite, available: https://www.sqlite.org/index.html.

Retrieving taxonomic, genomic or sequence data from NCBI:
   Cock, P.J.A., Antao, T., Chang, J.T., Chapman, B.A., Cox, C.J., Dalke, A., et al (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics, 25(11), 1422-1423.
   
   Wheeler,D.L., Benson,D.A., Bryant,S., Canese,K., Church,D.M., Edgar,R., Federhen,S., Helmberg,W., Kenton,D., Khovayko,O. et al (2005) Database resources of the National Centre for Biotechnology Information: Update, Nucleic Acid Research, 33, D39-D45

Retrieving data from UniProt:
   Cokelaer, T., Pultz, D., Harder, L. M., Serra-Musach, J., Saez-Rodriguez, J. (2013) BioServices: a common Python package to access biological Web Services programmatically, Bioinformatics, 19(24), 3241-3242.

Downloading protein structure files from RSCB PDB:
   Berman, H.M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T.N., Weissig, H., et al (2022) The Protein Data Bank, Nucleic Acids Research, 28(1), 235-242.
   
   Hamelryck, T., Manderick, B. (2003), PDB parser and structure class implemented in Python. Bioinformatics, 19 (17), 2308–2310

Retrieving and using taxonomic data from GTDB:
   Parks, D.H., Chuvochina, M., Rinke, C., Mussig, A.J., Chaumeil, P., Hugenholtz, P. (2022) GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy, Nucleic Acids Research, 50(D1), D785-D794.

----------------------
Development and issues
----------------------

If there are additional features you wish to be added, you have problems with the scraper, or would like to contribute please raise an issue at the `GitHub repository <https://github.com/HobnobMancer/cazy_webscraper>`_.

* Issues Page: `https://github.com/HobnobMancer/cazy_webscraper/issues <https://github.com/HobnobMancer/cazy_webscraper/issues>`_
