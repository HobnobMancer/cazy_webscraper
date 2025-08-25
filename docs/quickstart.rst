==============================
Quickstart
==============================

For all subcommands, when no options are specified, data is retrieved for all CAZymes in the local CAZyme database, and all available data is retrieved from external databases.

-----------------------
Downloading all of CAZy
-----------------------

To download the entire CAZy dataset, and save the data set to the current working directory with the file name 
``cazy_webscraper_<date>_<time>.db``, use the following command structure: 

.. code-block:: bash

   cazy_webscraper download <user_email> [options]

-------------------------------
Retrieve NCBI GenBank sequences
-------------------------------

To supplement an existing local CAZyme database with protein sequences from NCBI GenBank, use the following command structure:

.. code-block:: bash

   cazy_webscraper get_genbank_seqs <path to local cazyme db> <user_email> [options]


----------------------------
Retrieve NCBI taxonomic data 
----------------------------

To update an existing local CAZyme database with new taxonomic classifications from NCBI, use the following command structure:

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs <path to local cazyme db> <user_email> [options]

--------------------------------
Retrieve UniProt data
--------------------------------

To update an existing local CAZyme database with new UniProt data, use the following command structure:

.. code-block:: bash

   cazy_webscraper get_uniprot_data <path to local cazyme db> [options]

--------------------------------
Retrieve PDB structure files
--------------------------------

To download protein structure files from RCSB PDB for proteins in an existing local CAZyme database, use the following command structure:

.. code-block:: bash

   cazy_webscraper get_pdb_structures <path to local cazyme db> [options]

--------------------------------
Retrieve GTDB taxonomic data
--------------------------------

To update an existing local CAZyme database with new GTDB taxonomic classifications, use the following command structure:

.. code-block:: bash

   cazy_webscraper get_gtdb_taxs <path to local cazyme db> [options]

-----------------------------------------
Extract data from a local CAZyme database
-----------------------------------------

To query a local CAZyme database and output results in a variety of formats, use the following command structure:

.. code-block:: bash

   cazy_webscraper extract_data <path to local cazyme db> [options]
