===========================
Caching and Using the Cache
===========================

To facilitate the reproducibility of scrapping CAZy, ``cazy_webscraper`` logs each scrape within the compiled 
CAZyme datbase and writes out cache files.

This pages walks through the data logged in the CAZyme database and the different cache files written by 
``cazy_webscraper``, as well as how to define a different cache directory, instead of the default.

--------------
The Logs table
--------------

The database built by ``cazy_webscraper`` contains a table called 'Logs'. This table logs every 
scrape of CAZy, UniProt and GenBank which added data to the database.

The table contains the following columns and data:

* **log_id:** Autoincrement ID number
* **date:** Date scrape was initated (in ISO format)
* **time:** Time scrape was initated (in ISO format)
* **database:** Name of the external database from which data was retrieved (i.e. 'CAZy', 'UniProt' or 'GenBank')
* **retrieved_annotations:** List of annotation types retrieved (e.g. 'EC number, PDB accession, Sequence')
* **classes:** CAZy classes for which data was retrieved
* **families:** CAZy families for which data was retrievedclasses
* **kingdoms:** taxonomy Kingdoms filteres applied
* **genera_filter:** taxonomy Kingdoms filteres applied
* **species_filter:** taxonomy Kingdoms filteres applied
* **strains_filter:** taxonomy Kingdoms filteres applied
* **ec_filter:** EC fliters applied to retrieve data for CAZymes annotated with the specified EC numbers (only applies to retrieval of data from UniProt, GenBank and PDB)
* **cmd_line:** Reproduction of the command line arguments passed to ``cazy_webscraper``.

The 'Logs' table allows any one who uses the database to see how the dataset was compiled.

-----------
Cache files
-----------

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cache files when retrieving data from CAZy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``cazy_webscraper`` writes out 2 cache files. These are:

* ``cazy_db_<date-time>.zip`` which is the txt file downloaded from CAZy containing all CAZy Data
* ``<db_name>_<date-time>_connection_failures.log`` which contains a list of proteins for which data was unsuccessfully parsed, and CAZy families for which the family member population could not be retrieved (if the 'validation' option is used).

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cache files when retrieving data from UniProt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The dataframes retrieved with each query to UniProt are cached.

In addition two JSON files are created:
* ``uniprot_accessions_YYYY-MM-DD_HH-MM-SS.json``
* ``uniprot_data_YYYY-MM-DD_HH-MM-SS.json``

**UniProt accessions:**

The first file (``uniprot_accessions``) contains the UniProt accessions/IDs for each GenBank accession retrieved 
from the local CAZyme database that matches the provided criteria. These UniProt IDs are used to query 
UniProt and retrieve protein data. UniProt cannot be batch queried by GenBank accessions to retrieve protein 
data using ``bioservices``.

The JSON file is keyed by the UniProt accession and is valued by a Python dictionary like structure, 
containing the GenBank accession the corresponding ID of its record in the local CAZyme database. For example: 

.. code-block:: python
    {"A0A1S6JHP8": {"gbk_acc": "AQS71285.1", "db_id": 1225219}

This file can be used to skip the retrieval of UniProt IDs from UniProt (which is the first step performed by ``cw_get_uniprot_data``). To 
do this use the ``--skip_uniprot_accessions`` flag followed by the path to the corresponding ``uniprot_accessions_YYYY-MM-DD_HH-MM-SS.json`` file.

**UniProt data:**

The second json file (``uniprot_data``) contains all data retrieved from UniProt for all proteins in the local 
CAZyme database that match the specified criteria. The data retrieved from UniProt was parsed into a Python dictionary 
which is then dumped into the JSON file.

This file is used for mannually checking the parsing method employed by ``cw_get_uniprot_data`` is working, as well as skipping the 
retrieval of the same dataset from UniProt (for example, if you wanted to recreate a specific CAZyme proteome dataset).

To use the data cached in the ``uniprot_data`` file, using the ``--use_uniprot_cache`` flag, followed by the 
path pointing to the corresponding file. Using this flag, skips the retrieval of protein data from UniProt, and only adds 
data from the cache file into the local CAZyme database.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cache files when retrieving data from GenBank
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When retrieving protein sequences from GenBank three cache file are produced:

* ``genbank_sequence_retrieved_YYYY-MM-DD_HH-MM-SS.txt`` contains the GenBank protein accessions for which a protein sequence **was** successfully retrieved
* ``genbank_no_sequence_YYYY-MM-DD_HH-MM-SS.txt`` contains the GenBank protein accessions for which a protein sequence was not retrieved, as well as the possible reason
* ``gennbank_seqs_YYYY-MM-DD_HH-MM-SS.json``, which contains the GenBank accessions of proteins matching the provided criteria and the retrieved protein sequences.

.. NOTE::
    Future features for ``cazy_webscraper`` will include the ability to use this cached protein sequences, to (i) skip the 
    retrieval of data from GenBank and (ii) facilitate the reproduction of local CAZyme databases.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cache files when retrieving data from PDB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One cache file is created when using ``cazy_webscraper`` to retrieval protein structure files from PDB: ``pdb_retrieval_YYYY-MM-DD_HH-MM-SS.txt``, which 
lists the PDB accessions of all files that were successfully downloaded from PDB using ``cazy_webscraper`` and ``BioPython``.

---------------
Cache directory
---------------

By default ``cazy_webscraper`` creates the cache directory in the same directory that the datbase is created, and 
with the name ``.cazy_webscraper``.

To use a different cache directory instead add the ``--cache_dir`` flag, followed by the path to the cache directory.

.. NOTE::
    The cache directory does not need to already exist, ``cazy_webscraper`` will build the cache directory 
    and all it's parent directories.

If the target cache directory already exists, by default ``cazy_webscraper`` will delete the content already 
present in the already existing cache directory. To not delete the exsiting directory content add the 
``nodelete_cache`` flag.
