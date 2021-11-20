=========================
``cazy_webscraper`` cache
=========================

To facilitate reproducibility of scrapping CAZy, ``cazy_webscraper`` logs each scrape within the compiled 
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
* ``<db_name>_<date-time>_connection_failures.log`` which contains :

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cache files when retrieving data from UniProt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When retreiving data from UniProt several cache files are produced, and the number depends on the size of the 
dataset retrieval:



^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cache files when retrieving data from GenBank
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When retrieving protein sequences from GenBank two cache file are produced:

* ``genbank_sequence_retrieved.txt`` contains the GenBank protein accessions for which a protein sequence **was** successfully retrieved
* ``genbank_no_sequence.txt`` contains the GenBank protein accessions for which a protein sequence was not retrieved, as well as the possible reason

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cache files when retrieving data from PDB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




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
