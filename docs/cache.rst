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



------------------------------------------
Cache files when retrieving data from CAZy
------------------------------------------

``cazy_webscraper`` writes out 2 cache files. These are:

* ``cazy_db_<date-time>.zip`` which is the txt file downloaded from CAZy containing all CAZy Data
* ``<db_name>_<date-time>_connection_failures.log`` which contains a list of proteins for which data was unsuccessfully parsed, and CAZy families for which the family member population could not be retrieved (if the 'validation' option is used).


---------------------------------------------
Cache files when retrieving data from UniProt
---------------------------------------------

``get_uniprot_data`` writes two plain text files into the cache directory:

* ``uniprot_connection_errors_YYYY-MM-DD_HH-MM-SS.txt`` -- requests that failed with a connection
  error
* ``uniprot_failed_ids_YYYY-MM-DD_HH-MM-SS.txt`` -- accessions UniProt returned no usable record
  for

Both are diagnostic: they tell you which proteins were missed so that you can re-run for just
those, by passing the accessions to ``--genbank_accessions``.

.. note::
   Version 2 additionally wrote ``uniprot_accessions_*.json`` and ``uniprot_data_*.json`` cache
   files, which could be replayed with the ``--skip_uniprot_accessions`` and
   ``--use_uniprot_cache`` flags. Version 3 writes neither file, and both flags have been removed.
   See :doc:`migration`.

----------------------------------------------------------
Cache files when retrieving protein sequences from GenBank
----------------------------------------------------------

When retrieving protein sequences from GenBank, all cache files will be contained in the 
``genbank_seq_retrieval.tar`` compressed file. The cache files that will be generated are:

 three cache file are produced:

* ``cazy_webscraper get_ncbi_seqs <cache_time>.fasta`` contains the protein sequences downloaded from GenBank
* ``failed_retrieve_ids`` contains the GenBank accessions containing the IDs of **all** proteins for whom no protein sequence was retrieved
* ``failed_entrez_connection_accessions`` contains the GenBank accessions of proteins whom no protein sequence was retrieved owing to failure to connect to NCBI
* ``invalid_ids`` contains GenBank protein version accessions that were retrieved from CAZy, but are no longer listed in NCBI

.. NOTE::
    ``cazy_webscraper`` includes the option to use cached protein sequences, to (i) skip the 
    retrieval of data from GenBank and (ii) facilitate the reproduction of local CAZyme databases.

-----------------------------------------------------------
Cache files when retrieving taxonomic information from NCBI
-----------------------------------------------------------

When retrieving taxonomic classifications the following cache files are generated:
* ``failed_protein_accs.txt`` - list the version accession of all proteins for which not taxonomy data could be retrieved from NCBI (possibly the protein record has been removed)
* ``lineage_data.json`` - lists the lineage data and the associated protein sequences accessions
* ``ncbi_lineages.json`` - lists the lineage data retrieved from NCBI Taxonomy
* ``protein_ncbi_ids.out`` - the NCBI Protein IDs retrieved from NCBI when querying by protein version accession
* ``tax_ids.out`` - the NCBI Taxonomy IDs retrieved from NCBI when querying by NCBI Protein IDs

When retrieving genomic assembly data a file listing all protein sequence accessions for which no data was retrieved is cached (``failed_protein__accessions.txt``). 
Additionally, the downloaded genomic assemlby feature tables are retained in their compress format and cached.

-----------------------------------------
Cache files when retrieving data from PDB
-----------------------------------------

One cache file is created when using ``cazy_webscraper`` to retrieval protein structure files from PDB: ``pdb_retrieval_YYYY-MM-DD_HH-MM-SS.txt``, which 
lists the PDB accessions of all files that were successfully downloaded from PDB using ``cazy_webscraper`` and ``BioPython``.


---------------------------------------------------------
Cache files when retrieving genomic information from NCBI
---------------------------------------------------------

All genomic assembly feature tables are cached in the cache directory. In addition, two other cache files 
are generated to list proteins and genomes for whom data could not be retrieved from NCBI:

* ``no_assembly_urls.txt`` contains the NCBI Assembly Name of genomic assemlbies for whom a feature table could not be downloaded from the NCBI FTTP server - typically because a feature table is not available
* ``failed_protein_accessions.txt`` contains the NCBI protein version accessions for whom genomic assembly information

-----------------------------------------------------------
Cache files when retrieving taxonomic information from GTDB
-----------------------------------------------------------

The GTDB database dumbs are downloaded and written to the cache directory:

* ``archaea_data-{release_number}.gz``
* ``bacteria_data-{release_number}.gz``

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
