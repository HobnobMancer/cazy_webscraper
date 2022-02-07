====================================
Interrogating the data using the API
====================================

The data stored in the local CAZyme database can be interrogated using SQL. ``cazy_webscraper`` also 
includes an API, which can be used to interrogate the data in the local CAZyme database and write out the retrieved data 
in ``JSON`` and/or ``CSV`` format.

By default ``cazy_webscraper`` only includes the GenBank accessions of proteins matching the provided 
criteria, but the inclusion of additional data (such as protein squences, UniProt accessions, EC numbers, etc) 
is fully customisable.

-----------
Quick Start
-----------

To retrieve the GenBank accession of all CAZymes stored in the local CAZyme database, use the following command 
structure:

.. code-block:: bash

   cw_query_database <path to local CAZyme db> <desired file formats>

.. NOTE::
   The ``cw`` prefix on command is an abbreviation of ``cazy_webscraper``.
   
---------------------
Accepted file formats
---------------------

``cw_query_database`` can write the output to a csv or json file.

These are provided as the second arguments to ``cw_query_database``. To write out both a csv and json 
file use both ``csv`` and ``json`` after the path to the local CAZyme database, separted with a single space.

.. code-block:: bash

    cw_query_database <path to local CAZyme db> csv json

.. NOTE::
    The order ``csv`` and ``json`` are written does not matter.

.. WARNING::
    Both ``csv`` and ``json`` are case sensitive.

--------------------
Command line options
--------------------

``database`` - **REQUIRED** Path to a local CAZyme database to add UniProt data to.

``file_types`` - **REQUIRED** List of file formats to export the data in. Currently supported: ``csv`` and ``json``.

``--cache_dir`` - Path to cache dir to be used instead of default cache dir path.

``--cazy_synonyms`` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

``--config``, ``-c`` - Path to a configuration YAML file. Default: None.

``--classes`` - list of classes to retrieve data from.

``--ec_filter`` - List of EC numbers to limit the retrieval of protein data with at least one of the given EC numbers *in the local CAZyme database*.

``--families`` - List of CAZy (sub)families to retrieve UniProt protein data for.

``--genbank_accessions`` - Path to text file containing a list of GenBank accessions to retrieve protein data for. A unique accession per line.

``--genera`` - List of genera to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given genera.

``--include`` - List additional data to include in the output. Multiple fields can be named, separating each with a single space (' '). The accepted fields are:
* 'class' - Include the CAZy class annotations
* 'family' - Include the CAZy family annotations
* 'subfamily' - Include the subfamily class annotations
* 'kingdom' - Include the taxonomic kingdom of the source organism
* 'genus' - Include the genus of the source organism
* 'organism' - Include the scientific name of the source organism
* 'uniprot_acc' - Include the UniProt accession
* 'uniprot_name' - Include the protein name retrieved from UniProt
* 'ec' - Include the EC number annotations
* 'pdb' - Include the PDB accessions
* 'genbank_seq' - Include the GenBank protein sequence
* 'uniprot_seq' - Include the Uniprot protein sequence

``--kingdoms`` - List of taxonomy kingdoms to retrieve UniProt data for.

``--log``, ``-l`` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

``--nodelete``, ``-n`` - When called, content in the existing output  will **not** be deleted. Default: False (existing content is deleted).

``--nodelete_cache`` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

``--outdir``, ``-o`` - Define output directory to write out structure files. Default, write structure files to current working directory.

``--overwrite`` - Overwrite existing structure files with the same PDB accession as files being downloaded. Default false, do not overwrite existing files.

``--retries``, ``-r`` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

``--sql_echo`` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

``--species`` - List of species (organsim scientific names) to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given species.

``--strains`` - List of species strains to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given strains.

``--timeout``, ``-t`` - Connection timout limit (seconds). Default: 45.

``--uniprot_accessions`` - Path to text file containing a list of UniProt accessions to retrieve protein data for. A unique accession per line.

``--verbose``, ``-v`` - Enable verbose logging. This does **not** set the SQLite engine ``echo`` parameter to True. Default: False.



-----------
Basic Usage
-----------

The command-line options listed above can be used in combination to customise the retrieval of data for proteins of interest. Some options (e.g. ``--families`` and ``--classes``) define the broad group of proteins for which structure files are retrieved, others (e.g. ``--species``) are used to filter and fine-tune the protein dataset for which structure files are retrieved.

The ``--classes``, ``--families``, ``--kingdoms``, ``--genera``, ``--species``, and ``--strains`` filteres are applied 
in the exactly same for retrieving data from CAZy, UniProt, GenBank, and PDB as for the extraction of data via the API. Examples of using these flags 
can be found in the ``cazy_webscraper`` and ``cw_query_database`` tutorial in this documentation.

.. NOTE::
    To retrieve fata for members of specific CAZy subfamilies, list the subfamilies after the ``--families`` 
    flag.


--------------------------------------------
Retrieving data from a local CAZyme database
--------------------------------------------

The command for using ``cazy_webscraper`` for interrogating the local CAZyme database and extract data is ``cw_query_databsae``.
