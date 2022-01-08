===================================
Retrieving structure files from PDB
===================================

``cazy_webscraper`` can be used to retrieve protein structure files for PDB accessions in a local CAZyme database from `RSCB PDB database <https://www.rcsb.org/>`_. The downloading of the structure files is handled by the ``BioPython`` module ``Bio.PDB``. 

For specific information of the ``Bio.PDB`` module please see the 
`BioPython documentation <https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ>`_.

.. warning::
        If many PDB structure files are going to retrieved PDB (for example more than 100), it is expected practise to perform the
        opretional outside peak times.

.. note::
    PDB structure files are retrieved for the PDB accessions *in* a local CAZyme database created using ``cazy_webscraper``.

-----------
Quick Start
-----------

To download the protein structure file for all PDB accessions in a local CAZyme database, use the following command structure:

.. code-block:: bash

   cw_get_pdb_structures <path to local CAZyme db> <desired file formats>

.. NOTE::
   The ``cw`` prefix on command is an abbreviation of ``cazy_webscraper``.
   
----------------------
Structure file formats
----------------------

``cw_get_pdb_structures`` can retrieve protein structure files in a series of file formats. The options of file format are (as specified in the BioPython `documentation <https://biopython.org/docs/1.75/api/Bio.PDB.PDBList.html>`_):

* mmCif (default, PDBx/mmCif file),
* pdb (format PDB),
* xml (PDBML/XML format),
* mmtf (highly compressed),
* bundle (PDB formatted archive for large structure}

Any combination of file formats can be provided to ``cw_get_pdb_structures`` to download every file type for each PDB accession in the local CAZyme database. To list multiple file formats, separate each file format with a single comma. For example, to download the mmCif and xml files for every PDB accession in a local CAZyme database (located at ``cazy/cazyme_db.db``), use the following command:

.. code-block:: bash
    
    cw_get_pdb_structures cazy/cazyme_db.db mmCif,xml

.. WARNING::
    The file formats are case sensitive. For example, make sure to use 'mmCif' not 'mmcif'.

--------------------
Command line options
--------------------

``database`` - **REQUIRED** Path to a local CAZyme database to add UniProt data to.

``pdb`` - **REQUIRED** List of file formats to retrieve from PDB for each PDB accession.

``--batch_size`` - Size of an individual batch query of PDB accessions submitted to PDB. Default is 150.

``--cache_dir`` - Path to cache dir to be used instead of default cache dir path.

``--cazy_synonyms`` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

``--config``, ``-c`` - Path to a configuration YAML file. Default: None.

``--classes`` - list of classes to retrieve UniProt data for.

``--ec_filter`` - List of EC numbers to limit the retrieval of structure files to proteins with at least one of the given EC numbers *in the local CAZyme database*.

``--families`` - List of CAZy (sub)families to retrieve UniProt protein data for.

``--genera`` - List of genera to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given genera.

``--log``, ``-l`` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

``--nodelete``, ``-n`` - When called, content in the existing output  will **not** be deleted. Default: False (existing content is deleted).

``--nodelete_cache`` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

``--outdir``, ``-o`` - Define output directory to write out structure files. Default, write structure files to current working directory.

``--retries``, ``-r`` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

``--sql_echo`` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

``--species`` - List of species (organsim scientific names) to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given species.

``--strains`` - List of species strains to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given strains.

``--timeout``, ``-t`` - Connection timout limit (seconds). Default: 45.

``--update_seq`` - If a newer version of the protein sequence is available, overwrite the existing sequence for the protein in the database. Default is false, the protein sequence is **not** overwritten and updated.

``--verbose``, ``-v`` - Enable verbose logging. This does **not** set the SQLite engine ``echo`` parameter to True. Default: False.



-----------
Basic Usage
-----------

The command-line options listed above can be used in combination to customise the retrieval of protein structure files from PDB for proteins of interest. Some options (e.g. ``--families`` and ``--classes``) define the broad group of proteins for which structure files are retrieved, others (e.g. ``--species``) are used to filter and fine-tune the protein dataset for which structure files are retrieved.

The ``--classes``, ``--families``, ``--kingdoms``, ``--genera``, ``--species``, and ``--strains`` filteres are applied 
in the exactly same for retrieving data from CAZy and UniProt, as retrieving data from PDB. Examples of using these flags 
can be found in the ``cazy_webscraper`` and ``cw_get_uniprot_data`` tutorial in this documentation.

.. NOTE::
    To retrieve protein structures for members of specific CAZy subfamilies, list the subfamilies after the ``--families`` 
    flag.


---------------------------------
Structure file retrieval from PDB
---------------------------------

The command for using ``cazy_webscraper`` for retrieval of PDB structure files from PDB is ``cw_get_pdb_structures``.
