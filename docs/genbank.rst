=========================================
Retrieving Protein Sequences from GenBank
=========================================

``cazy_webscraper`` can be used to retrieve protein amino acid sequences from NCBI GenBank for user-specified data sets of CAZymes 
in the local CAZymes database. 

The retrieval of data from NCBI is performed by using the **BioPython** `Bio.entrez <https://biopython.org/docs/1.75/api/Bio.Entrez.html>_` module [Cock *et al.*, 2009].

    Cock, P. J. A, Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A. _et al._ (2009) 'Biopython: freely available Python tools for computaitonal molecular biology and bioinformatics', _Bioinformatics_, 25(11), pp. 1422-3.

.. note::
    For specific information of the ``Bio.entrez`` module please see the 
    `entrez documentation <https://biopython.org/docs/1.75/api/Bio.Entrez.html>`_.


-----------
Quick Start
-----------

To download protein sequences for all CAZymes in the local CAZyme database, and write them to the local CAZyme database, 
use the following command structure:

.. code-block:: console

    cw_get_genbank_seqs 'path to local CAZyme db' 'user email address'

For example:


.. code-block:: console
    
    cw_get_genbank_seqs cazy/cazyme.db myemail@domain.com

.. NOTE::
   The ``cw`` prefix is an abbreviation of ``cazy_webscraper``.


--------------------
Command line options
--------------------

``database`` - **REQUIRED** Path to a local CAZyme database to add UniProt data to.

``email`` - **REQUIRED** User email address, required by NCBI Entrez.

``--cache_dir`` - Path to cache dir to be used instead of default cache dir path.

``--cazy_data`` - Path to a txt file downloaded from CAZy containing a CAZy database dump

``--cazy_synonyms`` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

``--config``, ``-c`` - Path to a configuration YAML file. Default: None.

``--classes`` - list of classes to retrieve UniProt data for.

``--ec_filter`` - List of EC numbers to limit the retrieval of protein data for proteins annotated with at least one of the given EC numbers **in the local CAZyme database**.

``--entrez_batch_size`` - Change the query batch size submitted via [`Entrez`]() to retrieve protein sequences from GenBank data. Default is 150. `Entrez <https://www.ncbi.nlm.nih.gov/books/NBK179288/>_` recommands queries not larger than XXX objects in length.

``--families`` - List of CAZy (sub)families to retrieve UniProt protein data for.

``--genbank_accessions`` - Path to text file containing a list of GenBank accessions to retrieve protein data for. A unique accession per line.

``--genera`` - List of genera to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given genera.

``--log``, ``-l`` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

``--nodelete_cache`` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

``--nodelete_log`` - When called, content in the existing log dir will **not** be deleted. Default: False (existing content is deleted).

``--retries``, ``-r`` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

``--seq_dict``, - Path to a JSON file, keyed by GenBank accessions and valued by protein sequence. This file is created as part of the cache, after all protein sequences are retrieved from GenBank. This skips the retrieval of the protein sequences from GenBank.

``--seq_update`` - If a newer version of the protein sequence is available, overwrite the existing sequence for the protein in the database. Default is false, the protein sequence is **not** overwritten and updated.

``--sql_echo`` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

``--species`` - List of species (organsim scientific names) to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given species.

``--strains`` - List of species strains to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given strains.

``--timeout``, ``-t`` - Connection timout limit (seconds). Default: 45.

``--verbose``, ``-v`` - Enable verbose logging. This does **not** set the SQLite engine ``echo`` parameter to True. Default: False.

-----------
Basic Usage
-----------

The command-line options listed above can be used in combination to customise the scraping of CAZy. Some options (e.g. ``--families`` and ``--classes``) define the broad group of data that will be scraped, others (e.g. ``--species``) are used to filter and fine-tune the data that is scraped.

The ``--classes``, ``--families``, ``--kingdoms``, ``--genera``, ``--species``, and ``--strains`` filteres are applied 
in the exactly same for retrieving data from CAZy as retrieving data from UniProt. Examples of using these flags 
can be found in the ``cazy_webscraper`` tutorial in this documentation.

The ``--seq_update`` flag is used in the same way for retrieving protein sequences from UniProt and GenBank.

.. NOTE::
    To retrieve data for members of specific CAZy subfamilies, list the subfamilies after the ``--families`` 
    flag.

------------------------
Updating local sequences
------------------------

When using ``--sequence`` flag, ``cazy_webscraper`` will only add *new* protein sequences to the database, i.e.
it will only add protein sequences to records that do not have a sequence. Therefore, if a protein
already has a sequence in the local database, this sequence is **not** overwritten.

You may wish to update the protein sequences in your local CAZyme database. To do this use the ``--sequence``/``-s`` 
flag to tell ``cazy_webscraper`` to retrieve protein sequences, **and** use the ``--seq_update`` flag.

.. code-block:: console

    cw_get_genbank_seqs cazy_db.db -s --seq_update

This instructs ``cazy_webscraper`` to overwriting existing protein sequences in the local database *if* a newer version 
of the sequence is retrieved from UniProt. This is checked by comparing the 'last modified date' of the 
protein sequence in the local database against the sequence retrieved from UniProt.
