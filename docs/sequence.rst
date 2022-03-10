================================================
Extract protein squences from the local database
================================================

``cazy_webscraper`` can be used to extract protein sequences stored in the local CAZyme database.

The extracted protein sequences can be written to any combination of:
* A single FASTA file containing *all* extracted sequences
* One FASTA file per extracted sequence
* A BLASTp database

GenBank and/or UniProt protein sequences retrieved can be extracted.

-----------
Quick Start
-----------

To extract all protein sequences previously retrieved from GenBank and UniProt, use the following command structure:

.. code-block:: bash

   cw_extract_db_seqs <path to local CAZyme db> genbank uniprot

.. NOTE::
   The ``cw`` prefix on command is an abbreviation of ``cazy_webscraper``.

.. NOTE::
    'genbank' and 'uniprot' are not case sensitive, therefore, both GenBank and UniProt are also 
    accepted.

.. WARNING::
    At least one database (either GenBank or UniProt) **must** be provided.

--------------------
Command line options
--------------------

``database`` - **REQUIRED** Path to a local CAZyme database to add UniProt data to.

``source`` - **REQUIRED** Define source databases of protein sequences. Accepts 'genbank' and 'uniprot'. To list both, separate with a single space, e.g.   

.. code-block:: bash
`cw_extract_sequence cazy_database.db genbank uniprot

*The database names are not case sensitive, therefore, both GenBank and genbank are accepted.* 

``--blastdb``, ``-b`` - Create BLAST database of extracted protein sequences. Provide the path to the directory to store the BLAST database in.

``--fasta_dir`` - Write out each extracted sequence to a separate FASTA file in the provided dir. Provide a path to a directory to write out the FASTA files.

``--fasta_file`` - Write out all extracted sequences to a single FASTA file. Provide a path to write out the FASTA file.

.. WARNING::
    At least one of ``--blastdb``, ``--fasta_dir``, and ``--fasta_file`` must be called to inform `cazy_webscraper` where to write the output to. If none are called sequences will be extracted._

``--cache_dir`` - Path to cache dir to be used instead of default cache dir path.

``--cazy_synonyms`` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

``--config``, ``-c`` - Path to a configuration YAML file. Default: None.

``--classes`` - List of classes from which all families are to be scrape.

``--ec_filter`` - Limist retrieval of protein data to proteins annotated with a provided list of EC numbers. Separate the EC numbers bu single commas without spaces. Recommend to wrap the entire str in quotation marks, for example:
.. code-block:: bash
    cw_get_uniprot_data my_cazyme_db/cazyme_db.db --ec_filter 'EC1.2.3.4,EC2.3.1.-'

``--force``, ``-f`` - Force overwriting exsiting files and writing to existing output directory.

``--families``` - List of CAZy (sub)families to scrape.#

``--kingdoms`` - List of taxonomy kingdoms to retrieve UniProt data for.

``--genbank_accession`` - Path to text file containing a list of GenBanks accessions to extract protein sequences for. A unique accession per line.

``--genera`` - List of genera to restrict the scrape to. Default: None, filter not applied to scrape.

``--log``, ``-l`` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

``--nodelete`` - When called, content in the existing output dir will **not** be deleted. Default: False (existing content is deleted).

``--nodelete_cache`` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

``--sql_echo`` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

``--species`` - List of species written as Genus Species) to restrict the scraping of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

``--strains`` - List of specific species strains to restrict the scraping of CAZymes to.

``--uniprot_accessions`` - Path to text file containing a list of UniProt accessions to extract protein sequences for. A unique accession per line.

``--verbose``, ``-v`` - Enable verbose logging. This does not set the SQLite engine `echo` parameter to True. Default: False.

-----------
Basic Usage
-----------

The command-line options listed above can be used in combination to customise the retrieval the extraction of protein sequences 
to proteins of interest. Some options (e.g. ``--families`` and ``--classes``) define the broad group of proteins, 
others (e.g. ``--species``) are used to filter and fine-tune the protein dataset.

The ``--classes``, ``--families``, ``--kingdoms``, ``--genera``, ``--species``, and ``--strains`` filteres are applied 
in the exactly same for retrieving data from CAZy and UniProt. Examples of using these flags 
can be found in the ``cazy_webscraper`` and ``cw_get_uniprot_data`` tutorial in this documentation.

.. NOTE::
    To extract protein sequences for members of specific CAZy subfamilies, list the subfamilies after the ``--families`` 
    flag.
