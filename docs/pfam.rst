========================================
Retrieving Pfam domains from InterPro
========================================



``cazy_webscraper`` can be used to retrieve Pfam domain annotations for proteins in a local CAZyme database from the `InterPro <https://www.ebi.ac.uk/interpro/>`_ REST API.

.. note::
    Pfam data is retrieved by UniProt accession, not GenBank accession. A freshly built CAZyme database only contains NCBI protein accessions, taxonomic kingdoms, source organisms, and CAZy family annotations. Therefore, the ``cazy_webscraper get_uniprot_data`` command must be used to retrieve UniProt accessions **prior** to using the ``cazy_webscraper get_pfams`` command.

.. note::
    InterPro's per-protein Pfam lookup has no bulk/batch endpoint - one accession is queried per request. ``--batch_size`` controls how many accessions are processed and persisted to the database between progress updates, not how many are sent to InterPro in a single request.

-----------
Quick Start
-----------

To retrieve Pfam domain annotations for every protein with a UniProt accession in a local CAZyme database, use the following command structure:

.. code-block:: bash

   cazy_webscraper get_pfams <path to local CAZyme db>

.. note::
   ``get_pfams`` writes each Pfam family (accession, name, type and the InterPro release the
   data came from) to the ``Pfams`` table, and each domain match location to the
   ``Proteins_Pfams`` table, batch by batch as they are retrieved.

   A protein can have more than one match to the same Pfam family, at different positions
   along its sequence, so ``Proteins_Pfams`` stores one row per match location
   (``match_start``/``match_end``), along with the accession of the InterPro entry the Pfam
   family is integrated into (``interpro_accession``), where InterPro provides one.

   By default, only proteins with no Pfam matches yet are queried, so a plain re-run only
   fills in what is missing. Use ``--update`` to re-query proteins that already have matches.

.. NOTE::
   The ``cw`` prefix on command is an abbreviation of ``cazy_webscraper``.

--------------------
Command line options
--------------------

``database`` - **REQUIRED** Path to a local CAZyme database to add Pfam data to.

``--batch_size`` - Number of UniProt accessions to process, and persist to the database, before reporting progress. Default is 50. InterPro is queried one accession at a time regardless of this value.

``--cache_dir`` - Path to cache dir to be used instead of default cache dir path.

``--cazy_synonyms`` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

``--config``, ``-c`` - Path to a configuration YAML file. Default: None.

``--classes`` - list of classes to retrieve Pfam data for.

``--ec_filter`` - List of EC numbers to limit the retrieval of Pfam data to proteins with at least one of the given EC numbers *in the local CAZyme database*.

``--families`` - List of CAZy (sub)families to retrieve Pfam data for.

``--genbank_accessions`` - Path to text file containing a list of GenBank accessions to retrieve Pfam data for. A unique accession per line.

``--genera`` - List of genera to restrict the retrieval of Pfam data to proteins belonging to one of the given genera.

``--kingdoms`` - List of taxonomy kingdoms to retrieve Pfam data for.

``--log``, ``-l`` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

``--nodelete_cache`` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

``--retries``, ``-r`` - Define the number of times to retry making a connection to InterPro if the connection should fail. Default: 10.

``--sql_echo`` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

``--species`` - List of species (organsim scientific names) to restrict the retrieval of Pfam data to proteins belonging to one of the given species.

``--strains`` - List of species strains to restrict the retrieval of Pfam data to proteins belonging to one of the given strains.

``--timeout``, ``-t`` - Connection timout limit (seconds). Default: 45.

``--uniprot_accessions`` - Path to text file containing a list of UniProt accessions to retrieve Pfam data for. A unique accession per line.

``--update`` - Re-retrieve Pfam data for proteins that already have Pfam matches in the database. Default: False (only proteins with no Pfam matches yet are queried).

``--verbose``, ``-v`` - Enable verbose logging. This does **not** set the SQLite engine ``echo`` parameter to True. Default: False.



-----------
Basic Usage
-----------

The command-line options listed above can be used in combination to customise the retrieval of Pfam domain annotations from InterPro for proteins of interest. Some options (e.g. ``--families`` and ``--classes``) define the broad group of proteins for which Pfam data is retrieved, others (e.g. ``--species``) are used to filter and fine-tune the protein dataset for which Pfam data is retrieved.

The ``--classes``, ``--families``, ``--kingdoms``, ``--genera``, ``--species``, and ``--strains`` filteres are applied
in exactly the same way as for retrieving data from CAZy and UniProt. Examples of using these flags
can be found in the ``cazy_webscraper`` and ``cazy_webscraper get_uniprot_data`` tutorial in this documentation.

.. NOTE::
    To retrieve Pfam domains for members of specific CAZy subfamilies, list the subfamilies after the ``--families``
    flag.


-------------------------------------
Pfam domain retrieval from InterPro
-------------------------------------

The command for using ``cazy_webscraper`` for retrieval of Pfam domain annotations from InterPro is ``cazy_webscraper get_pfams``.
