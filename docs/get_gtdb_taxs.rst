==========================================
Retrieving GTDB taxonomic classifications
==========================================

The ``get_gtdb_taxs`` subcommand adds taxonomic classifications from the
`Genome Taxonomy Database <https://gtdb.ecogenomic.org/>`_ to the genomes in a local CAZyme
database.

.. note::
   GTDB classifies *genomes*, so ``get_gtdb_taxs`` needs genomic assembly accessions to be present
   already. Run ``get_ncbi_genomes`` first.

-----------
Quick start
-----------

.. code-block:: bash

   cazy_webscraper get_gtdb_taxs cazy/cazyme.db

By default this retrieves both the archaeal and bacterial GTDB releases.

------------------------
Selecting GTDB domains
------------------------

``--taxs`` chooses which GTDB domains to retrieve classifications from. It accepts ``archaea``,
``bacteria`` or both, and defaults to both.

.. code-block:: bash

   cazy_webscraper get_gtdb_taxs cazy/cazyme.db --taxs bacteria

.. tip::
   Each domain is published by GTDB as a separate release file, so naming only the domain you need
   avoids downloading the other one. If your dataset is entirely bacterial, ``--taxs bacteria``
   halves the download.

-------------------
Subcommand options
-------------------

``database``
  **Required.** Path to the local CAZyme database.

``--taxs``
  GTDB domains to retrieve. ``archaea``, ``bacteria``, or both. Default: both.

``--update``
  Re-check every genome and overwrite the GTDB classification already stored against it. Without
  this flag, only genomes with no classification yet are looked up. Default: ``False``.

``-n``, ``--nodelete``
  Keep the existing contents of the output directory. Default: ``False``.

``-t``, ``--timeout``
  Seconds before a connection to GTDB times out. Default: ``45``.

For the filtering flags and the shared housekeeping options, see
:doc:`Filtering and common options <filters>`. ``get_gtdb_taxs`` accepts both
``--genbank_accessions`` and ``--uniprot_accessions``; classifications are retrieved for the
genomes of the proteins listed.

-----------------
Worked examples
-----------------

**Retrieve bacterial classifications for two CAZy classes:**

.. code-block:: bash

   cazy_webscraper get_gtdb_taxs cazy/cazyme.db --taxs bacteria --classes GH,CE

**Refresh every stored classification:**

.. code-block:: bash

   cazy_webscraper get_gtdb_taxs cazy/cazyme.db --update

**Classify the genomes behind a named list of proteins:**

.. code-block:: bash

   cazy_webscraper get_gtdb_taxs cazy/cazyme.db --taxs bacteria \
       --genbank_accessions my_proteins.txt

See :ref:`filter-combining` for more on combining filters.
