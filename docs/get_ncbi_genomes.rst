=====================================================
Retrieving genomic assembly data from NCBI Assembly
=====================================================

The ``get_ncbi_genomes`` subcommand retrieves the genomic assembly each protein came from, and adds
the assembly name, GenBank and RefSeq version accessions, and assembly ID to a local CAZyme
database.

-----------
Quick start
-----------

.. code-block:: bash

   cazy_webscraper get_ncbi_genomes cazy/cazyme.db user@example.com

.. note::
   The email address is a requirement of NCBI Entrez. It is not stored by ``cazy_webscraper``. See
   the `NCBI Entrez <https://www.ncbi.nlm.nih.gov/books/NBK25497/>`_ documentation.

-------------------
Subcommand options
-------------------

``database``
  **Required.** Path to the local CAZyme database.

``email``
  **Required.** Your email address, for NCBI Entrez.

``--update``
  Re-check every protein and overwrite the assembly data already stored against it. Without this
  flag, only proteins with no assembly data are looked up. Default: ``False``.

  .. warning::
     ``--update`` overwrites existing assembly records in the database.

``--batch_size``
  Number of records per batch query sent to NCBI Entrez. Default: ``150``.

.. warning::
   ``get_ncbi_genomes`` also accepts ``-F``/``--file_only``, but nothing in version 3 reads it and
   it has no effect. Do not rely on it.

For the filtering flags and the shared housekeeping options, see
:doc:`Filtering and common options <filters>`.

-----------------
Worked examples
-----------------

**Retrieve assembly data for bacterial CAZymes in two classes:**

.. code-block:: bash

   cazy_webscraper get_ncbi_genomes cazy/cazyme.db user@example.com \
       --classes GH,CE --kingdoms bacteria

**Refresh assembly data for a named list of proteins:**

.. code-block:: bash

   cazy_webscraper get_ncbi_genomes cazy/cazyme.db user@example.com --update \
       --genbank_accessions my_proteins.txt

See :ref:`filter-combining` for more on combining filters.
