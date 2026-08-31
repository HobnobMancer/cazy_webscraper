============================================
Retrieving taxonomies from NCBI Taxonomy
============================================

The ``get_ncbi_taxs`` subcommand retrieves the latest taxonomic classifications from NCBI Taxonomy
for proteins in a local CAZyme database, including the NCBI Taxonomy ID and the full lineage.

-----------
Quick start
-----------

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs cazy/cazyme.db user@example.com

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
  Overwrite existing, outdated taxonomies in the local CAZyme database. Without this flag, only
  proteins with no NCBI taxonomy yet are looked up. Default: ``False``.

``--batch_size``
  Number of records per batch query sent to NCBI Entrez. Default: ``150``.

For the filtering flags and the shared housekeeping options, see
:doc:`Filtering and common options <filters>`. ``get_ncbi_taxs`` accepts both
``--genbank_accessions`` and ``--uniprot_accessions``.

-----------------
Worked examples
-----------------

**Retrieve taxonomies for two CAZy classes:**

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs cazy/cazyme.db user@example.com --classes GH,CE

**Refresh every stored taxonomy:**

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs cazy/cazyme.db user@example.com --update

**Retrieve taxonomies for a named list of proteins:**

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs cazy/cazyme.db user@example.com \
       --uniprot_accessions my_proteins.txt

See :ref:`filter-combining` for more on combining filters.
