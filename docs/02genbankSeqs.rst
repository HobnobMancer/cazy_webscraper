==========================================
Retrieving protein sequences from NCBI
==========================================

The ``get_ncbi_seqs`` subcommand retrieves protein amino acid sequences from NCBI GenBank for
proteins in a local CAZyme database, and stores them in the database along with the sequence's
last-updated date.

-----------
Quick start
-----------

.. code-block:: bash

   cazy_webscraper get_ncbi_seqs cazy/cazyme.db user@example.com

.. note::
   The email address is a requirement of NCBI Entrez. It is not stored by ``cazy_webscraper``. See
   the `NCBI Entrez <https://www.ncbi.nlm.nih.gov/books/NBK25497/>`_ documentation.

.. warning::
   Please do not retrieve sequences for the entire CAZy dataset unless it is genuinely necessary.
   Retrieve the subset you need, and run large jobs during quiet periods.

-------------------
Subcommand options
-------------------

``database``
  **Required.** Path to the local CAZyme database.

``email``
  **Required.** Your email address, for NCBI Entrez.

``--update``
  Overwrite a sequence already stored in the database when the sequence retrieved from NCBI
  differs. Without this flag, ``get_ncbi_seqs`` only fills in blanks: a protein that already has a
  sequence keeps it. Default: ``False``.

  .. code-block:: bash

     cazy_webscraper get_ncbi_seqs cazy/cazyme.db user@example.com --update

``--batch_size``
  Number of records per batch query sent to NCBI Entrez. Default: ``150``.

For the filtering flags and the shared housekeeping options, see
:doc:`Filtering and common options <filters>`.

.. note::
   ``get_ncbi_seqs`` accepts ``--genbank_accessions`` but not ``--uniprot_accessions``.

-----------------
Worked examples
-----------------

**Retrieve sequences for bacterial CAZymes in two classes:**

.. code-block:: bash

   cazy_webscraper get_ncbi_seqs cazy/cazyme.db user@example.com \
       --classes GH,CE --kingdoms bacteria

**Retrieve sequences for a functionally defined subset:**

.. code-block:: bash

   cazy_webscraper get_ncbi_seqs cazy/cazyme.db user@example.com \
       --classes GH --ec_filter "3.2.1.23,3.2.1.37,3.2.1.85"

**Refresh sequences for a named list of proteins:**

.. code-block:: bash

   cazy_webscraper get_ncbi_seqs cazy/cazyme.db user@example.com --update \
       --genbank_accessions my_proteins.txt

See :ref:`filter-combining` for more on combining filters.
