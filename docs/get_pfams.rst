======================================
Retrieving Pfam domains from InterPro
======================================

The ``get_pfams`` subcommand retrieves Pfam domain annotations for proteins in a local CAZyme
database from the `InterPro <https://www.ebi.ac.uk/interpro/>`_ REST API.

.. note::
   Pfam data is retrieved **by UniProt accession**, not GenBank accession. A freshly built CAZyme
   database contains only NCBI protein accessions, taxonomic kingdoms, source organisms and CAZy
   family annotations, so you must run ``get_uniprot_data`` to retrieve UniProt accessions
   **before** running ``get_pfams``.

-----------
Quick start
-----------

Retrieve Pfam domain annotations for every protein with a UniProt accession:

.. code-block:: bash

   cazy_webscraper get_pfams cazy/cazyme.db

------------------------
How the data is stored
------------------------

``get_pfams`` writes each Pfam family -- accession, name, type, and the InterPro release the data
came from -- to the ``Pfams`` table, and each domain match location to the ``Proteins_Pfams``
table, batch by batch as they are retrieved.

A protein can match the same Pfam family more than once at different positions along its sequence,
so ``Proteins_Pfams`` stores one row per match location (``match_start`` / ``match_end``), together
with the accession of the InterPro entry the Pfam family is integrated into
(``interpro_accession``) where InterPro provides one.

By default only proteins with no Pfam matches yet are queried, so a plain re-run just fills in what
is missing.

-------------------
Subcommand options
-------------------

``database``
  **Required.** Path to the local CAZyme database to add Pfam data to.

``--update``
  Re-query proteins that already have Pfam matches in the database, rather than only those with
  none. Default: ``False``.

``--batch_size``
  Number of UniProt accessions to process, and persist to the database, before reporting progress.
  Default: ``50``.

  .. note::
     InterPro's per-protein Pfam lookup has no bulk endpoint: one accession is queried per request.
     ``--batch_size`` therefore controls how often results are persisted and progress is reported,
     **not** how many accessions are sent to InterPro at once.

``-t``, ``--timeout``
  Seconds before a connection to InterPro times out. Default: ``45``.

For the filtering flags and the shared housekeeping options, see
:doc:`Filtering and common options <filters>`. ``get_pfams`` accepts both
``--genbank_accessions`` and ``--uniprot_accessions``.

-----------------
Worked examples
-----------------

**Retrieve Pfam annotations for two CAZy classes:**

.. code-block:: bash

   cazy_webscraper get_pfams cazy/cazyme.db --classes GH,CE

**Retrieve annotations for a functionally defined subset:**

.. code-block:: bash

   cazy_webscraper get_pfams cazy/cazyme.db \
       --classes GH --kingdoms bacteria --ec_filter "3.2.1.23,3.2.1.37"

**Refresh annotations for a named list of proteins:**

.. code-block:: bash

   cazy_webscraper get_pfams cazy/cazyme.db --update \
       --uniprot_accessions my_proteins.txt

See :ref:`filter-combining` for more on combining filters.
