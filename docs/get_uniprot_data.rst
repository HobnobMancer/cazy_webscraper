============================
Retrieving data from UniProt
============================

The ``get_uniprot_data`` subcommand supplements an existing local CAZyme database with data from
`UniProt <https://www.uniprot.org/>`_.

.. note::
   ``get_uniprot_data`` only adds data to proteins already in the local CAZyme database. It never
   adds new proteins. Build the database with ``download_cazy`` first.

-----------
Quick start
-----------

Retrieve UniProt accessions and protein names for every protein in a local CAZyme database:

.. code-block:: bash

   cazy_webscraper get_uniprot_data cazy/cazyme.db

That is the default behaviour: all proteins, accessions and names only. EC numbers, GO terms, PDB
accessions and sequences are each opted into with a flag:

.. code-block:: bash

   cazy_webscraper get_uniprot_data cazy/cazyme.db --ec --pdb --sequence

.. warning::
   Please do not retrieve data for the entire CAZy dataset unless it is genuinely necessary. Doing
   so takes several hours and may deny the service to others.

   For large datasets, break the work into subgroups and retrieve one at a time, and run during
   quiet periods where you can. Very large single retrievals (say, over 1,000,000 proteins, or all
   of GH) can cause UniProt to terminate the connection for sustained high bandwidth use.

---------------------
What can be retrieved
---------------------

Each of these flags turns on one category of data. None are on by default.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Flag
     - Retrieves
   * - ``-e``, ``--ec``
     - EC number annotations
   * - ``-g``, ``--go``
     - GO (Gene Ontology) terms
   * - ``-p``, ``--pdb``
     - PDB accessions
   * - ``-s``, ``--sequence``
     - Protein amino acid sequences, stored with UniProt's 'last modified' date

.. note::
   ``--ec`` retrieves **all** EC annotations for every protein matching your filters. It is not
   restricted by ``--ec_filter``, which selects *which proteins* to process rather than which
   annotations to keep. Used together, ``--ec_filter`` picks the protein set and ``--ec`` then
   retrieves every EC number those proteins have.

-------------------------
Subcommand options
-------------------------

``database``
  **Required.** Path to the local CAZyme database to add UniProt data to.

``--swissprot``
  Restrict retrieval to proteins in SwissProt, the reviewed subset of UniProt KB.
  Default: ``False``.

``--update``
  Overwrite existing values in the local database when the retrieved value differs. Without it,
  ``get_uniprot_data`` only fills in blanks: a protein that already has a sequence keeps it.
  Default: ``False``.

  .. code-block:: bash

     cazy_webscraper get_uniprot_data cazy/cazyme.db -s --update

  For sequences, "differs" is decided by comparing the local 'last modified' date against the one
  retrieved from UniProt, so an update only happens when UniProt genuinely has a newer record.

``--batch_size``
  Number of records per batch query submitted to the `UniProt REST API
  <https://www.uniprot.org/help/programmatic_access>`_. Default: ``150``.

For the filtering flags and the shared housekeeping options, see
:doc:`Filtering and common options <filters>`.

^^^^^^^^^^^^^^^^^^^^^^^^^
UniProt batch size limits
^^^^^^^^^^^^^^^^^^^^^^^^^

UniProt applies its own limits to ID mapping job submissions, which constrain how large
``--batch_size`` can usefully be:

========= =====================================================================================
Limit	  Details
========= =====================================================================================
100,000	  Total number of ids allowed in comma separated param ids in /idmapping/run api
500,000	  Total number of "mapped to" ids allowed
100,000	  Total number of "mapped to" ids allowed to be enriched by UniProt data
10,000	  Total number of "mapped to" ids allowed with filtering
========= =====================================================================================

-----------------
Worked examples
-----------------

**Retrieve PDB accessions for bacterial CAZymes in selected classes and families:**

.. code-block:: bash

   cazy_webscraper get_uniprot_data cazy/cazyme.db --pdb \
       --classes GH,CE --families CE1,CE5,CE8 --kingdoms bacteria

**Retrieve EC numbers and GO terms for two genera.** Several data flags can be combined in one
run, which is cheaper than repeating the query per data type:

.. code-block:: bash

   cazy_webscraper get_uniprot_data cazy/cazyme.db --ec --go \
       --classes GH --genera Aspergillus,Trichoderma

**Refresh sequences for a functionally defined subset.** ``--ec_filter`` narrows the protein set
using EC annotations already in the database, and ``--update`` refreshes sequences already stored:

.. code-block:: bash

   cazy_webscraper get_uniprot_data cazy/cazyme.db --sequence --update \
       --classes GH,CE,CBM --kingdoms bacteria \
       --ec_filter "3.2.1.23,3.2.1.37,3.2.1.85"

**Retrieve reviewed data only, for a named list of proteins:**

.. code-block:: bash

   cazy_webscraper get_uniprot_data cazy/cazyme.db --ec --swissprot \
       --genbank_accessions my_proteins.txt

See :ref:`filter-combining` for more on combining filters.
