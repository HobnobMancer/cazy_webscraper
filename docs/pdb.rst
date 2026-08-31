===================================
Retrieving structure files from PDB
===================================

The ``get_pdb_structures`` subcommand downloads protein structure files from
`RCSB PDB <https://www.rcsb.org/>`_ for proteins in a local CAZyme database, and records the
experimental method and resolution of each structure in the database.

.. note::
   Structure files are written to disk. Only the experimental method and resolution are stored in
   the local CAZyme database.

   ``get_pdb_structures`` works from PDB accessions already in your database, so run
   ``get_uniprot_data --pdb`` first to populate them.

-----------
Quick start
-----------

Download an ``mmCif`` file for every protein in the database that has a PDB accession:

.. code-block:: bash

   cazy_webscraper get_pdb_structures cazy/cazyme.db

To record structure metadata without downloading any files:

.. code-block:: bash

   cazy_webscraper get_pdb_structures cazy/cazyme.db --skip_download

------------------------
Structure file formats
------------------------

``--file_formats`` selects which formats to download. It accepts one or more of ``mmCif``, ``pdb``,
``xml`` and ``bundle``, separated by spaces, and defaults to ``mmCif``.

.. code-block:: bash

   cazy_webscraper get_pdb_structures cazy/cazyme.db --file_formats pdb mmCif

.. note::
   ``bundle`` is only published for very large structures, so requesting it alone will retrieve
   nothing for most proteins.

.. warning::
   The ``mmtf`` format offered by version 2 is gone: RCSB PDB has retired MMTF, and it is not an
   accepted value.

-------------------
Subcommand options
-------------------

``database``
  **Required.** Path to the local CAZyme database.

``--file_formats``
  Structure file format(s) to download. Default: ``mmCif``.

``--skip_download``
  Do not download structure files. Only retrieve the experimental method and resolution of each
  structure and add them to the local CAZyme database. Default: ``False``.

``--update``
  Overwrite the experimental method and resolution already stored in the database when the values
  retrieved from PDB differ. Default: ``False``.

``-o``, ``--outdir``
  Directory to write structure files to. Default: the current working directory.

``--overwrite``
  Overwrite existing structure files. Default: ``False``.

``-n``, ``--nodelete``
  Keep the existing contents of the output directory. Default: ``False``.

``--batch_size``
  Number of records per batch request to PDB. Default: ``150``.

``-t``, ``--timeout``
  Seconds before a connection to PDB times out. Default: ``45``.

For the filtering flags and the shared housekeeping options, see
:doc:`Filtering and common options <filters>`. ``get_pdb_structures`` accepts both
``--genbank_accessions`` and ``--uniprot_accessions``.

-----------------
Worked examples
-----------------

**Download two formats for bacterial CAZymes in selected classes:**

.. code-block:: bash

   cazy_webscraper get_pdb_structures cazy/cazyme.db --file_formats pdb mmCif \
       --classes GH,CE --kingdoms bacteria --outdir structures/

**Record metadata only for a functionally defined subset,** without downloading files:

.. code-block:: bash

   cazy_webscraper get_pdb_structures cazy/cazyme.db --skip_download \
       --classes GH --ec_filter "3.2.1.23,3.2.1.37"

**Refresh stored metadata for a named list of proteins:**

.. code-block:: bash

   cazy_webscraper get_pdb_structures cazy/cazyme.db --skip_download --update \
       --uniprot_accessions my_proteins.txt

See :ref:`filter-combining` for more on combining filters.
