==================================
Migrating from version 2 to 3
==================================

Version 3 of ``cazy_webscraper`` reorganises the command-line interface and changes the schema of
the local CAZyme database. This page lists everything that changed, and what you need to do about it.

There are two changes that affect every user:

1. The ten separate ``cw_*`` commands have been replaced by **subcommands of a single**
   ``cazy_webscraper`` **command**, with two of them merged into one new ``extract_data``
   subcommand.
2. The database schema has changed, so **CAZyme databases built with version 2 cannot be used with
   version 3**.

Everything else is a smaller, per-subcommand adjustment, covered in
:ref:`argument-changes` below.


.. _database-compatibility:

------------------------------------------
Databases built with version 2
------------------------------------------

.. warning::
   A local CAZyme database built by version 2 **is not compatible with version 3**, and version 3
   will not tell you so. Rebuild the database with ``download_cazy`` instead of reusing it.

In version 2 the central table of the database was ``Genbanks``, keyed on ``genbank_id`` /
``genbank_accession``. In version 3 that table is ``Proteins``, keyed on ``protein_id`` /
``protein_accession``, with an additional ``source`` column recording where the protein record came
from. Every association table was renamed to match (``Genbanks_CazyFamilies`` became
``Proteins_CazyFamilies``, and so on), two tables were added (``GoTerms`` and ``Pfams``, with the
association tables ``Proteins_GoTerms`` and ``Proteins_Pfams``), and the version 2 ``UniprotTaxs``
table is no longer written.

``Proteins_Pfams`` is populated by the ``get_pfams`` subcommand (see :doc:`get_pfams`), and unlike the
other ``Proteins_*`` association tables carries a row per Pfam domain *match location* rather than
per protein-family pair, since a protein can match the same Pfam family more than once at
different positions - it has its own primary key plus ``match_start``/``match_end`` and
``interpro_accession`` columns alongside the usual ``protein_id``/``pfam_id`` pair.

.. danger::
   **Pointing version 3 at a version 2 database fails silently.**

   ``cazy_webscraper`` opens a database with SQLAlchemy's ``create_all``, which adds any missing
   tables rather than checking the schema it was given. Running a version 3 subcommand against a
   version 2 database therefore succeeds, and *adds* the empty ``Proteins``, ``GoTerms`` and
   ``Pfams`` tables alongside your existing ``Genbanks`` tables. The version 2 data is left intact
   but is now orphaned: it is in ``Genbanks``, and every version 3 subcommand reads ``Proteins``,
   which is empty.

   The visible symptom is a subcommand that runs to completion and reports that it retrieved data
   for no proteins. There is no error message.

The recommended migration is simply to rebuild:

.. code-block:: bash

   # keep the old database, do not overwrite it
   mv cazy_webscraper_2024-01-01_12-00-00.db cazy_webscraper_v2_archive.db

   # build a fresh version 3 database
   cazy_webscraper download_cazy user@example.com -o cazy_webscraper_v3.db

Then re-run whichever ``get_*`` subcommands you had previously used to supplement the database.
Rebuilding is usually the better option regardless of compatibility, because it also picks up the
CAZy annotations that have changed since the old database was built.

.. note::
   There is no automated version 2 to version 3 database conversion tool. If you need to preserve
   data that cannot be re-retrieved, the version 2 tables remain readable with any SQLite client,
   and can be exported before you rebuild.


.. _command-mapping:

-----------------------
Command mapping
-----------------------

Version 2 installed ten console scripts. Version 3 installs one, ``cazy_webscraper``, and each of
the old scripts becomes a subcommand of it.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Version 2 command
     - Version 3 equivalent
   * - ``cazy_webscraper <email>``
     - ``cazy_webscraper download_cazy <email>``
   * - ``cw_get_genbank_seqs <db> <email>``
     - ``cazy_webscraper get_ncbi_seqs <db> <email>``
   * - ``cw_get_ncbi_taxs <db> <email>``
     - ``cazy_webscraper get_ncbi_taxs <db> <email>``
   * - ``cw_get_genomics <db> <email>``
     - ``cazy_webscraper get_ncbi_genomes <db> <email>``
   * - ``cw_get_uniprot_data <db>``
     - ``cazy_webscraper get_uniprot_data <db>``
   * - ``cw_get_pdb_structures <db> <formats>``
     - ``cazy_webscraper get_pdb_structures <db> --file_formats <formats>``
   * - ``cw_get_gtdb_taxs <db> <taxs>``
     - ``cazy_webscraper get_gtdb_taxs <db> --taxs <taxs>``
   * - ``cw_query_database <db> <file_types>``
     - ``cazy_webscraper extract_data <db> --file_types csv json``
   * - ``cw_extract_db_seqs <db> <source>``
     - ``cazy_webscraper extract_data <db> --file_types fasta --source <source>``
   * - ``cw_get_db_schema``
     - *Not yet migrated* (see :ref:`not-yet-migrated`)

Note that the subcommand names are not simply the old script names with the ``cw_`` prefix removed:
``cw_get_genbank_seqs`` became ``get_ncbi_seqs`` and ``cw_get_genomics`` became
``get_ncbi_genomes``.

Two version 2 commands were merged into one. ``cw_query_database`` wrote tabular data out as CSV or
JSON, and ``cw_extract_db_seqs`` wrote protein sequences out as FASTA or a BLAST database. Both are
now ``extract_data``, which chooses between them via ``--file_types``.

The positional arguments are unchanged apart from ``get_pdb_structures`` and ``get_gtdb_taxs``,
which are covered below. Where a version 2 command took ``<db> <email>``, the version 3 subcommand
takes ``<db> <email>`` in the same order.


.. _flag-order:

--------------------------------------------
Flag order: global arguments come first
--------------------------------------------

This is the change most likely to trip up an existing script or shell history.

Five arguments belong to ``cazy_webscraper`` itself rather than to any subcommand, and must be
given **before** the subcommand name:

* ``--version``
* ``--citation``
* ``-l`` / ``--log``
* ``--sql_echo``
* ``-v`` / ``--verbose``

In version 2 every script accepted these directly, so they could be placed anywhere:

.. code-block:: bash

   # version 2
   cw_get_ncbi_taxs cazy.db user@example.com --verbose --log run.log

In version 3 the same flags must precede the subcommand:

.. code-block:: bash

   # version 3
   cazy_webscraper -v -l run.log get_ncbi_taxs cazy.db user@example.com

Putting them after the subcommand is an error rather than a silent no-op, so this change surfaces
immediately when you run the command.

Two of these lost their short forms: version 2's ``-V``/``--version`` and ``-C``/``--citation`` are
now ``--version`` and ``--citation`` only. The short ``-v`` is ``--verbose``, as it was in version 2.


.. _argument-changes:

------------------------------------
Argument changes by subcommand
------------------------------------

Arguments not listed here are unchanged. The filtering arguments shared across subcommands
(``--classes``, ``--families``, ``--kingdoms``, ``--genera``, ``--species``, ``--strains``,
``--ec_filter``) and the common utility arguments (``--cache_dir``, ``--cazy_synonyms``,
``-f``/``--force``, ``--nodelete_cache``, ``-r``/``--retries``) keep their version 2 names and
meanings throughout.

A recurring theme: version 2 had a scattering of differently-named update flags
(``--seq_update``, ``--update_taxs``, ``--update_name``, ``--update_genome_lineage``, and so on).
In version 3 each subcommand has a single ``--update`` flag that means "overwrite what is already
stored in the database when the newly retrieved value differs".

``download_cazy``
=================

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Version 2
     - Version 3
   * - ``-V``/``--version``, ``-C``/``--citation``, ``-l``/``--log``, ``--sql_echo``,
       ``-v``/``--verbose``
     - Moved ahead of the subcommand; ``--version`` and ``--citation`` lost their short forms
       (see :ref:`flag-order`)
   * - ``--validate``
     - Removed
   * - ``-n``/``--nodelete``
     - Removed
   * - ``--nodelete_log``
     - Removed

``-o``/``--db_output``, ``-d``/``--database``, ``-s``/``--subfamilies``, ``--cazy_data``,
``--delete_old_relationships``, ``--skip_ncbi_tax``, ``--ncbi_batch_size`` and ``-t``/``--timeout``
are all unchanged.

``get_ncbi_seqs`` (was ``cw_get_genbank_seqs``)
================================================

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Version 2
     - Version 3
   * - ``--seq_update``
     - ``--update``
   * - ``--seq_dict``
     - Removed
   * - ``--seq_file``
     - Removed
   * - ``-F``/``--file_only``
     - Removed; version 3 always writes retrieved sequences to the local database

``get_ncbi_taxs``
=================

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Version 2
     - Version 3
   * - ``--update_gbk``, ``--update_taxs``
     - ``--update``
   * - ``--use_lineage_cache``
     - Removed
   * - ``--use_protein_ids``, ``--use_tax_ids``
     - Removed

``get_ncbi_genomes`` (was ``cw_get_genomics``)
===============================================

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Version 2
     - Version 3
   * - ``-n``/``--nodelete``
     - Removed
   * - ``--timeout``
     - Removed
   * - *(new)*
     - ``--ec_filter``, to restrict retrieval by EC number as the other subcommands do

``--update`` and ``--batch_size`` keep their version 2 names and behaviour.

.. note::
   ``get_ncbi_genomes`` also accepts ``-F``/``--file_only``, but the flag is not read anywhere in
   version 3 and has no effect. Do not rely on it.

``get_uniprot_data``
====================

This subcommand changed the most.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Version 2
     - Version 3
   * - ``--bioservices_batch_size``
     - ``--batch_size``
   * - ``--update_name``, ``--name_update``, ``--update_seq``
     - ``--update``
   * - ``-t``/``--taxonomy``
     - Removed; use ``get_ncbi_taxs`` or ``get_gtdb_taxs`` for taxonomic data
   * - ``--skip_download``
     - Removed
   * - ``--use_uniprot_cache``
     - Removed
   * - ``--delete_old_ecs``, ``--delete_old_ec_relationships``, ``--delete_old_pdbs``,
       ``--delete_old_pdb_relationships``
     - Removed
   * - ``--timeout``
     - Removed
   * - *(new)*
     - ``-g``/``--go``, to retrieve GO (Gene Ontology) terms
   * - *(new)*
     - ``--swissprot``, to restrict retrieval to the reviewed SwissProt subset of UniProt KB

The data-selection flags ``-e``/``--ec``, ``-p``/``--pdb`` and ``-s``/``--sequence`` are unchanged.

``get_pdb_structures``
======================

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Version 2
     - Version 3
   * - ``<get_pdb_structures>`` positional argument for the file formats
     - ``--file_formats``, an option that defaults to ``mmCif``
   * - *(new)*
     - ``--skip_download``, to record only the experimental method and resolution without
       downloading any structure files
   * - *(new)*
     - ``--update``, to overwrite a stored experimental method and resolution when the retrieved
       values differ
   * - *(new)*
     - ``--ec_filter``, ``-r``/``--retries``, ``-t``/``--timeout``

So a version 2 invocation such as:

.. code-block:: bash

   # version 2
   cw_get_pdb_structures cazy.db pdb mmCif --outdir structures/

becomes:

.. code-block:: bash

   # version 3
   cazy_webscraper get_pdb_structures cazy.db --file_formats pdb mmCif --outdir structures/

Because ``--file_formats`` now has a default, ``cazy_webscraper get_pdb_structures cazy.db`` on its
own downloads ``mmCif`` files, where the version 2 equivalent would have refused to run without the
positional argument.

.. note::
   The ``mmtf`` format that version 2 offered is no longer available: RCSB PDB retired MMTF, and it
   is not one of the version 3 ``--file_formats`` choices, which are ``mmCif``, ``pdb``, ``xml`` and
   ``bundle``.

``get_gtdb_taxs``
=================

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Version 2
     - Version 3
   * - ``<taxs>`` positional argument
     - ``--taxs``, an option that defaults to both ``archaea`` and ``bacteria``
   * - ``--archaea_file``, ``--bacteria_file``
     - Removed
   * - ``--update_genome_lineage``
     - ``--update``
   * - *(new)*
     - ``--ec_filter``, ``-n``/``--nodelete``, ``-t``/``--timeout``

So:

.. code-block:: bash

   # version 2
   cw_get_gtdb_taxs cazy.db bacteria

becomes:

.. code-block:: bash

   # version 3
   cazy_webscraper get_gtdb_taxs cazy.db --taxs bacteria

Listing only the domain you need still matters, because each domain is published by GTDB as a
separate release file and the default retrieves both.


``extract_data`` (was ``cw_query_database`` and ``cw_extract_db_seqs``)
=======================================================================

Version 2 had two separate commands for getting data back out of a local database. Version 3 merges
them into ``extract_data``, where the output format is chosen with ``--file_types`` rather than by
which command you run. The accepted values are ``csv``, ``json``, ``fasta``, ``fasta_dir`` and
``blastdb``, and ``--file_types`` defaults to ``csv``.

From ``cw_query_database``:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Version 2
     - Version 3
   * - ``<file_types>`` positional argument (``csv``, ``json``)
     - ``--file_types``, which also accepts the three sequence formats
   * - ``-p``/``--prefix``
     - ``--prefix`` (the short ``-p`` form was dropped)
   * - ``--include``
     - Unchanged, with ``genus`` added to the accepted values

From ``cw_extract_db_seqs``:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Version 2
     - Version 3
   * - ``<source>`` positional argument (``genbank``, ``uniprot``)
     - ``--source``, which defaults to ``genbank``
   * - ``--fasta_file <path>``
     - ``--file_types fasta``, written into ``-o``/``--output_dir``
   * - ``--fasta_dir <path>``
     - ``--file_types fasta_dir``, written into ``-o``/``--output_dir``
   * - ``-b``/``--blastdb <path>``
     - ``--file_types blastdb``, written into ``-o``/``--output_dir``

The three version 2 sequence arguments each took the output path directly. In version 3 they select
an output *format*, and the destination for all of them is ``-o``/``--output_dir``, which defaults
to the current working directory. So:

.. code-block:: bash

   # version 2
   cw_extract_db_seqs cazy.db genbank --fasta_file all_seqs.fasta

   # version 3
   cazy_webscraper extract_data cazy.db --file_types fasta --source genbank -o .

Because both halves are now one subcommand, a single invocation can write several formats at once:

.. code-block:: bash

   cazy_webscraper extract_data cazy.db --file_types csv fasta --include family organism

``--include`` still only affects the ``csv`` and ``json`` output, and ``--source`` still only
affects the ``fasta``, ``fasta_dir`` and ``blastdb`` output.

.. note::
   ``--file_types blastdb`` requires ``makeblastdb`` from BLAST+ to be available on your ``PATH``.


.. _not-yet-migrated:

------------------------------------
Features not yet migrated
------------------------------------

One version 2 command has no version 3 equivalent yet, and is not registered as a subcommand:

* ``cw_get_db_schema`` — printing the schema of a local database.

If you depend on it, stay on version 2 for that part of your workflow. Note that this means keeping
a version 2 database as well, since the command expects the version 2 schema described in
:ref:`database-compatibility`.

.. note::
   The version 2 ``api`` and ``sequence`` pages have been replaced by :doc:`extract_data`, which
   documents the merged ``extract_data`` subcommand.


-----------------------------
Migration checklist
-----------------------------

#. Install version 3, and check that ``cazy_webscraper --version`` reports it.
#. Replace each ``cw_*`` command in your scripts with the matching subcommand
   (:ref:`command-mapping`).
#. Move ``-v``, ``-l``, ``--sql_echo``, ``--version`` and ``--citation`` ahead of the subcommand
   name (:ref:`flag-order`).
#. Rename the arguments listed in :ref:`argument-changes`, most commonly by collapsing a
   version 2 ``--update_*`` or ``--*_update`` flag into ``--update``.
#. Move the ``get_pdb_structures`` and ``get_gtdb_taxs`` positional arguments to ``--file_formats``
   and ``--taxs``.
#. Replace ``cw_query_database`` and ``cw_extract_db_seqs`` with ``extract_data``, selecting the
   output format with ``--file_types``.
#. Rebuild your CAZyme database with ``download_cazy`` rather than reusing a version 2 database
   (:ref:`database-compatibility`), then re-run the ``get_*`` subcommands you rely on.
#. Check whether you depend on any of the commands in :ref:`not-yet-migrated`.
