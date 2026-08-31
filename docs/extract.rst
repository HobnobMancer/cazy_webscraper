============================================
Extracting data from a local CAZyme database
============================================

The ``extract_data`` subcommand writes data held in a local CAZyme database out to files. It
covers both tabular export (CSV, TSV, JSON, JSON Lines) and sequence export (FASTA files and BLAST
databases).

.. note::
   In version 2 these were two separate commands, ``cw_query_database`` and ``cw_extract_db_seqs``.
   They are now one subcommand, and the output format is chosen with ``--file_types`` rather than by
   which command you run. See :doc:`migration`.

The database can of course also be queried directly with SQL, through the ``sqlite3`` command line
or a browser such as `DB Browser for SQLite <https://sqlitebrowser.org/>`_. ``extract_data`` exists
to make the dataset easy to feed into downstream pipelines, and to make it usable without SQL.

-----------
Quick start
-----------

Write every protein in the database to a CSV file in the current directory:

.. code-block:: bash

   cazy_webscraper extract_data cazy/cazyme.db

Write all GenBank sequences to a single FASTA file instead:

.. code-block:: bash

   cazy_webscraper extract_data cazy/cazyme.db --file_types fasta

-------------------
Output file types
-------------------

``--file_types`` selects one or more output formats, separated by spaces. It defaults to ``csv``.

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Value
     - Output
   * - ``csv``
     - One comma-separated file of all matching proteins.
   * - ``tsv``
     - The same table, tab-separated.
   * - ``json``
     - One JSON file containing all matching proteins.
   * - ``jsonl``
     - JSON Lines: one JSON object per line, convenient for streaming into other tools.
   * - ``fasta``
     - One FASTA file containing all extracted sequences.
   * - ``fasta_dir``
     - One FASTA file per protein.
   * - ``blastdb``
     - A BLAST protein database.

Several formats can be produced in a single run, which is cheaper than querying the database
repeatedly:

.. code-block:: bash

   cazy_webscraper extract_data cazy/cazyme.db --file_types csv fasta blastdb

.. note::
   ``--file_types blastdb`` requires ``makeblastdb`` from BLAST+ to be available on your ``PATH``.

.. note::
   Output files are written into ``-o``/``--output_dir``, which defaults to the current working
   directory. Any parent directories that do not exist are created for you.

   In version 2, ``--fasta_file``, ``--fasta_dir`` and ``--blastdb`` each took a path directly.
   They now select a *format*, and the destination for all of them is ``--output_dir``.

^^^^^^^^^^^^^^^^^^
FASTA file format
^^^^^^^^^^^^^^^^^^

FASTA records written by ``extract_data`` have a deliberately simple ID line, containing only the
accession and the name of the source database:

.. code-block:: text

   >AIP21820.1 GenBank
   MPVALAVAAALGACSGDDDATLESRADAIVERMTTRQKVGQKLMMAFRYWCPDGQPACTT
   GMTEFPDAARDALRENGIGGVILFSNNLTGIEQTRRLIDGIRAAPAADSPLGLMIGIDEE
   GGNVFRLPRVEATAFAGNMALGAAYEATRDDRLAYDMGRVLAAEIAAVGFNVNFAPDVDV

-----------------------------
Choosing the sequence source
-----------------------------

``--source`` decides which protein sequences go into the ``fasta``, ``fasta_dir`` and ``blastdb``
outputs. It accepts ``genbank``, ``uniprot`` or both, and defaults to ``genbank``.

.. code-block:: bash

   cazy_webscraper extract_data cazy/cazyme.db --file_types fasta --source genbank uniprot

.. tip::
   The order does not matter and the values are not case sensitive.

.. note::
   ``--source`` has no effect on the ``csv``, ``tsv``, ``json`` and ``jsonl`` outputs. To include
   sequences in those, use ``--include genbank_seq`` and/or ``--include uniprot_seq``.

--------------------------
Customising the output
--------------------------

By default the tabular outputs contain only the GenBank accession of each matching protein.
``--include`` adds further columns, named with spaces between them:

.. code-block:: bash

   cazy_webscraper extract_data cazy/cazyme.db --include family organism ec

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Value
     - Adds
   * - ``class``
     - CAZy class annotations
   * - ``family``
     - CAZy family annotations
   * - ``subfamily``
     - CAZy subfamily annotations
   * - ``kingdom``
     - Taxonomic kingdom of the source organism
   * - ``genus``
     - Genus of the source organism
   * - ``organism``
     - Scientific name of the source organism
   * - ``ec``
     - EC number annotations
   * - ``pdb``
     - PDB accessions
   * - ``pfam``
     - Pfam domain accessions
   * - ``uniprot_acc``
     - UniProt accession
   * - ``uniprot_name``
     - Protein name retrieved from UniProt
   * - ``genbank_seq``
     - Protein sequence retrieved from GenBank
   * - ``uniprot_seq``
     - Protein sequence retrieved from UniProt

.. note::
   Columns always appear in the fixed order listed above, whatever order you name them in. Both
   ``--include kingdom genus organism`` and ``--include organism genus kingdom`` produce the
   columns Kingdom, Genus, Organism, in that order.

.. note::
   ``--include`` has no effect on the ``fasta``, ``fasta_dir`` and ``blastdb`` outputs, whose
   content is governed by ``--source``.

-------------------
Subcommand options
-------------------

``database``
  **Required.** Path to the local CAZyme database.

``--file_types``
  Output format(s) to write. Default: ``csv``.

``--include``
  Additional columns for the tabular outputs. Default: ``None``.

``--source``
  Which sequences to write to the sequence outputs. Default: ``genbank``.

``-o``, ``--output_dir``
  Directory to write output files to. Default: the current working directory.

``--prefix``
  String to prefix every output file name with. Default: ``None``.

``--overwrite``
  Overwrite existing output files. Default: ``False``.

``-n``, ``--nodelete``
  Keep the existing contents of the output directory. Default: ``False``.

``--batch_size``
  Number of proteins assembled and written at a time. Default: ``1000``.

``-t``, ``--timeout``
  Seconds before a connection times out. Default: ``45``.

For the filtering flags and the shared housekeeping options, see
:doc:`Filtering and common options <filters>`. ``extract_data`` accepts both
``--genbank_accessions`` and ``--uniprot_accessions``.

-----------------
Worked examples
-----------------

**Export bacterial CAZymes from selected classes and families, with family and organism columns:**

.. code-block:: bash

   cazy_webscraper extract_data cazy/cazyme.db --file_types csv \
       --classes GH,CE --families CE1,CE5,CE8 --kingdoms bacteria \
       --include family organism

**Export a table and a BLAST database in one run,** for two genera:

.. code-block:: bash

   cazy_webscraper extract_data cazy/cazyme.db --file_types csv blastdb \
       --classes GH --genera Aspergillus,Trichoderma \
       --include class ec pdb --output_dir results/

**Export sequences from both sources for a functionally defined subset:**

.. code-block:: bash

   cazy_webscraper extract_data cazy/cazyme.db --file_types fasta \
       --source genbank uniprot \
       --classes GH,CE,CBM --ec_filter "3.2.1.23,3.2.1.37,3.2.1.85"

**Export a named list of proteins with every available column:**

.. code-block:: bash

   cazy_webscraper extract_data cazy/cazyme.db --file_types json \
       --genbank_accessions my_proteins.txt \
       --include class family subfamily kingdom genus organism ec pdb pfam \
                 uniprot_acc uniprot_name

See :ref:`filter-combining` for more on combining filters.
