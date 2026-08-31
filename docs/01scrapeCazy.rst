=====================================
Creating a local CAZyme database
=====================================

The ``download_cazy`` subcommand retrieves data from `CAZy <http://www.cazy.org/>`_ and compiles it
into a local SQLite3 database. It is the only subcommand that *builds* a database; every ``get_*``
subcommand adds data to one that already exists.

Scraping the entirety of CAZy typically takes 5-15 minutes.

-----------
Quick start
-----------

Download all of CAZy into a new database in the current directory, named
``cazy_webscraper_<date>_<time>.db``:

.. code-block:: bash

   cazy_webscraper download_cazy user@example.com

.. note::
   The email address is a requirement of NCBI Entrez, which is queried to resolve the source
   organism when CAZy lists more than one for a single protein. It is not stored by
   ``cazy_webscraper``. See the `NCBI Entrez
   <https://www.ncbi.nlm.nih.gov/books/NBK25497/>`_ documentation.

-------------------------------
Choosing where data is written
-------------------------------

``-o``/``--db_output`` builds a **new** database at the path you give:

.. code-block:: bash

   cazy_webscraper download_cazy user@example.com -o cazy/cazyme.db

``-d``/``--database`` adds data to an **existing** database instead:

.. code-block:: bash

   cazy_webscraper download_cazy user@example.com -d cazy/cazyme.db

.. warning::
   ``-o`` and ``-d`` are mutually exclusive: one creates, the other appends. Use ``-f``/``--force``
   to overwrite an existing file at the ``-o`` path.

-------------------
Subcommand options
-------------------

``email``
  **Required.** Your email address, for NCBI Entrez.

``-o``, ``--db_output``
  Path at which to build a new database. Default: ``None``.

``-d``, ``--database``
  Path to an existing database to add data to. Default: ``None``.

``-s``, ``--subfamilies``
  Retrieve CAZy subfamily annotations as well as families. Default: ``False``.

``--cazy_data``
  Path to an already-downloaded CAZy text file, used instead of downloading afresh.
  Default: ``None``.

``--delete_old_relationships``
  Delete protein-to-family annotations that are in the local database but no longer in CAZy, for
  example where CAZy has moved a protein between families. Default: ``False``.

``--skip_ncbi_tax``
  Do not query NCBI Taxonomy for proteins that CAZy lists under several taxa. The first taxonomy
  CAZy lists is used instead. Default: ``False``.

``--ncbi_batch_size``
  Number of records per batch query sent to NCBI Entrez. Default: ``200``.

``-f``, ``--force``
  Overwrite an existing database at the ``--db_output`` path. Default: ``False``.

``-t``, ``--timeout``
  Seconds before a connection times out. Default: ``45``.

For the ``--classes``, ``--families``, ``--kingdoms``, ``--genera``, ``--species`` and
``--strains`` filters, the YAML configuration file and the shared housekeeping options, see
:doc:`Filtering and common options <filters>`.

.. note::
   ``download_cazy`` builds the dataset, so the filters that read an existing database do not apply
   to it: it has no ``--ec_filter``, ``--genbank_accessions`` or ``--uniprot_accessions``.

-----------------
Worked examples
-----------------

**Scrape only bacterial members of two classes:**

.. code-block:: bash

   cazy_webscraper download_cazy user@example.com -o cazy/cazyme.db \
       --classes GH,CE --kingdoms bacteria

**Scrape specific families, including their subfamilies:**

.. code-block:: bash

   cazy_webscraper download_cazy user@example.com -o cazy/cazyme.db \
       --families PL1,PL2,PL3 --subfamilies

**Restrict to named genera and species:**

.. code-block:: bash

   cazy_webscraper download_cazy user@example.com -o cazy/cazyme.db \
       --classes GH --genera Aspergillus,Trichoderma \
       --species Pythium ultimum

**Use a configuration file** instead of a long command line:

.. code-block:: bash

   cazy_webscraper download_cazy user@example.com -o cazy/cazyme.db \
       -c configuration_files/my_config.yaml

See :ref:`filter-combining` for more on combining filters, and :ref:`filter-yaml` for the
configuration file format.
