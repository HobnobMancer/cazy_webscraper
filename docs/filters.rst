.. _filtering:

==============================
Filtering and common options
==============================

Every ``cazy_webscraper`` subcommand that reads an existing local CAZyme database shares the same
system for choosing *which* proteins to act on, and the same set of housekeeping options. Those are
documented here once. Each subcommand's own page covers only the options specific to it.

.. note::
   The filters select proteins that are **already in your local CAZyme database**. They never cause
   new proteins to be added -- only ``download_cazy`` does that.

.. contents:: On this page
   :local:
   :depth: 1

.. _filter-overview:

-----------------------
How filtering works
-----------------------

With no filters, a subcommand acts on every protein in the database. Each filter you add selects a
group of proteins, and a protein is processed if it matches **at least one** criterion in **each**
category you have used.

There are four filter categories:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Category
     - Flags
   * - CAZy classification
     - ``--classes``, ``--families``
   * - Taxonomy
     - ``--kingdoms``, ``--genera``, ``--species``, ``--strains``
   * - Function
     - ``--ec_filter``
   * - Explicit list
     - ``--genbank_accessions``, ``--uniprot_accessions``

The order in which flags are written never matters.

.. warning::
   An explicit accession list overrides everything else. See :ref:`filter-accessions`.

.. _filter-cazy:

-------------------------------
CAZy class and family filters
-------------------------------

``--classes`` selects every family within the listed CAZy classes. ``--families`` selects specific
families or subfamilies. Separate multiple values with a single comma and no spaces.

.. code-block:: bash

   --classes GH,CE
   --families PL1,PL2,PL3

Both the class abbreviation and its full name are accepted, so these are equivalent:

.. code-block:: bash

   --classes GH,CE
   --classes Glycoside Hydrolases,Carbohydrate Esterases

.. warning::
   Families must use CAZy's own syntax. ``GH1`` is accepted; ``gh1`` and ``GlycosideHydrolases1``
   are not.

Used together, the two flags are **additive**: you get everything in the listed classes *plus*
everything in the listed families. To act on all of GH and CE, plus PL1, PL2 and PL3:

.. code-block:: bash

   --classes GH,CE --families PL1,PL2,PL3

.. tip::
   To select a CAZy subfamily, list it under ``--families``.

.. _filter-taxonomy:

-------------------
Taxonomic filters
-------------------

``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` restrict the selection by taxonomy. A
protein is selected if it matches at least one of the values given, so listing more taxa widens the
selection.

``--kingdoms`` accepts ``archaea``, ``bacteria``, ``eukaryota``, ``viruses`` and ``unclassified``.

.. warning::
   Kingdoms must be spelt the way CAZy spells them: ``eukaryot``\ **a**, not ``eukaryot``\ **e**.

.. note::
   Kingdom names are not case sensitive and may be listed in any order, so ``bacteria,eukaryota``
   and ``eukaryota,Bacteria`` are equivalent.

.. warning::
   Use standard scientific name formatting: capitalise the genus, lowercase the species.

   * ``Aspergillus niger`` -- correct
   * ``aspergillus niger`` -- incorrect
   * ``ASPERGILLUS NIGER`` -- incorrect

.. note::
   Naming a species selects **all** strains of that species. Use ``--strains`` to target individual
   strains, written ``Genus species strain``.

The taxonomic flags combine freely. To select all viral proteins, everything from *Aspergillus*,
two named *Layia* species and two specific *Trichoderma reesei* strains:

.. code-block:: bash

   --kingdoms viruses --genera Aspergillus \
       --species Layia carnosa,Layia chrysanthemoides \
       --strains Trichoderma reesei QM6a,Trichoderma reesei QM9414

.. _filter-ec:

----------------------
The EC number filter
----------------------

``--ec_filter`` restricts the selection to proteins already annotated with at least one of the
given EC numbers, letting you target a functional subset.

.. code-block:: bash

   --ec_filter "EC1.2.3.4,EC2.3.4.5"

.. note::
   Give complete EC numbers. Dashes (``-``) and asterisks (``*``) both stand in for missing digits.

   * ``EC1.2.3.-`` and ``EC1.2.3.*`` are accepted
   * ``EC1.2.3.`` and ``EC 1.2.3`` are not

   The ``EC`` prefix is optional: ``EC1.2.3.4`` and ``1.2.3.4`` are both fine.

.. warning::
   Quote the whole list when using dashes. Some shells read ``EC1.2.-.-`` as an attempt to invoke
   an option, so write ``--ec_filter "EC1.2.-.-"``.

.. note::
   Proteins matching **at least one** of the listed EC numbers are selected. Listing several EC
   numbers does not restrict the selection to proteins carrying all of them.

``--ec_filter`` reads the EC annotations *already stored in your local CAZyme database*, not those
held by any external database. If protein A is annotated with EC1.2.3.4 in UniProt but that
annotation has never been retrieved into your database, ``--ec_filter EC1.2.3.4`` will not select
protein A, because your database has no record of it.

.. tip::
   Run ``get_uniprot_data --ec`` first to populate EC annotations, then use ``--ec_filter`` on
   subsequent commands.

.. _filter-combining:

----------------------
Combining filters
----------------------

Filters from different categories can be combined freely to define a precise subset. These
examples use ``get_uniprot_data``, but the filters behave identically for every subcommand.

**Bacterial proteins from selected classes and families.** Everything in GH and CE, plus CE1, CE5
and CE8, restricted to bacteria:

.. code-block:: bash

   cazy_webscraper get_uniprot_data cazy/cazyme.db --pdb \
       --classes GH,CE --families CE1,CE5,CE8 --kingdoms bacteria

**Two genera within one class:**

.. code-block:: bash

   cazy_webscraper get_uniprot_data cazy/cazyme.db --ec --go \
       --classes GH --genera Aspergillus,Trichoderma

**A functionally defined subset.** Bacterial proteins from three classes that already carry at
least one of three EC annotations:

.. code-block:: bash

   cazy_webscraper get_uniprot_data cazy/cazyme.db --sequence \
       --classes GH,CE,CBM --kingdoms bacteria \
       --ec_filter "3.2.1.23,3.2.1.37,3.2.1.85"

.. tip::
   A trailing backslash (``\``) breaks a long command across several lines, which makes complex
   filter combinations much easier to read and to keep in a script.

.. _filter-accessions:

--------------------------------
Providing a list of accessions
--------------------------------

Instead of describing the protein set with filters, you can name its members. Write one accession
per line in a plain text file and pass its path to ``--genbank_accessions``, or to
``--uniprot_accessions`` where the subcommand supports it.

.. code-block:: bash

   cazy_webscraper get_uniprot_data cazy/cazyme.db --genbank_accessions my_proteins.txt

.. warning::
   An accession list takes precedence over the filter flags. When ``--genbank_accessions`` is
   given, the filtering criteria are ignored entirely: combining it with ``--classes`` silently
   drops the ``--classes`` selection and processes only the listed accessions.

Not every subcommand accepts both list flags. ``--uniprot_accessions`` is available on
``get_ncbi_taxs``, ``get_pdb_structures``, ``get_gtdb_taxs`` and ``extract_data``.

.. _filter-yaml:

---------------------------------
Configuration using a YAML file
---------------------------------

The same filters can be written into a YAML file and passed with ``-c``/``--config``. A config file
keeps long commands readable and makes a run reproducible and documentable.

A template lives in the repository at ``configuration_files/template-get_data_config.yaml``:

.. code-block:: yaml

    # Under 'classes', list the classes from which all proteins will be retrieved.
    # Under each class's own name, list the specific families/subfamilies wanted.
    # Write the FULL family name, e.g. 'GH1', not just its number, e.g. '1'.
    # List each family on a new line, indented once relative to the parent class,
    # with the name in quotation marks.
    classes:
    - "GH"
    - "CE"
    Glycoside Hydrolases (GHs):
    GlycosylTransferases (GTs):
    - "GT1"
    - "GT5"
    - "GT6"
    Polysaccharide Lyases (PLs):
    Carbohydrate Esterases (CEs):
    Auxiliary Activities (AAs):
    Carbohydrate-Binding Modules (CBMs):
    genera:
    - "Trichoderma"
    - "Aspergillus"
    species:
    - "Pythium ultimum"
    strains:
    kingdoms:  # Archaea, Bacteria, Eukaryota, Viruses, Unclassified
    ec:
    - "EC1.2.3.4"

.. attention::
   All of the tags above must be present in the file, even where the value is left empty:
   ``classes``, the six class names, ``genera``, ``species``, ``strains``, ``kingdoms`` and ``ec``.

.. warning::
   The ``ec`` tag is the one exception. Leaving it present but empty currently raises
   ``TypeError: 'NoneType' object is not iterable`` rather than being read as "no EC filter".
   Either list at least one EC number under it, or omit the tag and use ``--ec_filter`` instead.

Each value must be on its own line, indented, with the name in single or double quotation marks:

.. code-block:: yaml

    classes:
        - "GT"
        - "pl"
    Glycoside Hydrolases (GHs):
        - "GH1"
        - "GH2"

.. note::
   Filters given at the command line and in a config file are **combined**, not overridden. A
   protein is selected if it matches the criteria from either source. Passing ``--classes PL``
   alongside a config file listing ``CE`` selects both.

The ``ec`` tag is equivalent to ``--ec_filter``.

.. _filter-synonyms:

---------------------------
Synonyms for CAZy classes
---------------------------

Several synonyms are accepted for each CAZy class: both ``GH`` and ``Glycoside-Hydrolases`` resolve
to ``Glycoside Hydrolases (GHs)``, the name CAZy itself records. The accepted alternatives are
defined in ``cazy_dictionary.json`` in the ``cazy_webscraper`` repository. To supply your own set,
pass a JSON file with ``--cazy_synonyms``.

.. _common-options:

------------------
Common options
------------------

These behave identically across every subcommand that accepts them.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Description
   * - ``-c``, ``--config``
     - Path to a YAML configuration file holding the filters above. Default: ``None``.
   * - ``--cache_dir``
     - Path to use as the cache directory instead of the default. Default: ``None``.
   * - ``--cazy_synonyms``
     - Path to a JSON file of accepted CAZy class synonyms, if the defaults are not sufficient.
   * - ``-f``, ``--force``
     - Write into an existing cache directory rather than refusing to. Default: ``False``.
   * - ``--nodelete_cache``
     - Keep the existing contents of the cache directory. Default: ``False``, contents are deleted.
   * - ``-r``, ``--retries``
     - Number of times to retry a failed connection. Default: ``10``.

Some subcommands additionally accept:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Description
   * - ``-t``, ``--timeout``
     - Seconds before a connection times out. Default: ``45``. Available on ``download_cazy``,
       ``get_pdb_structures``, ``get_gtdb_taxs`` and ``extract_data``.
   * - ``-n``, ``--nodelete``
     - Keep the existing contents of the output directory. Default: ``False``. Available on
       ``get_pdb_structures``, ``get_gtdb_taxs`` and ``extract_data``.
   * - ``--batch_size``
     - Number of records per batch request. Default: ``150``, except ``extract_data`` where it is
       ``1000``. ``download_cazy`` calls its equivalent ``--ncbi_batch_size`` (default ``200``).

.. _global-options:

------------------
Global options
------------------

These five belong to ``cazy_webscraper`` itself rather than to any subcommand, so they must be
given **before** the subcommand name.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Description
   * - ``--version``
     - Print the version number and exit.
   * - ``--citation``
     - Print citation information and exit.
   * - ``-l``, ``--log``
     - Path to write a log file to. Default: ``None``, no log file is written.
   * - ``--sql_echo``
     - Set the SQLite engine ``echo`` parameter to ``True``, so SQLite prints its log messages.
   * - ``-v``, ``--verbose``
     - Enable verbose logging. This does *not* set the SQLite ``echo`` parameter.

.. code-block:: bash

   cazy_webscraper -v -l run.log get_uniprot_data cazy/cazyme.db

Placing them after the subcommand is an error. See :doc:`migration` if you are converting a
version 2 command, where these were accepted anywhere on the line.
