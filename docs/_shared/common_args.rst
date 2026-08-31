Every ``get_*`` subcommand accepts the following, with the same meaning throughout.

**Filtering the protein set**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Argument
     - Description
   * - ``--classes``
     - CAZy classes to select. Comma-separated, no spaces.
   * - ``--families``
     - CAZy families or subfamilies to select. Comma-separated, no spaces.
   * - ``--kingdoms``
     - Taxonomic kingdoms to select: ``archaea``, ``bacteria``, ``eukaryota``, ``viruses``,
       ``unclassified``.
   * - ``--genera``
     - Genera to select. Comma-separated, no spaces.
   * - ``--species``
     - Species, written ``Genus species``. Selects all strains of the species.
   * - ``--strains``
     - Strains, written ``Genus species strain``.
   * - ``--ec_filter``
     - Restrict to proteins already annotated in the local database with at least one of these
       EC numbers.
   * - ``-c``, ``--config``
     - Path to a YAML configuration file holding the filters above. Default: ``None``.
   * - ``--genbank_accessions``
     - Path to a text file of GenBank accessions, one per line. Overrides the filters above.

**Utility arguments**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Argument
     - Description
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
