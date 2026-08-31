^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CAZy class and family filters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``--classes`` restricts the operation to every family within the listed CAZy classes, and
``--families`` restricts it to specific families or subfamilies. Separate multiple values with a
single comma and no spaces.

Both the CAZy class abbreviation and its full name are accepted, so ``--classes GH,CE`` and
``--classes Glycoside Hydrolases,Carbohydrate Esterases`` are equivalent. The accepted spellings
are defined in ``cazy_dictionary.json`` in the ``cazy_webscraper`` repository, and can be
overridden with ``--cazy_synonyms``.

.. warning::
   Families must be written in CAZy's own syntax. ``GH1`` is accepted; ``gh1`` and
   ``GlycosideHydrolases1`` are not.

``--classes`` and ``--families`` combine additively: given both, the protein set is everything in
the listed classes **plus** everything in the listed families. The order of the flags does not
matter.

.. tip::
   To select a CAZy subfamily, list it under ``--families``.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Taxonomic filters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` restrict the operation by taxonomy. A
protein is included if it matches **at least one** of the criteria given, so these flags widen the
selection rather than narrowing it against each other.

``--kingdoms`` accepts ``archaea``, ``bacteria``, ``eukaryota``, ``viruses`` and ``unclassified``.

.. warning::
   Kingdoms must be spelt the way CAZy spells them: ``eukaryot``\ **a**, not ``eukaryot``\ **e**.

.. note::
   Kingdom names are not case sensitive, and may be listed in any order. ``bacteria,eukaryota`` and
   ``eukaryota,Bacteria`` are equally valid.

.. warning::
   Use standard scientific name formatting for ``--genera``, ``--species`` and ``--strains``:
   capitalise the genus, lowercase the species.

   * ``Aspergillus niger`` -- correct
   * ``aspergillus niger`` -- incorrect
   * ``ASPERGILLUS NIGER`` -- incorrect

.. note::
   Naming a species selects **all** strains of that species. Use ``--strains`` to target
   individual strains.
