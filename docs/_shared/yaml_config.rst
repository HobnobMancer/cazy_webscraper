The same filters available at the command line can be written into a YAML file and passed with
``-c``/``--config``. A config file keeps long commands readable and makes a run reproducible and
documentable.

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
   The ``ec`` tag is the one exception: leaving it present but empty currently raises
   ``TypeError: 'NoneType' object is not iterable`` instead of being treated as "no EC filter".
   Either list at least one EC number under it, or omit the tag and supply EC numbers with
   ``--ec_filter`` instead.

Each value must be on its own line, indented, with the name in single or double quotation marks:

.. code-block:: yaml

    classes:
        - "GT"
        - "pl"
    Glycoside Hydrolases (GHs):
        - "GH1"
        - "GH2"

.. note::
   Filters given at the command line and in a config file are combined, not overridden. A protein
   is selected if it matches the criteria from either source.
