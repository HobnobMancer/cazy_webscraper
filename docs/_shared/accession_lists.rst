^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Providing a list of accessions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of selecting proteins with the filter flags, you can name them explicitly. Write one
accession per line in a plain text file and pass its path to ``--genbank_accessions`` (or, where
the subcommand supports it, ``--uniprot_accessions``).

.. warning::
   An accession list takes precedence over the filter flags. When ``--genbank_accessions`` is
   given, the filtering criteria are ignored entirely -- combining it with ``--classes``, for
   example, silently drops the ``--classes`` selection and processes only the listed accessions.
