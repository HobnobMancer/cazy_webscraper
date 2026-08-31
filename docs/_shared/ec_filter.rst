^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The EC number filter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``--ec_filter`` restricts the operation to proteins already annotated with at least one of the
given EC numbers, letting you target a functional subset of a local CAZyme database.

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

``--ec_filter`` reads the EC annotations *already stored in the local CAZyme database*, not the
annotations held by any external database. If protein A is annotated with EC1.2.3.4 in UniProt but
that annotation has never been retrieved into your database, ``--ec_filter EC1.2.3.4`` will not
select protein A, because the local database has no record of it.

.. tip::
   Run ``get_uniprot_data --ec`` first to populate EC annotations, then use ``--ec_filter`` on
   subsequent commands.
