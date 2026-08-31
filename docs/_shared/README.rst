.. This directory holds reusable documentation fragments.

   They are pulled into pages with ``.. include:: _shared/<name>.rst`` and are
   excluded from the Sphinx build (see ``exclude_patterns`` in ``conf.py``), so
   they never render as pages of their own.

   Editing a fragment updates every page that includes it. That is the point:
   before adding a flag description or a caveat to a page, check whether it
   belongs in a fragment instead.

   Heading convention, so fragments slot into a page correctly:

     ===  page title      (pages only, never in a fragment)
     ---  section         (pages only, never in a fragment)
     ^^^  subsection      (fragments start here)
     """  sub-subsection

   Fragments describe *semantics* -- what a flag means, its syntax rules and its
   gotchas -- which are identical for every subcommand. Worked examples stay on
   the page that includes the fragment, because they differ per subcommand.
