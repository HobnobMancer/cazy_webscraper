==========================================
Retrieving Sequences from GenBank Tutorial
==========================================

Tutorial coming soon!

------------------------------
Providing a list of accessions
------------------------------

Instead of retrieving protein sequences for all CAZymes matching a defined set of criteria, 
``cw_get_genbank_seqs`` can retrieve protein sequences a set of CAZymes defined by their 
GenBank and/or UniProt accession.

The flag ``--genbank_accessions`` can be used to provide ``cw_get_genbank_seqs`` a list of GenBank accessions 
to identify the specific set of CAZymes to retrieve protein sequences for.

The flag ``--uniprot_accessions`` can be used to provide ``cw_get_genbank_seqs`` a list of UniProt accessions 
to identify the specific set of CAZymes to retrieve protein sequences for.

In both instances (for ``--genbank_accessions`` and ``--uniprot_accessions``) the list of respective accessions 
are provided via a plain text file, with a unique protein accession of each line. The path to this file is 
then passed to ``cw_get_genbank_seqs`` via the respective ``--genbank_accessions`` and ``--uniprot_accessions`` flag.

``--genbank_accessions`` and ``--uniprot_accessions`` can be used at the same time to define all 
CAZymes of interest.

.. WARNING::
   ``--genbank_accessions`` and ``--uniprot_accessions`` take president over the filter flags.

   When either ``--genbank_accessions`` or ``--uniprot_accessions`` is used, ``cw_get_genbank_seqs`` will 
   **not** retrieve any CAZymes from the local database matching a set of criteria.

   Therefore, if ``--genbank_accessions`` and ``--classes`` are used, ``cw_get_genbank_seqs`` will ignore 
   the ``--classes`` flag and only retrieve protein sequences for the proteins listed in the file provided via 
   the ``--genbank_accessions``.
