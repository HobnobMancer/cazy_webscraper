============================
``cazy_webscraper`` Database
============================

To facilitate thorough interrogation of data retrieved from CAZy, minimise storing duplicate and redundant data, data retrieved from CAZy is stored in a local SQL database. Every CAZyme scraped from CAZy has the following data:

- Protein name
- CAZy (sub)family
- GenBank accession(s)

Each CAZyme may or may not have the following data, depending on the entry:

- EC number(s)
- UniProt acession(s)
- PDB accession(s)

------------------------------------------
Primary and non-primary GenBank accessions
------------------------------------------

Multiple GenBank accession numbers may be listed for a single CAZyme record at CAZy. CAZy records the "best model" for each CAZyme in bold (see `here for details <http://www.cazy.org/Help.html>`_). The bold accession is interpreted by ``cazy_webscraper`` as the **primary GenBank accession**. Where multiple GenBank accessions are written in bold, only the _first_ listed bold GenBank accession is listed as the **primary GenBank accession**. All other listed GenBank accessions are listed as **non-primary GenBank accessions**.

Each CAZyme in the local SQLite3 database is identified uniquely by its **primary GenBank accession**. If a CAZyme is associated with multiple families, it is represented as a single CAZyme record in the local SQLite database.

.. NOTE::
    GenBank protein sequences are not recorded in the local SQLite3 database. They are written to disk in a user-specified directory.

------------------------------------------
Primary and non-primary UniProt accessions
------------------------------------------

Multiple UniProt accessions may be listed for a single CAZyme record at CAZy. ``cazy_webscraper`` records accessions listed in bold as "**primary UniProt accessions**", and other accessions as **non-primary UniProt accessions** in the local database. 

--------------
PDB accessions
--------------

All PDB accessions associated with a CAZyme record at CAZy are recorded in the local SQLite3 database, and no differentiation is made between primary and non-primary accessions.

.. NOTE::
    Not all PDB accessions represented in a CAZyme record at CAZy are necessarily present in PDB. For example, some accessions are placeholders while structures are under embargo.

.. NOTE::
    PDB/RCSB protein structures are not recorded in the local SQLite3 database. They are written to disk in a user-specified directory.
