====================================
The local CAZyme database structure
====================================

``cazy_webscraper`` stores everything it retrieves in a single SQLite3 file. The schema is
normalised so that shared values -- a CAZy family, an EC number, a source organism, a genome -- are
stored once and referenced, rather than repeated against every protein.

The central table is ``Proteins``. Every other table either describes a protein
(``Uniprots``, ``NcbiTaxs``), or is linked to proteins through an association table
(``Proteins_CazyFamilies``, ``Proteins_Ecs``, and so on).

.. warning::
   The schema changed substantially in version 3: the central table was renamed from ``Genbanks``
   to ``Proteins``, and every association table with it. A database built with version 2 cannot be
   used with version 3. See :doc:`migration`.

-------------------------------
Which subcommand fills what
-------------------------------

Only ``download_cazy`` creates a database. Each ``get_*`` subcommand populates a different part of
it, so the tables present and populated in your database depend on which subcommands you have run.

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Subcommand
     - Tables it writes to
   * - ``download_cazy``
     - ``Proteins``, ``CazyFamilies``, ``Proteins_CazyFamilies``, ``Taxs``, ``Kingdoms``,
       ``NcbiTaxs``, ``Logs``
   * - ``get_ncbi_seqs``
     - ``Proteins`` (``sequence``, ``seq_update_date``)
   * - ``get_ncbi_taxs``
     - ``NcbiTaxs``, and ``Proteins.ncbi_tax_id``
   * - ``get_ncbi_genomes``
     - ``Genomes``, ``Proteins_Genomes``
   * - ``get_uniprot_data``
     - ``Uniprots``, ``Ecs``, ``Proteins_Ecs``, ``Pdbs``, ``Proteins_Pdbs``, ``GoTerms``,
       ``Proteins_GoTerms``, and ``Proteins.uniprot_id``
   * - ``get_pfams``
     - ``Pfams``, ``Proteins_Pfams``
   * - ``get_pdb_structures``
     - ``Pdbs`` (``method``, ``resolution``). Structure files go to disk, not the database.
   * - ``get_gtdb_taxs``
     - ``GtdbTaxs``, and ``Genomes.gtdb_tax_id``

-----------------
Schema diagram
-----------------

Arrows point from a foreign key to the table it references. Association tables are shown in grey.

.. graphviz::

   digraph cazyme_db {
     graph [rankdir=LR, splines=ortho, nodesep=0.35, ranksep=0.9, fontname="Helvetica"];
     node  [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=10,
            fillcolor="#eef4fb", color="#4a7ab5"];
     edge  [color="#7a7a7a", arrowsize=0.7];

     Proteins [label="Proteins\nprotein_id (PK)\nprotein_accession\nsequence\nsource",
               fillcolor="#dcebfa", color="#1f4e79", penwidth=2];

     Taxs      [label="Taxs\ntaxonomy_id (PK)\ngenus, species"];
     Kingdoms  [label="Kingdoms\nkingdom_id (PK)\nkingdom"];
     NcbiTaxs  [label="NcbiTaxs\nncbi_tax_id (PK)\ndomain..strain"];
     Uniprots  [label="Uniprots\nuniprot_id (PK)\nuniprot_accession\nswissprot, md5"];
     Genomes   [label="Genomes\ngenome_id (PK)\nassembly_name\ngbk/refseq accession"];
     GtdbTaxs  [label="GtdbTaxs\ngtdb_tax_id (PK)\nkingdom..species\nrelease"];
     CazyFamilies [label="CazyFamilies\nfamily_id (PK)\nfamily, subfamily"];
     Ecs       [label="Ecs\nec_id (PK)\nec_number"];
     Pdbs      [label="Pdbs\npdb_id (PK)\npdb_accession\nmethod, resolution"];
     Pfams     [label="Pfams\npfam_id (PK)\naccession, name\nrelease"];
     GoTerms   [label="GoTerms\ngo_id (PK)\ngoterm_id\ndescription"];
     Logs      [label="Logs\nlog_id (PK)\nrun provenance", fillcolor="#f6f2e7", color="#a08b52"];

     node [fillcolor="#f0f0f0", color="#909090", fontsize=9];
     P_fam  [label="Proteins_CazyFamilies"];
     P_ec   [label="Proteins_Ecs"];
     P_pdb  [label="Proteins_Pdbs"];
     P_go   [label="Proteins_GoTerms"];
     P_gen  [label="Proteins_Genomes"];
     P_pfam [label="Proteins_Pfams\nmatch_start, match_end\ninterpro_accession"];

     Proteins -> Taxs;
     Proteins -> NcbiTaxs;
     Proteins -> Uniprots;
     Taxs -> Kingdoms;
     Genomes -> GtdbTaxs;

     P_fam  -> Proteins; P_fam  -> CazyFamilies;
     P_ec   -> Proteins; P_ec   -> Ecs;
     P_pdb  -> Proteins; P_pdb  -> Pdbs;
     P_go   -> Proteins; P_go   -> GoTerms;
     P_gen  -> Proteins; P_gen  -> Genomes;
     P_pfam -> Proteins; P_pfam -> Pfams;

     {rank=same; Logs;}
   }

.. note::
   ``Proteins_Pfams`` is not a plain association table. A protein can match the same Pfam family at
   several positions along its sequence, so it carries its own primary key and stores the match
   coordinates (``match_start``, ``match_end``) plus the InterPro entry the family is integrated
   into (``interpro_accession``).

-------------------
Table reference
-------------------

^^^^^^^^^^
Proteins
^^^^^^^^^^

The central table. One row per protein, keyed on ``protein_id``, with the accession CAZy lists in
``protein_accession``. ``source`` records the database CAZy attributes that accession to (for
example ``ncbi``).

``sequence`` and ``seq_update_date`` are filled by ``get_ncbi_seqs``. The update date is what lets
``--update`` decide whether a newly retrieved sequence is actually newer than the stored one.

.. warning::
   Unlike version 2's ``Genbanks.genbank_accession``, ``Proteins.protein_accession`` has an index
   but **no UNIQUE constraint**. Uniqueness is enforced in application code when data is inserted,
   not by the database, so a direct write with ``sqlite3`` can introduce duplicate accessions
   without error.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CazyFamilies and Proteins_CazyFamilies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``CazyFamilies`` lists the CAZy families retrieved. When subfamilies are retrieved (``download_cazy
--subfamilies``) each subfamily is stored against its parent family.
``Proteins_CazyFamilies`` records which protein belongs to which family.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Taxs, Kingdoms and NcbiTaxs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Taxs`` holds the genus and species of each source organism as CAZy lists it, linked to
``Kingdoms``. ``NcbiTaxs`` holds the fuller lineage retrieved by ``get_ncbi_taxs``, from ``domain``
down to ``strain``, keyed on the NCBI Taxonomy ID.

A protein therefore has two taxonomic links: ``taxonomy_id`` (from CAZy) and ``ncbi_tax_id`` (from
NCBI). They can disagree, and the NCBI lineage is the more authoritative of the two.

^^^^^^^^^^^^^^^^^^^^^^^^^^
Genomes and GtdbTaxs
^^^^^^^^^^^^^^^^^^^^^^^^^^

``Genomes`` holds the genomic assemblies retrieved by ``get_ncbi_genomes``, with both the GenBank
and RefSeq version accessions and NCBI IDs. ``GtdbTaxs`` holds GTDB classifications, which attach
to the *genome* rather than the protein -- which is why ``get_ncbi_genomes`` must run before
``get_gtdb_taxs``. ``release`` records which GTDB release a classification came from.

^^^^^^^^^^
Uniprots
^^^^^^^^^^

Data retrieved by ``get_uniprot_data``: the accession, UniParc ID, protein name, whether the record
is in the reviewed SwissProt subset, the sequence with its MD5 and molecular weight, and the date
UniProt last updated the sequence.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Ecs, Pdbs, GoTerms and Pfams
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each of these holds one kind of annotation, linked to proteins through its association table:

* ``Ecs`` -- EC numbers, stored **without** the ``EC`` prefix.
* ``Pdbs`` -- PDB accessions from UniProt. ``get_pdb_structures`` adds the experimental
  ``method`` and ``resolution``.
* ``GoTerms`` -- GO terms and their descriptions, from ``get_uniprot_data --go``.
* ``Pfams`` -- Pfam families from InterPro, with the InterPro ``release`` the data came from.

.. note::
   Not every PDB accession CAZy or UniProt lists is necessarily present in PDB -- some are
   placeholders for structures still under embargo.

.. note::
   Structure files themselves are **not** stored in the database. They are written to disk in the
   directory given by ``--outdir``.

^^^^^^
Logs
^^^^^^

Every run that adds data appends a row, recording the date and time, the external database queried,
which annotation types were retrieved, every filter applied, and a reproduction of the command
line. This is what lets anyone using the database see how the dataset was compiled.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Contents
   * - ``log_id``
     - Autoincrementing ID
   * - ``date``, ``time``
     - When the run was initiated, in ISO format
   * - ``database``
     - External database queried, e.g. ``CAZy``, ``UniProt``, ``GenBank``
   * - ``retrieved_annotations``
     - Annotation types retrieved, e.g. ``EC number, PDB accession, Sequence``
   * - ``classes``, ``families``
     - CAZy class and family filters applied
   * - ``kingdoms``, ``genera_filter``, ``species_filter``, ``strains_filter``
     - Taxonomic filters applied
   * - ``ec_filter``
     - EC number filter applied
   * - ``cmd_line``
     - The command line as invoked

-----------------------------------------------
Inspecting the schema of your own database
-----------------------------------------------

.. warning::
   ``cw_get_db_schema`` was a version 2 command and has **not yet been migrated** to the version 3
   subcommand interface. Use ``sqlite3`` instead. See :doc:`migration`.

.. code-block:: bash

    sqlite3 <path to local CAZyme database> .schema

To list just the table names:

.. code-block:: bash

    sqlite3 <path to local CAZyme database> ".tables"

.. tip::
   Checking for a ``Genbanks`` table is a quick way to tell a version 2 database from a version 3
   one. Version 3 databases have ``Proteins`` instead.

------------------
The full schema
------------------

The schema of a database created by ``cazy_webscraper`` version 3:

.. code-block:: sql

    CREATE TABLE IF NOT EXISTS "Kingdoms" (
    	kingdom_id INTEGER NOT NULL,
    	kingdom VARCHAR,
    	PRIMARY KEY (kingdom_id),
    	UNIQUE (kingdom)
    );
    CREATE TABLE IF NOT EXISTS "GoTerms" (
    	go_id INTEGER NOT NULL,
    	goterm_id VARCHAR,
    	description VARCHAR,
    	PRIMARY KEY (go_id),
    	UNIQUE (go_id)
    );
    CREATE INDEX "ix_GoTerms_goterm_id" ON "GoTerms" (goterm_id);
    CREATE TABLE IF NOT EXISTS "GtdbTaxs" (
    	gtdb_tax_id INTEGER NOT NULL,
    	kingdom VARCHAR,
    	phylum VARCHAR,
    	tax_class VARCHAR,
    	tax_order VARCHAR,
    	family VARCHAR,
    	genus VARCHAR,
    	species VARCHAR,
    	release VARCHAR,
    	PRIMARY KEY (gtdb_tax_id),
    	UNIQUE (kingdom, phylum, tax_class, tax_order, family, genus, species, release)
    );
    CREATE TABLE IF NOT EXISTS "CazyFamilies" (
    	family_id INTEGER NOT NULL,
    	family VARCHAR NOT NULL,
    	subfamily VARCHAR,
    	PRIMARY KEY (family_id),
    	UNIQUE (family, subfamily)
    );
    CREATE INDEX fam_index ON "CazyFamilies" (family, subfamily);
    CREATE TABLE IF NOT EXISTS "NcbiTaxs" (
    	ncbi_tax_id INTEGER NOT NULL,
    	domain VARCHAR,
    	kingdom VARCHAR,
    	phylum VARCHAR,
    	tax_class VARCHAR,
    	tax_order VARCHAR,
    	family VARCHAR,
    	genus VARCHAR,
    	species VARCHAR,
    	strain VARCHAR,
    	PRIMARY KEY (ncbi_tax_id),
    	UNIQUE (ncbi_tax_id)
    );
    CREATE INDEX ncbi_index ON "NcbiTaxs" (ncbi_tax_id, genus, species);
    CREATE TABLE IF NOT EXISTS "Uniprots" (
    	uniprot_id INTEGER NOT NULL,
    	uniprot_accession VARCHAR,
    	uniparc_id VARCHAR,
    	name VARCHAR,
    	swissprot BOOLEAN,
    	sequence VARCHAR,
    	md5 VARCHAR,
    	molecular_weight INTEGER,
    	seq_update_date VARCHAR,
    	PRIMARY KEY (uniprot_id),
    	UNIQUE (uniprot_accession)
    );
    CREATE INDEX uniprot_option ON "Uniprots" (uniprot_id, uniprot_accession);
    CREATE TABLE IF NOT EXISTS "Ecs" (
    	ec_id INTEGER NOT NULL,
    	ec_number VARCHAR,
    	PRIMARY KEY (ec_id),
    	UNIQUE (ec_number)
    );
    CREATE INDEX "ix_Ecs_ec_number" ON "Ecs" (ec_number);
    CREATE TABLE IF NOT EXISTS "Pdbs" (
    	pdb_id INTEGER NOT NULL,
    	pdb_accession VARCHAR,
    	method VARCHAR,
    	resolution FLOAT,
    	PRIMARY KEY (pdb_id),
    	UNIQUE (pdb_accession)
    );
    CREATE INDEX pdb_idx ON "Pdbs" (pdb_accession);
    CREATE TABLE IF NOT EXISTS "Pfams" (
    	pfam_id INTEGER NOT NULL,
    	accession VARCHAR,
    	name VARCHAR,
    	description VARCHAR,
    	annotation_type VARCHAR,
    	release VARCHAR,
    	PRIMARY KEY (pfam_id),
    	UNIQUE (pfam_id)
    );
    CREATE INDEX pfam_idx ON "Pfams" (accession);
    CREATE TABLE IF NOT EXISTS "Logs" (
    	log_id INTEGER NOT NULL,
    	date VARCHAR,
    	time VARCHAR,
    	"database" VARCHAR,
    	retrieved_annotations VARCHAR,
    	classes VARCHAR,
    	families VARCHAR,
    	kingdoms VARCHAR,
    	genera_filter VARCHAR,
    	species_filter VARCHAR,
    	strains_filter VARCHAR,
    	ec_filter VARCHAR,
    	cmd_line VARCHAR,
    	PRIMARY KEY (log_id)
    );
    CREATE TABLE IF NOT EXISTS "Taxs" (
    	taxonomy_id INTEGER NOT NULL,
    	genus VARCHAR,
    	species VARCHAR,
    	kingdom_id INTEGER,
    	PRIMARY KEY (taxonomy_id),
    	UNIQUE (genus, species),
    	FOREIGN KEY(kingdom_id) REFERENCES "Kingdoms" (kingdom_id)
    );
    CREATE INDEX organism_option ON "Taxs" (taxonomy_id, genus, species);
    CREATE TABLE IF NOT EXISTS "Genomes" (
    	genome_id INTEGER NOT NULL,
    	assembly_name VARCHAR,
    	gbk_version_accession VARCHAR,
    	gbk_ncbi_id INTEGER,
    	refseq_version_accession VARCHAR,
    	refseq_ncbi_id INTEGER,
    	gtdb_tax_id INTEGER,
    	PRIMARY KEY (genome_id),
    	UNIQUE (assembly_name, gbk_version_accession, refseq_version_accession),
    	FOREIGN KEY(gtdb_tax_id) REFERENCES "GtdbTaxs" (gtdb_tax_id)
    );
    CREATE INDEX genome_options ON "Genomes" (assembly_name, gbk_version_accession, refseq_version_accession);
    CREATE TABLE IF NOT EXISTS "Proteins" (
    	protein_id INTEGER NOT NULL,
    	protein_accession VARCHAR,
    	sequence VARCHAR,
    	seq_update_date VARCHAR,
    	taxonomy_id INTEGER,
    	ncbi_tax_id INTEGER,
    	uniprot_id INTEGER,
    	source VARCHAR,
    	PRIMARY KEY (protein_id),
    	FOREIGN KEY(taxonomy_id) REFERENCES "Taxs" (taxonomy_id),
    	FOREIGN KEY(ncbi_tax_id) REFERENCES "NcbiTaxs" (ncbi_tax_id),
    	FOREIGN KEY(uniprot_id) REFERENCES "Uniprots" (uniprot_id)
    );
    CREATE INDEX "ix_Proteins_protein_accession" ON "Proteins" (protein_accession);
    CREATE TABLE IF NOT EXISTS "Proteins_Genomes" (
    	protein_id INTEGER NOT NULL,
    	genome_id INTEGER NOT NULL,
    	PRIMARY KEY (protein_id, genome_id),
    	FOREIGN KEY(protein_id) REFERENCES "Proteins" (protein_id),
    	FOREIGN KEY(genome_id) REFERENCES "Genomes" (genome_id)
    );
    CREATE TABLE IF NOT EXISTS "Proteins_CazyFamilies" (
    	protein_id INTEGER NOT NULL,
    	family_id INTEGER NOT NULL,
    	PRIMARY KEY (protein_id, family_id),
    	FOREIGN KEY(protein_id) REFERENCES "Proteins" (protein_id),
    	FOREIGN KEY(family_id) REFERENCES "CazyFamilies" (family_id)
    );
    CREATE TABLE IF NOT EXISTS "Proteins_Ecs" (
    	protein_id INTEGER NOT NULL,
    	ec_id INTEGER NOT NULL,
    	PRIMARY KEY (protein_id, ec_id),
    	FOREIGN KEY(protein_id) REFERENCES "Proteins" (protein_id),
    	FOREIGN KEY(ec_id) REFERENCES "Ecs" (ec_id)
    );
    CREATE TABLE IF NOT EXISTS "Proteins_Pdbs" (
    	protein_id INTEGER NOT NULL,
    	pdb_id INTEGER NOT NULL,
    	PRIMARY KEY (protein_id, pdb_id),
    	FOREIGN KEY(protein_id) REFERENCES "Proteins" (protein_id),
    	FOREIGN KEY(pdb_id) REFERENCES "Pdbs" (pdb_id)
    );
    CREATE TABLE IF NOT EXISTS "Proteins_GoTerms" (
    	protein_id INTEGER NOT NULL,
    	go_id INTEGER NOT NULL,
    	PRIMARY KEY (protein_id, go_id),
    	FOREIGN KEY(protein_id) REFERENCES "Proteins" (protein_id),
    	FOREIGN KEY(go_id) REFERENCES "GoTerms" (go_id)
    );
    CREATE TABLE IF NOT EXISTS "Proteins_Pfams" (
    	protein_pfam_id INTEGER NOT NULL,
    	protein_id INTEGER,
    	pfam_id INTEGER,
    	interpro_accession VARCHAR,
    	match_start INTEGER,
    	match_end INTEGER,
    	PRIMARY KEY (protein_pfam_id),
    	UNIQUE (protein_id, pfam_id, match_start, match_end),
    	FOREIGN KEY(protein_id) REFERENCES "Proteins" (protein_id),
    	FOREIGN KEY(pfam_id) REFERENCES "Pfams" (pfam_id)
    );
