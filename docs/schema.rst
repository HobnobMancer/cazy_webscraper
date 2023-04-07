===================================
Retrieving structure files from PDB
===================================

The schema of a local CAZyme database can be retrieved using ``cazy_webscraper``:

.. code-block:: bash
    cw_get_db_schema <path to local CAZyme database>


Alternatively, `sqlite3` can be used to retrieve the schema:

.. code-block:: bash
    sqlite3 <path to local CAZyme database> .schema


A visual representation of the db schema when using `cazy_webscraper` version >= 2.3.0 can be found `here <https://hobnobmancer.github.io/cazy_webscraper/database_schema.pdf>`_.

As of `cazy_webscraper` version >= 2.3.0, the schema of a local CAZyme database will be:

.. code-block:: bash

    CREATE TABLE IF NOT EXISTS "Kingdoms" (
            kingdom_id INTEGER NOT NULL, 
            kingdom VARCHAR, 
            PRIMARY KEY (kingdom_id), 
            UNIQUE (kingdom)
    );
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
            uniprot_name VARCHAR, 
            sequence VARCHAR, 
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
            PRIMARY KEY (pdb_id), 
            UNIQUE (pdb_accession)
    );
    CREATE INDEX pdb_idx ON "Pdbs" (pdb_accession);
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
    CREATE TABLE IF NOT EXISTS "Genbanks" (
            genbank_id INTEGER NOT NULL, 
            genbank_accession VARCHAR, 
            sequence VARCHAR, 
            seq_update_date VARCHAR, 
            taxonomy_id INTEGER, 
            ncbi_tax_id INTEGER, 
            uniprot_id INTEGER, 
            PRIMARY KEY (genbank_id), 
            UNIQUE (genbank_accession), 
            FOREIGN KEY(taxonomy_id) REFERENCES "Taxs" (taxonomy_id), 
            FOREIGN KEY(ncbi_tax_id) REFERENCES "NcbiTaxs" (ncbi_tax_id), 
            FOREIGN KEY(uniprot_id) REFERENCES "Uniprots" (uniprot_id)
    );
    CREATE INDEX "ix_Genbanks_genbank_accession" ON "Genbanks" (genbank_accession);
    CREATE TABLE IF NOT EXISTS "Genbanks_Genomes" (
            genbank_id INTEGER NOT NULL, 
            genome_id INTEGER NOT NULL, 
            PRIMARY KEY (genbank_id, genome_id), 
            FOREIGN KEY(genbank_id) REFERENCES "Genbanks" (genbank_id), 
            FOREIGN KEY(genome_id) REFERENCES "Genomes" (genome_id)
    );
    CREATE TABLE IF NOT EXISTS "Genbanks_CazyFamilies" (
            genbank_id INTEGER NOT NULL, 
            family_id INTEGER NOT NULL, 
            PRIMARY KEY (genbank_id, family_id), 
            FOREIGN KEY(genbank_id) REFERENCES "Genbanks" (genbank_id), 
            FOREIGN KEY(family_id) REFERENCES "CazyFamilies" (family_id)
    );
    CREATE TABLE IF NOT EXISTS "Genbanks_Ecs" (
            genbank_id INTEGER NOT NULL, 
            ec_id INTEGER NOT NULL, 
            PRIMARY KEY (genbank_id, ec_id), 
            FOREIGN KEY(genbank_id) REFERENCES "Genbanks" (genbank_id), 
            FOREIGN KEY(ec_id) REFERENCES "Ecs" (ec_id)
    );
    CREATE TABLE IF NOT EXISTS "Genbanks_Pdbs" (
            genbank_id INTEGER NOT NULL, 
            pdb_id INTEGER NOT NULL, 
            PRIMARY KEY (genbank_id, pdb_id), 
            FOREIGN KEY(genbank_id) REFERENCES "Genbanks" (genbank_id), 
            FOREIGN KEY(pdb_id) REFERENCES "Pdbs" (pdb_id)
    );
