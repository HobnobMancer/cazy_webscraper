# cazy_webscraper

-------------------------------

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6343936.svg)](https://doi.org/10.5281/zenodo.6343936)
[![licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/cazy_webscraper/blob/master/LICENSE)  
[![CircleCI](https://circleci.com/gh/HobnobMancer/cazy_webscraper.svg?style=shield)](https://circleci.com/gh/HobnobMancer/cazy_webscraper)
[![codecov](https://codecov.io/gh/HobnobMancer/cazy_webscraper/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/cazy_webscraper)
[![Documentation Status](https://readthedocs.org/projects/cazy-webscraper/badge/?version=latest)](https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest)  
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cazy_webscraper/badges/version.svg)](https://anaconda.org/bioconda/cazy_webscraper)
[![Anaconda-Update Badge](https://anaconda.org/bioconda/cazy_webscraper/badges/latest_release_date.svg)](https://anaconda.org/bioconda/cazy_webscraper)  
[![pyani PyPi version](https://img.shields.io/pypi/v/cazy_webscraper "PyPI version")](https://pypi.python.org/pypi/cazy_webscraper)  
[![Downloads](https://static.pepy.tech/badge/cazy-webscraper)](https://pepy.tech/project/cazy-webscraper)  
[![CITATION.cff](https://github.com/HobnobMancer/cazy_webscraper/actions/workflows/main.yml/badge.svg)](https://github.com/HobnobMancer/cazy_webscraper/actions/workflows/main.yml)

--------------------------------

> [!CAUTION]
> This is the development branch. It is not stable.

> [!IMPORTANT]
> This branch is a significantly faster and more lightweight `cazy_webscraper`, capable of
> compiling a local CAZyme database containing every CAZyme record in CAZy in about 2 minutes
> on a standard office laptop. Data is written to the database as it is retrieved, so an
> interrupted run can be resumed by re-issuing the same command.

`cazy_webscraper` is an application for the automated retrieval of protein and annotation data from the [CAZy](http://wwww.cazy.org/) database and builds a local SQLite3 database, enabling users to integrate the dataset into analytical pipelines in a manner unachievable through the CAZy website. Successive queries can be collated into a single local database, and a log of each query is recorded in the database for transparency, reproducibility and shareability.

The full documentation can be found at [Read the Docs](https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest).

For full details see our publication in [Microbial Genomics](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001086), which includes example analyses.

> Hobbs, E. E. M, Gloster, T. M., Pritchard, L. (2023) cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets, _Microbial Genomics_, 9(8). [https://doi.org/10.1099/mgen.0.001086](https://doi.org/10.1099/mgen.0.001086)

--------------------------------

**Retrieve for user-defined datasets of interest.** `cazy_webscraper` can recover specified CAZy Classes and/or CAZy families and these queries can be further filtered by taxonomy at the Kingdom, genus, species or strain level. 

**Enhance the CAZy dataset.** `cazy_webscraper` can retrieve and add to a local CAZyme database additional data from the following external databases :

* **[GenBank](https://www.ncbi.nlm.nih.gov/genbank/):**  
    - Protein sequences 
    - Complete taxonomic lineages (domain, superkingdom, phylum, class, order and family)
    - Latest genomic assembly version accessions and IDs - GenBank and RefSeq (when available)
- **[UniProt](https://www.uniprot.org/):** 
    - UniProt ID/accession
    - UniParc sequence ID
    - Protein name
    - EC numbers
    - PDB accessions
    - GO terms
    - Protein sequence (and date sequence was last updated)

- **[RCSB Protein Data Bank (RSCB PDB)](https://www.rcsb.org/):**
    - Protein structure files - *Structure files are written to disk, **not** stored in the local CAZyme database*
- **[Genome Taxonomy Database (GTDB)](https://gtdb.ecogenomic.org/):**
    - Complete archaeal and bacterial lineages
- **[Pfam](https://www.ebi.ac.uk/interpro/):**
    - Pfam family annotations, which includes non-CAZyme domain annotations *(work in progress)*

**Get the data back out.** The `extract_data` subcommand writes any subset of a local CAZyme
database to CSV, JSON, FASTA files or a BLAST database.

## New in version 3.0.0

* A single `cazy_webscraper` command with subcommands, replacing the separate `cw_*` scripts.
* Significantly faster, and far less memory hungry - data is written to the database as it is
  retrieved rather than accumulated in memory or in cache files.
* Interrupted runs can be resumed by re-issuing the same command.
* Integration of [GO (GeneOntology) terms](https://www.geneontology.org/).
* Integration of Pfam family and domain annotations (work in progress).

## Migrating from version 2

The separate `cw_*` scripts are now subcommands of one `cazy_webscraper` command:

| Version 2 | Version 3 |
| --- | --- |
| `cazy_webscraper <email>` | `cazy_webscraper download_cazy <email>` |
| `cw_get_genbank_seqs <db> <email>` | `cazy_webscraper get_ncbi_seqs <db> <email>` |
| `cw_get_ncbi_taxs <db> <email>` | `cazy_webscraper get_ncbi_taxs <db> <email>` |
| `cw_get_genomics <db> <email>` | `cazy_webscraper get_ncbi_genomes <db> <email>` |
| `cw_get_uniprot_data <db>` | `cazy_webscraper get_uniprot_data <db>` |
| `cw_get_pdb_structures <db> <formats>` | `cazy_webscraper get_pdb_structures <db> --file_formats <formats>` |
| `cw_get_gtdb_taxs <db> <taxs>` | `cazy_webscraper get_gtdb_taxs <db> --taxs <taxs>` |
| `cw_query_database <db> <file_types>` | `cazy_webscraper extract_data <db> --file_types csv json` |
| `cw_extract_db_seqs <db> <source>` | `cazy_webscraper extract_data <db> --file_types fasta --source <source>` |
| `cw_get_db_schema` | Not yet migrated |

The two version 2 commands for getting data back *out* of a database were merged into
`extract_data`, selected with `--file_types csv|tsv|json|jsonl|fasta|fasta_dir|blastdb`.

`--version`, `--citation`, `-l`/`--log`, `--sql_echo` and `-v`/`--verbose` belong to
`cazy_webscraper` itself, so they go **before** the subcommand:
`cazy_webscraper -v get_ncbi_seqs ...`.

See the [migration guide](https://cazy-webscraper.readthedocs.io/en/latest/migration.html) for the
full mapping and the renamed arguments for each subcommand.

## Installation

### Conda

```bash
conda install -c bioconda cazy_webscraper
```

Please see the [`conda` documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) and [`bioconda` documentation](https://bioconda.github.io/) for further details.

### Pip

```bash
pip install cazy_webscraper
```

Please see the [`pip` documentation](https://pypi.org/project/pip/) for further details.

### From source

```bash
git clone --branch dev --sing-branch https://github.com/HobnobMancer/cazy_webscraper.git
cd cazy_webscraper
pip install .
```

## Quickstart

Build a local CAZyme database containing all of CAZy, then add NCBI protein sequences to it:

```bash
cazy_webscraper download_cazy user@email.com -o cazy.db
cazy_webscraper get_ncbi_seqs cazy.db user@email.com
```

Your email address is required by NCBI Entrez. It is not stored by `cazy_webscraper`.

## Subcommands

`download_cazy` builds a local CAZyme database; every `get_*` subcommand adds data to a database
that already exists; `extract_data` writes data back out to files.

| Subcommand | Description | Status |
| --- | --- | --- |
| `download_cazy` | Download CAZy and build a local SQLite database | Implemented |
| `get_ncbi_seqs` | Add protein sequences from NCBI | Implemented |
| `get_ncbi_taxs` | Add the latest NCBI taxonomic classifications | Implemented |
| `get_ncbi_genomes` | Add the genomic assembly each protein came from | Implemented |
| `get_uniprot_data` | Add UniProt accessions, names, EC numbers, PDB accessions and GO terms | Implemented |
| `get_pdb_structures` | Download PDB structure files; add each structure's method and resolution | Implemented |
| `get_pfams` | Add Pfam domain annotations from InterPro | Implemented |
| `get_gtdb_taxs` | Add GTDB taxonomic classifications for the genomes in the database | Implemented |
| `extract_data` | Write data out to CSV, TSV, JSON, JSON Lines, FASTA or a BLAST database | Implemented |
| Database schema printing | Print the schema of a local database | Not yet migrated to v3 |

**Every flag a subcommand takes is listed by its own help**, which is generated from the code and
so is always up to date:

```bash
cazy_webscraper <subcommand> --help
```

## Selecting a subset of proteins

The subcommands share a common set of filters, which can be combined to define a dataset of
interest. Without them a subcommand acts on the whole database.

| Filter | Selects | Available on |
| --- | --- | --- |
| `--classes` | CAZy classes, comma separated | all |
| `--families` | CAZy families and subfamilies, e.g. `GH1,GH2` | all |
| `--kingdoms` | `archaea`, `bacteria`, `eukaryota`, `viruses`, `unclassified` | all |
| `--genera`, `--species`, `--strains` | Taxonomy, e.g. `--species "Pythium ultimum"` | all |
| `-c`, `--config` | A YAML configuration file (see below) | all |
| `--ec_filter` | Proteins annotated with the given EC numbers | all except `download_cazy` |
| `--genbank_accessions` | A text file of GenBank accessions, one per line | all except `download_cazy` |
| `--uniprot_accessions` | A text file of UniProt accessions, one per line | `get_ncbi_taxs`, `get_pdb_structures`, `get_pfams`, `get_gtdb_taxs`, `extract_data` |

`download_cazy` has no EC or accession filters because it builds the database from CAZy, before
there are any EC annotations or accessions in it to filter against.

## Building a database

```bash
# the entire CAZy dataset
cazy_webscraper download_cazy user@email.com -o cazy.db

# only selected families, and only bacteria
cazy_webscraper download_cazy user@email.com -o cazy.db --families GH1,GH2,CE1 --kingdoms Bacteria

# driven by a configuration file
cazy_webscraper download_cazy user@email.com -o cazy.db --config my_config.yaml
```

Command line filters and a configuration file can be used together; the two are combined, so
`--families GH1 --config conf.yaml` retrieves GH1 *and* everything the file asks for.

## Adding data to a database

Each subcommand below writes to the database as it retrieves, so an interrupted run can be resumed
by re-issuing the same command. By default only proteins missing the data are fetched; add
`--update` to refresh what is already stored.

### Protein sequences, taxonomy and genomes from NCBI

```bash
cazy_webscraper get_ncbi_seqs cazy.db user@email.com
cazy_webscraper get_ncbi_taxs cazy.db user@email.com
cazy_webscraper get_ncbi_genomes cazy.db user@email.com
```

`get_ncbi_genomes` links each protein to the specific assembly it came from, and records both the
GenBank (`GCA_`) and RefSeq (`GCF_`) accession of that assembly.

### UniProt

```bash
# accessions, names and sequences
cazy_webscraper get_uniprot_data cazy.db

# also EC numbers, PDB accessions and GO terms
cazy_webscraper get_uniprot_data cazy.db --ec --pdb --go

# restrict to reviewed (Swiss-Prot) entries
cazy_webscraper get_uniprot_data cazy.db --swissprot
```

### PDB structures

PDB accessions come from UniProt, so run `get_uniprot_data --pdb` first.

```bash
# download structure files (default format: mmCif)
cazy_webscraper get_pdb_structures cazy.db --file_formats mmCif pdb -o structures/

# add each structure's experimental method and resolution, downloading no files
cazy_webscraper get_pdb_structures cazy.db --skip_download
```

### Pfam domains

Pfam matches are looked up by UniProt accession, so run `get_uniprot_data` first.

```bash
cazy_webscraper get_pfams cazy.db
```

### GTDB taxonomy

GTDB classifies genomes rather than proteins, so run `get_ncbi_genomes` first.

```bash
cazy_webscraper get_gtdb_taxs cazy.db --taxs bacteria
```

## Getting data back out

`extract_data` writes any subset of the database to files. The outputs are chosen with
`--file_types`:

| `--file_types` | Output |
| --- | --- |
| `csv` | One row per protein; columns chosen with `--include` |
| `tsv` | Same as `csv`, tab separated |
| `json` | One object per protein, keyed by protein accession |
| `jsonl` | JSON Lines - one self-contained JSON object per line, one per protein |
| `fasta` | A single FASTA file of all extracted sequences |
| `fasta_dir` | One FASTA file per protein |
| `blastdb` | A BLAST protein database (needs `makeblastdb` from BLAST+ on your `PATH`) |

```bash
# a table of proteins with their families, taxonomy and EC numbers
cazy_webscraper extract_data cazy.db --file_types csv --include class family kingdom organism ec

# every sequence as FASTA, plus a BLAST database
cazy_webscraper extract_data cazy.db --file_types fasta blastdb

# UniProt sequences for one family only
cazy_webscraper extract_data cazy.db --file_types fasta --source uniprot --families GH13
```

`--include` accepts `class`, `family`, `subfamily`, `kingdom`, `genus`, `organism`, `ec`, `pdb`,
`pfam`, `uniprot_acc`, `uniprot_name`, `genbank_seq` and `uniprot_seq`. Sequences are written with
the accession as the ID and the source database (`GenBank` or `UniProt`) as the description.

Output files are named `<prefix>_<database name>.<ext>` and written to `--output_dir` (default: the
working directory). Existing files are never silently replaced - the run stops and lists them
unless `--overwrite` is given.

The database is an ordinary SQLite3 file, so it can equally be queried directly with SQL, or
through a browser such as [DB Browser for SQLite](https://sqlitebrowser.org/).

## Configuring `cazy_webscraper` using a YAML file

The retrieval of data from CAZy, UniProt, GenBank and PDB can be configured at the command-line **and** via a YAML file.

The YAML file must have the following structure, specifically the YAML file must have the exact keys presented below and the values can be customised to configure the behaviour of `cazy_webscraper`:
```yaml
classes:  # classes from which all proteins will be retrieved
  - "GH"
  - "CE"
Glycoside Hydrolases (GHs):
GlycosylTransferases (GTs):
Polysaccharide Lyases (PLs):
  - "GT1"
  - "GT5"
  - "GT6"
Carbohydrate Esterases (CEs):
Auxiliary Activities (AAs):
Carbohydrate-Binding Modules (CBMs):
genera:  # list genera to be scraped
 - "Trichoderma"
 - "Aspergillus"
species:  # list species, this will scrape all strains under the species
- "Pythium ultimum"
strains:  # list specific strains to be scraped
kingdoms:  # Archaea, Bacteria, Eukaryota, Viruses, Unclassified
```

For configuring the retrieval of data from UniProt, GenBank and PDB (_but not CAZy) the additional `ec` tag can be included to limit the retrieval of data to proteins annotated with specific EC numbers.

When listing EC numbers, the 'EC' prefix can be included or excluded. For example, 'EC1.2.3.4' and '1.2.3.4' are accepted. Additionally, both dashes ('-') and astrixes ('*') can be used to represent missing digits, both '1.2.3.-' and '1.2.3.\*' are accepted.

`cazy_webscraper` performs a direct EC number comparison. Therefore, supplying `cazy_webscraper` with the EC number EC1.2.3.- will only retrieve protein specifically annotated with EC1.2.3.-. `cazy_webscraper` will **not** retrieve proteins will all completed EC numbers under EC1.2.3.-, thus, `cazy_webscraper` will **not** retrieve data for proteins annotated with EC1.2.3.1, EC1.2.3.2, EC1.2.3.3, etc.

Example configuration files, and an empty configuraiton file template are located in the `configuration_files/` directory of this repo.

## Database Schema

The central table is `Proteins`, keyed on `protein_id`/`protein_accession`. Shared values — CAZy
families, EC numbers, source organisms, genomes — are stored once and linked through association
tables (`Proteins_CazyFamilies`, `Proteins_Ecs`, `Proteins_Pdbs`, `Proteins_GoTerms`,
`Proteins_Pfams`, `Proteins_Genomes`).

For the entity-relationship diagram, a table-by-table reference, and the full `CREATE TABLE`
schema, see [the database page in the
documentation](https://cazy-webscraper.readthedocs.io/en/latest/database.html).

To inspect the schema of a database you already have:

```bash
sqlite3 <path to local CAZyme database> .schema
```

> [!WARNING]
> The schema changed in version 3: version 2's `Genbanks` table became `Proteins`, and every
> association table was renamed with it. A version 2 database cannot be used with version 3 — see
> [Migrating from version 2](#migrating-from-version-2). The diagram in `assets/` is the version 2
> schema and no longer reflects the database.


## Contributions

We welcome contributions and suggestions. You can raise issues at this repository, or fork the repository and submit pull requests, at the links below:

- [Issues](https://github.com/HobnobMancer/cazy_webscraper/issues)
- [Pull Requests](https://github.com/HobnobMancer/cazy_webscraper/pulls)

## License and copyright

MIT License

Copyright (c) 2025 University of St Andrews  
Copyright (c) 2025 University of Strathclyde  
Copyright (c) 2025 James Hutton Institute  

## Citation

If you use `cazy_webscraper`, please cite the following publication:

> Hobbs, E. E. M, Gloster, T. M., Pritchard, L. (2023) cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets, _Microbial Genomics_, 9(8). [https://doi.org/10.1099/mgen.0.001086](https://doi.org/10.1099/mgen.0.001086)

The supplementary information for this manuscript is available via the BioRxiv server, and in the `manuscript` directory in this repository.

cazy_webscraper depends on a number of tools. To recognise the contributions that the 
authors and developers have made, please also cite the following:

When making an SQLite database:
> Hipp, R. D. (2020) SQLite, available: https://www.sqlite.org/index.html.

Retrieving taxonomic, genomic or sequence data from NCBI:
> Cock, P.J.A., Antao, T., Chang, J.T., Chapman, B.A., Cox, C.J., Dalke, A., et al (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics, 25(11), 1422-1423.
> Wheeler,D.L., Benson,D.A., Bryant,S., Canese,K., Church,D.M., Edgar,R., Federhen,S., Helmberg,W., Kenton,D., Khovayko,O. et al (2005) Database resources of the National Centre for Biotechnology Information: Update, Nucleic Acid Research, 33, D39-D45

Retrieving data from UniProt:
> Cokelaer, T., Pultz, D., Harder, L. M., Serra-Musach, J., Saez-Rodriguez, J. (2013) BioServices: a common Python package to access biological Web Services programmatically, Bioinformatics, 19(24), 3241-3242.

Downloading protein structure files from RSCB PDB:
> Berman, H.M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T.N., Weissig, H., et al (2022) The Protein Data Bank, Nucleic Acids Research, 28(1), 235-242.
> Hamelryck, T., Manderick, B. (2003), PDB parser and structure class implemented in Python. Bioinformatics, 19 (17), 2308–2310

Retrieving and using taxonomic data from GTDB:
> Parks, D.H., Chuvochina, M., Rinke, C., Mussig, A.J., Chaumeil, P., Hugenholtz, P. (2022) GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy, Nucleic Acids Research, 50(D1), D785-D794.
