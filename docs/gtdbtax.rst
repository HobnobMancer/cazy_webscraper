=========================================
Retrieving GTDB Taxonomic Classifications
=========================================

``cazy_webscraper`` can be used to retrieve the latest taxonomic classification from the `GTDB <https://gtdb.ecogenomic.org/>`_ taxonomy database.
for a set of proteins of interest in a local CAZyme database.

.. Note::
    As in the GTDB database, GTDB taxonomic classifications are retrieved and associated with genomes stored 
    in the local CAZyme database. To retrieve GTDB taxonomic classifications the genomic data for the 
    proteins of interest **must** be listed in the local CAZyme database.

-----------
Quick Start
-----------

To download the GTDB taxonomic classifications for all proteins in a local CAZyme datbase, use the following command structure:

.. code-block:: bash

   cw_get_gtdb_taxs <path to local CAZyme db> <kingdoms>

.. NOTE::
   The ``cw`` prefix on command is an abbreviation of ``cazy_webscraper``.

.. NOTE::
    GTDB only provides data for bacteria and archaeal genomes.

   
-------------------------------------
Storing the taxonomic classifications
-------------------------------------

The GTDB taxonomy classifications retrieved from the GTDB Taxonomy database are stored in the 
``GtdbTaxs`` table in the local CAZyme database. 

Each unique organism strain retrieved from GTDB is stored as a unique record in the ``GtdbTaxs`` table, which lists for each record the:
* Superkingdom (referred to as the kingdom)
* Phylum
* Class (called tax_class in the database due to keyword clash with Python)
* Order (called tax_order in the database due to keyword clash with SLQ)
* Family
* Genus
* Species
* Strain
* Release - the GTDB release from which the lineage was retrieved

The child prteins for each taxonomy record in the ``GtdbTaxs`` table is identified by the 
including a ``gtdb_tax_id`` from the ``GtdbTaxs`` table in the respecitve ``Genomes`` table records.

--------------------
Command line options
--------------------

``database`` - **REQUIRED** Path to a local CAZyme database to add UniProt data to.

``taxs`` - **REQUIRED** Kingdoms to get lineages from. Accepts 'archaea' and/or 'bacteria'. Separate with a single space. Order does not matter.
Determines which datafiles are retrieved from GTDB.

``--archaea_file`` - Path to GTDB archaea data file. Default: None, download latest dataset from GTDB.

``--bacteria_file`` - Path to GTDB bacteria data file. Default: None, download latest dataset from GTDB.

.. NOTE::
    The filenames of provided GTDB data files must match the filename format used by GTDB, to allow 
    ``cazy_webscraper`` to retrieve the release number of the dataset.

``--cache_dir`` - Path to cache dir to be used instead of default cache dir path.

``--cazy_synonyms`` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

``--config``, ``-c`` - Path to a configuration YAML file. Default: None.

``--classes`` - list of classes to retrieve UniProt data for.

``--ec_filter`` - List of EC numbers to limit the retrieval of structure files to proteins with at least one of the given EC numbers *in the local CAZyme database*.

``--families`` - List of CAZy (sub)families to retrieve UniProt protein data for.

``--genbank_accessions`` - Path to text file containing a list of GenBank accessions to retrieve protein data for. A unique accession per line.

``--genera`` - List of genera to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given genera.

``--kingdoms`` - List of taxonomy kingdoms to retrieve UniProt data for.

``--log``, ``-l`` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

``--nodelete_cache`` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

``--retries``, ``-r`` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

``--sql_echo`` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

``--species`` - List of species (organsim scientific names) to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given species.

``--strains`` - List of species strains to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given strains.

``--uniprot_accessions`` - Path to text file containing a list of UniProt accessions to retrieve protein data for. A unique accession per line.

``--update_genome_lineage_gbk`` - Update Genome GTDB lineage. Default: only add lineages to Genomes without a lineage.

``--verbose``, ``-v`` - Enable verbose logging. This does **not** set the SQLite engine ``echo`` parameter to True. Default: False.



-----------
Basic Usage
-----------

The command-line options listed above can be used in combination to customise the retrieval of the GTDB 
taxonomic classifications for proteins of interest. Some options (e.g. ``--families`` and ``--classes``) define 
the broad group of proteins for which taxonomic data will be retrieved. Others filters (e.g. ``--species``) are used to filter and fine-tune the protein dataset for which data is retrieved.

The ``--classes``, ``--families``, ``--kingdoms``, ``--genera``, ``--species``, and ``--strains`` filteres are applied 
in the exactly same for retrieving data from CAZy, UniProt, and PDB. Examples of using these flags 
can be found in the ``cazy_webscraper`` and ``cw_get_uniprot_data`` tutorial in this documentation.

.. NOTE::
    To retrieve taxonomic information for members of specific CAZy subfamilies, list the subfamilies after the ``--families`` 
    flag.


-------------------------------------------
Retrieval of GTDB taxonomic classifications
-------------------------------------------

The command for using ``cazy_webscraper`` for retrieving taxonomic classifications 
from the GTDB Taxonomy database is ``cw_get_gtdb_taxs``.
