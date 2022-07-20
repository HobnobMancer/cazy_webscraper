===================================================
Retrieving genomic assembly data from NCBI Assembly
===================================================

``cazy_webscraper`` can be used to retrieve the latest genomic assembly data from the NCBI Assembly database 
for a set of proteins of interest in a local CAZyme database.

Querying NCBI is handled by the `BioPython <https://biopython.orgQ>`_ module ``Bio.Entrez``. 

-----------
Quick Start
-----------

To download the genomic assembly data for all proteins in a local CAZyme datbase, use the following command structure:

.. code-block:: bash

   cw_get_genomics <path to local CAZyme db> <user email address>

.. NOTE::
   The ``cw`` prefix on command is an abbreviation of ``cazy_webscraper``.

.. NOTE::
    The user email address is a requirement of NCBI.Entrez.
   
-------------------------------------
Storing the taxonomic classifications
-------------------------------------

The NCBI genomic assembly data retrieved from the NCBI Assembly database are stored in the 
``Genomes`` table in the local CAZyme database. 

Specifically, ``cazy_webscraper`` retrieves the sequence version accession from the ``Genbanks`` table in 
the local CAZyme database that meet the user's specified criteria. The sequence accessions are then split 
into GenBank version accessions and RefSeq version accessions.

For GenBank version accessions the following data is retrieved from NCBI Assembly:
* Assembly name
* GenBank version accession
* GenBank assembly ID
* RefSeq version accession (if available)
* RefSeq assembly ID (if available)

For RefSeq version accessions, the following data is retrieved from NCBI Assembly:
* Assembly name
* Refseq version accession
* RefSeq assembly ID

Protein sequence accessions in the ``Genbanks`` table are linked to their source genomic assembly via 
the ``Genbanks_Genomes`` table.

To link protein sequence accessions with their source genomic assembly, ``cazy_webscraper`` queries NCBI Assembly 
to retrieve the genomic assembly data for all genomes linked to any of the protein sequence accessions retrieved from 
the local CAZyme database. ``cazy_webscraper`` then downloads the feature table for each genomic assembly, and parses 
the feature table to identify the CAZymes associated with the genomic assembly. The downloaded feature table is left in 
its compressed format and cached.

--------------------
Command line options
--------------------

``database`` - **REQUIRED** Path to a local CAZyme database to add UniProt data to.

``email`` - **REQUIRED** User email address. Required by NCBI.

``--batch_size`` - Size of an individual batch query of NCBI sequence version accessions to NCBI. Default is 150.

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

``--update`` - Update the existing genomic assembly data in records in the ``Genomes`` table. By default, only NCBI assembly data not already in the local CAZyme database is added and existing data is not updated. Note updating data in the local CAZyme database ``Genomes`` table will overwrite existing data.

``--timeout`` - Connection timeout limit (s) when downloading feature tables from NCBI Assembly. Default 45s)

``--verbose``, ``-v`` - Enable verbose logging. This does **not** set the SQLite engine ``echo`` parameter to True. Default: False.


-----------
Basic Usage
-----------

The command-line options listed above can be used in combination to customise the retrieval of the NCBI 
taxonomic classifications for proteins of interest. Some options (e.g. ``--families`` and ``--classes``) define 
the broad group of proteins for which taxonomic data will be retrieved. Others filters (e.g. ``--species``) are used to filter and fine-tune the protein dataset for which data is retrieved.

The ``--classes``, ``--families``, ``--kingdoms``, ``--genera``, ``--species``, and ``--strains`` filteres are applied 
in the exactly same for retrieving data from CAZy, UniProt, and PDB. Examples of using these flags 
can be found in the ``cazy_webscraper`` and ``cw_get_uniprot_data`` tutorial in this documentation.

.. NOTE::
    To retrieve taxonomic information for members of specific CAZy subfamilies, list the subfamilies after the ``--families`` 
    flag.


-------------------------------------------
Retrieval of NCBI genomic assembly data
-------------------------------------------

The command for using ``cazy_webscraper`` for retrieving GenBank and RefSeq genomic assembly data
from the NCBI Assembly database is ``cw_get_genomics``.
