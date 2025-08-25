=========================================
Retrieving NCBI Taxonomic Classifications
=========================================

``cazy_webscraper`` can be used to retrieve the latest taxonomic classification from the NCBI Taxonomy database 
for a set of proteins of interest in a local CAZyme database.

Querying NCBI is handled by the `BioPython <https://biopython.orgQ>`_ module ``Bio.Entrez``. 

To download the NCBI taxonomic classifications for all proteins in a local CAZyme datbase, use the following command structure:

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs <path to local CAZyme db> <user email address>

.. NOTE::
    The user email address is a requirement of NCBI.Entrez.
   
-------------------------------------
Storing the taxonomic classifications
-------------------------------------

The NCBI taxonomy classifications retrieved from the NCBI Taxonomy database are stored in the 
``NcbiTaxs`` table in the local CAZyme database. 

Each unique organism strain retrieved from NCBI is stored as a unique record in the ``NcbiTaxs`` table, which lists for each record the:
* Superkingdom (referred to as the kingdom)
* Phylum
* Class (called tax_class in the database due to keyword clash with Python)
* Order (called tax_order in the database due to keyword clash with SLQ)
* Family
* Genus
* Species
* Strain

The child prteins for each taxonomy record in the ``NcbiTaxs`` table is identified by the 
including a ``ncbi_tax_id`` from the ``NcbiTaxs`` table in the respecitve ``Genbanks`` table records.

--------------------
Command line options
--------------------

REQUIRED arguments 
    ``database``              Path to local CAZy database
    ``email``                 User email address, requirement of NCBI-Entrez. The address is not stored by cazy_webscraper

OPTIONAL arguments

Filtering arguments:
    ``--classes`` CLASSES     Classes from which all families are to be scraped. Separate classes by ',' (default: None)
    ``--families`` FAMILIES   Families to scrape. Separate families by commas 'GH1,GH2' (default: None)
    ``--kingdoms`` KINGDOMS   Kingdoms to scrape. Separate by a single comma. Options= archaea, bacteria, eukaryota, viruses, unclassified (not case
                                                sensitive) (default: None)
    ``--genera`` GENERA       Genera to restrict the scrape to (default: None)
    ``--species`` SPECIES     Species (written as Genus Species) to restrict the scrape to (default: None)
    ``--strains`` STRAINS     Specific strains of organisms to restrict the scrape to (written as Genus Species Strain) (default: None)
    ``--ec_filter`` EC_FILTER
                                                Limit retrieval to proteins annotated with the provided EC numbers. Separate EC numbers with single commas (default: None)

Operational arguments:
    ``-c`` Config file, ``--config`` Config file
                                                Path to configuration file. Default: None, scrapes entire database (default: None)
    ``--genbank_accessions`` GENBANK_ACCESSIONS
                                                Path to a text file containing a list of GenBank accessions to retrieve data for. When used, accessions will not be retrieved
                                                from the local db using the filtering criteria (default: None)
    ``--uniprot_accessions`` UNIPROT_ACCESSIONS
                                                Path to a text file containing a list of UniProt accessions to retrieve data for. When used, accessions will not be retrieved
                                                from the local db using the filtering criteria (default: None)
    ``--lineage_cache`` LINEAGE_CACHE
                                                Path to previously cached lineage dict containing lineages and protein seq acc (called lineage_data.json) (default: None)
    ``-F``, ``--file_only``       Only add seqs provided via JSON and/or FASTA file. Do not retrieve data from NCBI (default: False)

Utility arguments:
    ``--batch_size`` BATCH_SIZE
                                                Batch size for queries sent to NCBI.Entrez (default: 150)
    ``--cache_dir`` CACHE_DIR
                                                Path to cache directory (default: None)
    ``--cazy_synonyms`` CAZY_SYNONYMS
                                                Path to JSON file containing CAZy class synoymn names (default: None)
    ``-f``, ``--force``           Force writing in existing cache directory (default: False)
    ``--nodelete_cache``      Do not delete content in existing cache dir (default: False)
    ``-r`` RETRIES, ``--retries`` RETRIES
                                                Number of times to retry scraping a family or class page if error encountered (default: 10)

------------------------------------------------------------
Retrieving taxonomies for specific CAZy classes and families
------------------------------------------------------------

The ``--classes`` and ``--families`` flags from scraping data from CAZy are applied in the extact same way 
for retrieving taxonomy data from NCBI.


For example, if you want to retrieve protein data for all CAZymes from Glycoside Hydrolase and Carbohydrate Esterases then use the command:

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs cazy/cazyme.db dummyEmail@domain.com --classes GH,CE

For 
example, to retrieve protein data for all proteins in PL1, PL2 and PL3 in the local CAZyme database, use the 
following command:

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs cazy/cazyme.db dummyEmail@domain.com --families PL1,PL2,PL3


As with scraping data from CAZy, the ``--classes`` and ``--families`` flags can be combined. To retrieve 
protein data for all CAZymes in PL1, PL2, PL3 and *all* of GH and CE both:

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs cazy/cazyme.db dummyEmail@domain.com --families PL1,PL2,PL3 --classes GH,CE

--------------------------
Applying taxonomic filters
--------------------------


The ``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` flags can be used to refine the dataset 
of proteins to retrieve protein data by taxonomy. These flags are applied in the exact same way as they 
are used for the scraping of data from CAZy. Only proteins in the local CAZyme database and 
matching at least on of the provided taxonomy criteria will have data retrieved from NCBI taxonomy.


For example, if you want to retrieve data for all CAZymes in a local CAZyme database from bacterial and eukaryotic species, then use the command 

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs cazy/cazyme.db dummyEmail@domain.com --kingdoms bacteria,eukaryota

-------------------------
Applying EC number filter
-------------------------

The retrieval of taxonomic data from NCBI can also be limited to proteins in a local CAZyme database that are
annotated with specific EC numbers.

Having previously retrieved EC number annotations from UniProt and adding them to the local CAZyme database, you may 
wish to retrieve protein data for CAZymes annotated with specific EC numbers. To do this add the 
``--ec_filter`` flag to the command, follwed by a list of EC numbers.

.. code-block:: bash
   
   cazy_webscraper get_ncbi_taxs cazy/cazyme.db dummyEmail@domain.com --ec_filter "EC1.2.3.4,EC2.3.4.5"

.. NOTE::
    Provide complete EC numbers.  
    EC1.2.3.- and EC1.2.3.* are accepted.  
    EC1.2.3. and EC 1.2.3 are **not** accepted.  
    If using dashes to represent missing digits in EC numbers, it is recommended to bookend the entire in quotation marks.  
    ``--ec_filter`` is based upon EC number annotations stored within the local CAZyme database. 

---------------------
Combining all filters
---------------------

The ``--classes``, ``--families``, ``--ec_filter``, ``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` flags can 
be used in any combination to define a specific subset of proteins in the local CAZyme database for whom
taxonomic data will be retrieved from NCBI.

Below we run through 3 example commands of combining these flags, and the resulting behaviour.

**Example 1:**
To taxonomic data for all CAZymes in GH, GT, CE1, CE5 and CE8, and which are derived from baceterial species, we use the command:

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs cazy/cazyme.db dummyEmail@domain.com --classes GH,CE --families CE1,CE5,CE8 --kingdoms bacteria


**Example 2:**
To taxonomic data for all CAZymes in GH and which are derived from *Aspegillus* and *Trichoderma* species, we use the command:

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs cazy/cazyme.db dummyEmail@domain.com --classes GH --genera Aspegillus,Trichoderma


**Example 3:**
To taxonomic classifications for all CAZymes in GH,CE and CBM which are derived from baceterial species and are annotated with at least one of 
EC3.2.1.23, EC3.2.1.37 and EC3.2.1.85, we use the command:

.. code-block:: bash

   cazy_webscraper get_ncbi_taxs cazy/cazyme.db dummyEmail@domain.com --classes GH,CE,CBM --kingdoms bacteria --ec_filter "3.2.1.23,3.2.1.37,3.2.1.85"


------------------------------
Providing a list of accessions
------------------------------

Instead of retrieving taxonomic data for all CAZymes matching a defined set of criteria, 
``cazy_webscraper get_ncbi_taxs`` can retrieve taxonomic data for a set of CAZymes defined by their 
GenBank and/or UniProt accession.

The flag ``--genbank_accessions`` can be used to provide ``cazy_webscraper get_ncbi_taxs`` a list of GenBank accessions 
to identify the specific set of CAZymes to retrieve taxonomic data for.

The flag ``--uniprot_accessions`` can be used to provide ``cazy_webscraper get_ncbi_taxs`` a list of UniProt accessions 
to identify the specific set of CAZymes to retrieve taxonomic data for.

In both instances (for ``--genbank_accessions`` and ``--uniprot_accessions``) the list of respective accessions 
are provided via a plain text file, with a unique protein accession of each line. The path to this file is 
then passed to ``cazy_webscraper get_ncbi_taxs`` via the respective ``--genbank_accessions`` and ``--uniprot_accessions`` flag.

``--genbank_accessions`` and ``--uniprot_accessions`` can be used at the same time to define all 
CAZymes of interest.

.. WARNING::
   ``--genbank_accessions`` and ``--uniprot_accessions`` take president over the filter flags.

   When either ``--genbank_accessions`` or ``--uniprot_accessions`` is used, ``cazy_webscraper get_ncbi_taxs`` will 
   **not** retrieve any CAZymes from the local database matching a set of criteria.

   Therefore, if ``--genbank_accessions`` and ``--classes`` are used, ``cazy_webscraper get_ncbi_taxs`` will ignore 
   the ``--classes`` flag and only taxonomic classifications for the proteins listed in the file provided via 
   the ``--genbank_accessions``.
