=========================================
Retrieving Protein Sequences from GenBank
=========================================

``cazy_webscraper`` can be used to retrieve protein amino acid sequences from NCBI GenBank for user-specified data sets of CAZymes 
in the local CAZymes database, and store these sequences in an existing local CAZyme database.

The retrieval of data from NCBI is performed by using the **BioPython** `Bio.entrez <https://biopython.org/docs/1.75/api/Bio.Entrez.html>_` module [Cock *et al.*, 2009].

    Cock, P. J. A, Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A. _et al._ (2009) 'Biopython: freely available Python tools for computaitonal molecular biology and bioinformatics', _Bioinformatics_, 25(11), pp. 1422-3.

.. note::
    For specific information of the ``Bio.entrez`` module please see the 
    `entrez documentation <https://biopython.org/docs/1.75/api/Bio.Entrez.html>`_.


-----------
Quick Start
-----------

To download protein sequences for all CAZymes in the local CAZyme database, and write them to the local CAZyme database, 
use the following command structure:

.. code-block:: console

    cazy_webscraper get_genbank_seqs <path to local CAZyme db> <user email address>


--------------------
Command line options
--------------------


REQUIRED arguments:
    database              Path to local CAZy database
    email                 User email address, requirement of NCBI-Entrez. The address is not stored by cazy_webscraper

OPTIONAL arguments:

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
                                                Path to a text file containing a list of GenBank accessions to retrieve data for. Note, accessions from the local DB will not
                                                be retrieved, only the accessions in this this file will be used. (default: None)
    ``--seq_dict`` SEQ_DICT   Path to a JSON file, keyed by GenBank accessions and valued by protein sequence Add seqs in file to the local CAZyme
                                                database. This is created by cazy_webscraper during get_gbk_seqs (default: None)
    ``--seq_file`` SEQ_FILE   Path to a FASTA file of protein sequences Add seqs in file to the local CAZyme database (default: None)
    ``-F``, ``--file_only``       Only add seqs provided via JSON and/or FASTA file. Do not retrieve data from NCBI (default: False)
    ``--seq_update``          Enable overwriting sequences in the database if the retrieved sequence is different (default: False)

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

--------------------------------
GenBank sequence retrieval cache
--------------------------------

``cazy_webscraper get_genbank_seqs`` produces three cache files, which are written to the cache dir:
1. ``no_seq_retrieved.txt`` which lists the GenBank accessions for which no sequence could be retrieved from GenBank
2. ``seq_retrieved.txt`` which list GenBank accessions for which a sequence was retrieved from GenBank
3. JSON file keyed by GenBank accessions and valued by the retrieved protein sequence

-----------
Basic Usage
-----------

The command-line options listed above can be used in combination to customise the scraping of CAZy. Some options (e.g. ``--families`` and ``--classes``) 
define the broad group of data that will be scraped, others (e.g. ``--species``) are used to filter and fine-tune the data that is scraped.

The ``--classes``, ``--families``, ``--kingdoms``, ``--genera``, ``--species``, and ``--strains`` filteres are applied 
in the exactly same for retrieving data from CAZy as retrieving protein sequences from GenBank and protein data from UniProt. Examples of using these flags 
can be found in the `tutorial <https://cazy-webscraper.readthedocs.io/en/latest/genbanktutorial.html>`_.

The ``--seq_update`` flag is used in the same way for retrieving protein sequences from UniProt and GenBank.

.. NOTE::
    To retrieve data for members of specific CAZy subfamilies, list the subfamilies after the ``--families`` 
    flag.

------------------------
Updating local sequences
------------------------

When using ``--sequence`` flag, ``cazy_webscraper get_genbank_seqs`` will only add *new* protein sequences to the database, i.e.
it will only add protein sequences to records that do not have a sequence. Therefore, if a protein
already has a sequence in the local database, this sequence is **not** overwritten.

To update existing protein sequences in the local CAZyme database, use the ``--seq_update`` flag:

.. code-block:: console

    cazy_webscraper get_genbank_seqs my_databases/cazy_db.db --seq_update

``cazy_webscraper get_genbank_seqs`` will overwrite existing protein sequences in the local database *if* a newer version 
of the sequence is retrieved from UniProt. This is checked by comparing the 'last modified date' of the 
protein sequence in the local database against the sequence retrieved from UniProt.

-------------------------------------------------------------------
Retrieving protein sequences for specific CAZy classes and families
-------------------------------------------------------------------

The ``--classes`` and ``--families`` flags from scraping data from CAZy are applied in the extact same way 
for retrieving protein sequences from GenBanks.

For example, if you want to retrieve protein sequences for all CAZymes from Glycoside Hydrolase and Carbohydrate Esterases then use the command:

.. code-block:: bash

   cazy_webscraper get_genbank_seqs cazy/cazyme.db dummy.email@domain.co.uk --classes GH,CE

To retrieve protein sequences for all proteins in PL1, PL2 and PL3 in the local CAZyme database use the 
following command:

.. code-block:: bash

   cazy_webscraper get_genbank_seqs cazy/cazyme.db dummy.email@domain.co.uk --families PL1,PL2,PL3

As with scraping data from CAZy, the ``--classes`` and ``--families`` flags can be combined. To retrieve 
protein sequences for all CAZymes in PL1, PL2, PL3 and *all* of GH and CE both:

.. code-block:: bash

   cazy_webscraper get_genbank_seqs cazy/cazyme.db dummy.email@domain.co.uk --families PL1,PL2,PL3 --classes GH,CE

--------------------------
Applying taxonomic filters
--------------------------


The ``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` flags can be used to refine the dataset 
of proteins to retrieve protein sequences by taxonomy. These flags are applied in the exact same way as they 
are used for the scraping of data from CAZy. Only proteins in the local CAZyme database and matching at least on of the provided taxonomy 
criteria will have protein sequences retrieved from GenBank and added to the local CAZyme datbase.

For example, if you want to retrieve protein sequences for all CAZymes in a local CAZyme database from bacterial and eukaryotic species then use the command 

.. code-block:: bash

   cazy_webscraper get_genbank_seqs cazy/cazyme.db dummy.email@domain.co.uk \
       --kingdoms bacteria,eukaryota


You can combine any combination of the optional flags, including combining the taxonomic filters. For example,
you may wish to retrieve protein sequences for all CAZymes in a local CAZyme database that are derived from all viral species, Aspergillus species, Layia carnosa, Layia chrysanthemoides, Trichoderma reesei QM6a and 
Trichoderma reesei QM9414. To do this we would combine the respective flags for a single ``cazy_webscraper get_genbank_seqs`` command. The command 
we would use would be:

.. code-block:: bash

   cazy_webscraper get_genbank_seqs cazy/cazyme.db dummy.email@domain.co.uk \
       --kingdoms viruses \
       --genera Aspergillus \
       --species Layia carnosa,Layia chrysanthemoides \
       --strains Trichoderma reesei QM6a,Trichoderma reesei QM9414

.. warning::
   Use the standard scientific name formating. Captialise the first letter of *genus* and write a lower 
   case letter for the first letter of the species.

   Aspergillus niger is **correct**

   asepergillus niger is **incorrect**

   ASPERGILLUS NIGER is **incorrect**

   -------------------------
Applying EC number filter
-------------------------

The retrieval of protein sequences from GenBank can also be limited to proteins in a local CAZyme database that are
annotated with specific EC numbers.


Having previously retrieved EC number annotations and added them to the local CAZyme database, you  may 
wish to retrieve protein sequences for CAZymes annotated with specific EC numbers. To do this add the 
``--ec_filter`` flag to the command, follwed by a list of EC numbers.

``--ec_filter`` is based upon EC number annotations stored within the local CAZyme database.

.. code-block:: bash
   
   cazy_webscraper get_genbank_seqs cazy/cazyme.db dummy.email@domain.co.uk \
       --ec_filter "EC1.2.3.4,EC2.3.4.5"

Provide complete EC numbers. Both dashes ('-') and asterixes ('*') are accepted for missing digits in EC numbers. However,
If using dashes to represent missing digits in EC numbers, it is recommended to bookend the entire list in quotation marks.

The 'EC' prefix is not necessary. EC1.2.3.4 and 1.2.3.4 are accepted.

.. NOTE::
    EC1.2.3.- and EC1.2.3.* are accepted.
    EC1.2.3. and EC 1.2.3 are **not** accepted.

.. NOTE::
    ``cazy_webscraper`` will retrieve the specified UniProt data for all proteins in the local CAZyme 
    database that are annotated with **at least one** of the given EC numbers. Therefore, if multiple 
    EC numbers are given this **does not mean** only CAZymes annotaed with all provided EC numbers will have data retrieved
    from UniProt for them.


---------------------
Combining all filters
---------------------

The ``--classes``, ``--families``, ``--ec_filter``, ``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` flags can 
be used in any combination to define a specific subset of proteins in the local CAZyme database for whom
protein sequences from GenBank will be retrieved.

Below we run through 3 example commands of combining these flags, and the resulting behaviour.

**Example 1:**
To retrieve protein sequences for all CAZymes in GH, GT, CE1, CE5 and CE8, and which are derived from baceterial species we use the command:

.. code-block:: bash

   cazy_webscraper get_genbank_seqs cazy/cazyme.db dummy.email@domain.co.uk --classes GH,CE --families CE1,CE5,CE8 --kingdoms bacteria


**Example 2:**
To protein sequences for all CAZymes in GH and which are derived from *Aspegillus* and *Trichoderma* species we use the command:

.. code-block:: bash

   cazy_webscraper get_genbank_seqs cazy/cazyme.db dummy.email@domain.co.uk -classes GH --genera Aspegillus,Trichoderma


**Example 3:**
To retrieve protein sequences for all CAZymes in GH,CE and CBM which are derived from baceterial species and are annotated with at least one of 
EC3.2.1.23, EC3.2.1.37 and EC3.2.1.85, we use the command:

.. code-block:: bash

   cazy_webscraper get_genbank_seqs cazy/cazyme.db dummy.email@domain.co.uk --ec --sequences --classes GH,CE,CBM --kingdoms bacteria --ec_filter "3.2.1.23,3.2.1.37,3.2.1.85"


------------------------------
Providing a list of accessions
------------------------------

Instead of retrieving protein sequences for all CAZymes matching a defined set of criteria, 
``cazy_webscraper get_genbank_seqs`` can retrieve protein sequences a set of CAZymes defined by their 
GenBank and/or UniProt accession.

The flag ``--genbank_accessions`` can be used to provide ``cazy_webscraper get_genbank_seqs`` a list of GenBank accessions 
to identify the specific set of CAZymes to retrieve protein sequences for.

The flag ``--uniprot_accessions`` can be used to provide ``cazy_webscraper get_genbank_seqs`` a list of UniProt accessions 
to identify the specific set of CAZymes to retrieve protein sequences for.

In both instances (for ``--genbank_accessions`` and ``--uniprot_accessions``) the list of respective accessions 
are provided via a plain text file, with a unique protein accession of each line. The path to this file is 
then passed to ``cazy_webscraper get_genbank_seqs`` via the respective ``--genbank_accessions`` and ``--uniprot_accessions`` flag.

``--genbank_accessions`` and ``--uniprot_accessions`` can be used at the same time to define all 
CAZymes of interest.

.. WARNING::
   ``--genbank_accessions`` and ``--uniprot_accessions`` take president over the filter flags.

   When either ``--genbank_accessions`` or ``--uniprot_accessions`` is used, ``cazy_webscraper get_genbank_seqs`` will 
   **not** retrieve any CAZymes from the local database matching a set of criteria.

   Therefore, if ``--genbank_accessions`` and ``--classes`` are used, ``cazy_webscraper get_genbank_seqs`` will ignore 
   the ``--classes`` flag and only retrieve protein sequences for the proteins listed in the file provided via 
   the ``--genbank_accessions``.


-------------------------------
Providing sequences from a file
-------------------------------

While ``cazy_webscraper get_genbank_seqs`` is retrieving protein sequences from NCBI, the retrieved protein sequences 
are written to a FASTA file in the cache directory.

To add sequences from a cached FASTA file (e.g. to continue a download that was previously interrupted) and/or 
add GenBank sequences from a previous download (e.g. by a colleage), use the ``--seq_file`` flag followed by 
the path to the FASTA containing the protein sequences to be added to the database. The ID for 
each sequence **must** be the NCBI protein version accession.

``cazy_webscraper get_genbank_seqs`` also generates a JSON file of the cached sequences. To add sequences from the 
cached JSON file to the local CAZyme database, use the ``--seq_dict`` flag followed by the path to the 
JSON file.

By default ``cazy_webscraper get_genbank_seqs`` will add sequences retrieved from the FASTA and/or JSON file **and** will retrieve 
protein sequences from NCBI for proteins matching the provided criteria to define proteins of interest. 

To add **only** the sequences from a FASTA and/or JSON file, and **not** download any sequences from NCBI, use 
the ``--file_only`` flag.
