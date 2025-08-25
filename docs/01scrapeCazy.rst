==================
How to scrape CAZy
==================

``cazy_webscraper download`` can be used to retrieve user-specified data sets from the CAZy database. The ``cazy_webscraper`` application can be invoked via the command line

To download the entire CAZy dataset, and save the data set to the current working directory with the final name 
``cazy_webscraper_<date>_<time>.db``, use the following command structure:  

.. code-block:: bash
  
   cazy_webscraper download <email> <user_email>

.. NOTE::
   The user email address is a requirement of NCBI. NCBI is queried to identify the currect source organism 
   for a given protein, when multiple source organisms are retrieved from CAZy for a single protein. 
   For more information please see the `NCBI Entrez <https://www.ncbi.nlm.nih.gov/books/NBK25497/>`_ documentation.

.. NOTE::
  Typically, downloading the entire CAZy dataset takes 5-15 minutes, although this is dependent on the amount of available memory.

--------------------
Command line options
--------------------

Listed below are the required and optional command-line options when using ``cazy_webscraper download <email>`` 
to download data from CAZy.

REQUIRED arguments:
  ``email``                 User email address. Required NCBI Entrez - used to get source organsism data. The email address is not stored be
                            cazy_webscraper.

OPTIONAL arguments:

Filtering arguments:
  ``--classes CLASSES``     Classes from which all families are to be scraped. Separate classes by ',' (default: None)
  ``--families FAMILIES``   Families to scrape. Separate families by commas 'GH1,GH2' (case sensitive) (default: None)
  ``--kingdoms KINGDOMS``   Tax Kingdoms to restrict the scrape to (default: None)
  ``--genera GENERA``       Genera to restrict the scrape to (default: None)
  ``--species SPECIES``     Species (written as Genus Species) to restrict the scrape to (default: None)
  ``--strains STRAINS``     Specific strains of organisms to restrict the scrape to (written as Genus Species Strain) (default: None)

Operational arguments:
  ``-o DB_OUTPUT``, ``--db_output DB_OUTPUT``
                            Build a NEW database. Provide the path to build new SQL database (default: None)
  ``-d DATABASE``, ``--database DATABASE``
                            Path to an EXISTING local CAZy SQL database. Add data to this database (default: None)
  ``-s``, ``--subfamilies`` Enable retrieval of subfamilies from CAZy (default: False)
  ``-c config file``, ``--config config file``
                            Path to configuration file. Default: None, scrapes entire database (default: None)
  ``-f``, ``--force``       Force over writting an EXISTING database (default: False)
  ``--cazy_data CAZY_DATA``
                            Path predownloaded CAZy txt file (default: None)
  ``--delete_old_relationships``
                            Delete old GenBank accession - CAZy family relationships (annotations) that are in the local db but are not in CAZy, e.g.
                            when CAZy has moved a protein from one fam to another, delete the old family annotation. (default: False)
  ``--skip_ncbi_tax``       Skip retrieving the latest tax classification from the NCBI Taxonomy db for proteins listed with multiple taxs in CAZy. For
                            these proteins the first taxonomy listed in CAZy is added to the local CAZyme db (default: False)

Utility arguments:
  ``--cache_dir CACHE_DIR``
                            Target path for cache dir to be used instead of default path (default: None)
  ``--cazy_synonyms CAZY_SYNONYMS``
                            Path to JSON file containing CAZy class synoymn names (default: None)
  ``--ncbi_batch_size NCBI_BATCH_SIZE``
                            Number of genbank accessions in each NCBI Taxonomy db batch query (default: 200)
  ``--nodelete_cache``
                            When called, content in the existing cache dir is NOT deleted (default: False)
  ``-r RETRIES``, ``--retries RETRIES``
                            Number of times to retry scraping a family or class page if error encountered (default: 10)
  ``-t TIMEOUT``, ``--timeout TIMEOUT``
                            Connection timeout limit (seconds) (default: 45)

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Defining CAZy families and classes to scrape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``--classes`` and ``--families`` arguments indicate which groups of data will be selected for scraping from CAZy, e.g.

.. code-block:: bash

  cazy_webscraper download <email> --families GH169 -o GH169.db
  cazy_webscraper download <email> --classes AA -o AA.db

will download all CAZymes from the GH169 family, and the AA class, respectively. More than one class or family can be specified, e.g.

.. code-block:: bash

  cazy_webscraper download <email> --families GH169,GH1,GH2,GH3 -o GH_families.db
  cazy_webscraper download <email> --classes AA,CBM -o other_classes.db

and members of distinct families and classes can be selected simultaneously, e.g.

.. code-block:: bash

  cazy_webscraper download <email> --families GH169,GH1,GH2,GH3 --classes AA,CBM -o complex_query.db

.. NOTE::
  CAZy families should be named using the standard CAZy syntax.
  GH1 is **accepted** (case-insensitive). "Glycoside hydrolase 1" is **not** accepted.

^^^^^^^^^^^^^^^^^
Taxonomic Filters
^^^^^^^^^^^^^^^^^

CAZy classifies source organisms under one of 5 kingdoms:

* Archaea
* Bacteria
* Eukaryota
* Viruses
* Unclassified

To restrict the scraping of CAZy to retrieve proteins that are only derived from taxonomic lineages
of interest use the ``-kingdoms``, ``--genera``, ``--species`` and ``--strains`` flags, in any combination:


Scrape by kingdom:

.. code-block:: bash

   cazy_webscraper download --kingdoms bacteria,eukaryota

Scrape by genus:

.. code-block:: bash

   cazy_webscraper download --genera Aspergillus,Trichoderma

Scrape by species (use quotes if there are spaces):

.. code-block:: bash

   cazy_webscraper download --species "Aspergillus niger,Trichoderma reesei"

Scrape by strain:

.. code-block:: bash

   cazy_webscraper download --strains "Aspergillus niger ATCC 1015"

You can combine these filters as needed. For example, to retrieve all CAZymes from viral species, Aspergillus genus, and specific Layia species:

.. code-block:: bash

   cazy_webscraper download \
       --kingdoms viruses \
       --genera Aspergillus \
       --species "Layia carnosa,Layia chrysanthemoides"

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Retrieving Subfamily Annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To store CAZy subfamily annotations in the local CAZyme database, use the ``--subfamilies`` flag:

.. code-block:: bash

   cazy_webscraper download --families GH3 --subfamilies

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Build a new local CAZyme database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default ``cazy_webscraper`` writes out a SQL database file to the current working directory, with a 
filename with the following structure ``cazy_webscraper_<date>_<time>.db``, where the date and time mark 
the time ``cazy_webscraper`` was called.

To specify the location of the output database the ``--db_output`` / ``-o`` option can be used:

.. code-block:: bash

  cazy_webscraper download <email> --families GH169 -o GH169_output.db

will write an SQL database file to ``GH169_output.db``.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Overwrite an existing local CAZyme database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the target output file already exists, ``cazy_webscraper`` by default will not overwrite the existing file and will terminate. To 
overwrite an existing file use the ``--force`` / ``-f`` options:

.. code-block:: bash

  cazy_webscraper download <email> --families GH169 -o GH169_output.db -f

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Adding data to an existing local CAZyme database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you already have an existing CAZy database, then specifying this database with the ``-d`` / ``--database`` option will cause the scraper to use the existing database rather than creating a new one:

.. code-block:: bash

  cazy_webscraper download <email> --families GH169 -d GH169/GH169_output.db

^^^^^^^^^^^^^^^^^^^^
Other Useful Options
^^^^^^^^^^^^^^^^^^^^

- Use ``--help`` for a summary of all options
- Use ``--cazy_data <file>`` to scrape from a downloaded CAZy text file.
- Use ``--log <file>`` to write terminal output to a log file.
- Use ``--verbose`` for detailed logging.
- Use ``--timeout <seconds>`` to set the connection timeout.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Configuration using a YAML file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All command-line options to control CAZy scraping can be provided instead *via* a YAML configuration file. This supports reproducible documentation of ``cazy_webscraper`` usage.

An template YAML file is provided in the ``cazy_webscraper`` repository (``scraper/scraper_config.yaml``):

.. code-block:: yaml

  # Under 'classes' list class from which all proteins will retrieved
  # Under each families respective name, list the specific families/subfamilies to be scraped
  # Write the FULL family name, e.g. 'GH1', NOT only its number, e.g. '1'
  # To list multiple families, each familiy must be on a new line starting indented once
  # relative to the parent class name, and the name written within quotation marks.
  # For more information on writing lists in Yaml please see:
  # https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html 
  classes:  # classes from which all proteins will be retrieved
  Glycoside Hydrolases (GHs):
  GlycosylTransferases (GTs):
  Polysaccharide Lyases (PLs):
    - "PL28"
  Carbohydrate Esterases (CEs):
  Auxiliary Activities (AAs):
  Carbohydrate-Binding Modules (CBMs):
  genera:  # list genera to be scraped
   - "Trichoderma"
  species:  # list species, this will scrape all strains under the species
  strains:  # list specific strains to be scraped
  kingdoms:  # Archaea, Bacteria, Eukaryota, Viruses, Unclassified
   - "Bacteria"

.. ATTENTION::
  The YAML configuration file must contain all tags/headings indicated in the example configuration file found in the repository:

  * classes
  * Glycoside Hydrolases (GHs)
  * GlycosylTransferases (GTs)
  * Polysaccharide Lyases (PLs)
  * Carbohydrate Esterases (CEs)
  * Auxiliary Activities (AAs)
  * Carbohydrate-Binding Modules (CBMs)
  * genera
  * species
  * strains
  * kingoms

Each value in the YAML mappings for these arguments must be listed on a separate line, indented by 4 spaces, and the class name encapsulated with single or double quotation marks. For example:

.. code-block:: yaml

    classes:
        - "GT"
        - "pl"
    Glycoside Hydrolases (GHs):
        - "GH1"
        - "GH2"
