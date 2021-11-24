=========================
Using ``cazy_webscraper``
=========================

``cazy_webscraper`` can be used to retrieve user-specified data sets from the CAZy database. The ``cazy_webscraper`` application can be invoked *via* the command line

----------------------
Quick Start
----------------------

To download the entire CAZy dataset, and save the data set to the current working directory with the final name 
``cazy_webscraper_<date>_<time>.db``, use the following command structure:  

.. code-block:: bash
   cazy_webscraper <user_email>

.. NOTE::
   The user email address is a requirement of NCBI. NCBI is queried to identify the currect source organism 
   for a given protein, when multiple source organisms are retrieved from CAZy for a single protein.

--------------------
Command line options
--------------------


``email`` - **REQUIRED** User email address. This is required by NCBI Entrez for querying the Entrez server.

``--cache_dir`` - Path to cache dir to be used instead of default cache dir path.

``--cazy_data`` - Path to a txt file downloaded from CAZy containing a CAZy database dump

``--cazy_synonyms`` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

``--config``, ``-c`` - Path to a configuration YAML file. Default: None.

``--citation``, ``-C`` - Print the `cazy_webscraper` citation. When called, the program terminates after printng the citation and CAZy is **not** scraped.

``--classes`` - list of classes from which all families are to be scrape.

``--database``, ``-d`` - Path to an **existings** local CAZyme database to add newly scraped too. Default: None.

``--db_output``, ``-o`` - Path to write out a **new** local CAZyme database to. Include the name of the new database, including the `.db` extension. Default: None.

.. WARNING::
  **Do not use ``--db_output`` and ``--database`` at the same time.**

.. NOTE::
  If ``--db_output`` **and** ``--database`` are **not** called,
  ``cazy_webscraper`` will write out a local CAZyme database to the cwd with the standardised name ``cazy_webscraper_<date>_<time>.db``

``--delete_old_relationships`` - Detele old CAZy family annotations of GenBank accessions. These are CAZy family annotations of a given GenBank accession are in the local database but the accession is not longer associated with those CAZy families, so delete old accession-family relationships.

``--families`` - List of CAZy (sub)families to scrape.

``--force``, ``-f`` - force overwriting existing output file. Default: False.

.. WARNING::
  If a specified output directory already exists, if ``--force`` is not called, ``cazy_webscraper`` 
  will not overwrite the output and terminate.

``--genera`` - List of genera to restrict the scrape to. Default: None, filter not applied to scrape.

``--log``, ``-l`` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

``--nodelete``, ``-n`` - When called, content in the existing output dir will **not** be deleted. Default: False (existing content is deleted).

``--nodelete_cache`` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

``--nodelete_log`` - When called, content in the existing log dir will **not** be deleted. Default: False (existing content is deleted).

``--retries``, ``-r`` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

``--sql_echo`` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

``--subfamilies``, ``-s`` - Enable retrival of CAZy subfamilies, otherwise **only** CAZy family annotations will be retrieved. Default: False.

``--species`` - List of species written as Genus Species) to restrict the scraping of CAZymes to. CAZymes will be retrieved for **all** strains of each given species.

``--strains`` - List of specific species strains to restrict the scraping of CAZymes to.

``--timeout``, ``-t`` - Connection timout limit (seconds). Default: 45.

``--validate``, - Retrieve CAZy family population sizes from the CAZy website and check against the number of family members added to the local CAZyme database, as a method for validating the complete retrieval of CAZy data.

``--verbose``, ``-v`` - Enable verbose logging. This does **not** set the SQLite engine ``echo`` parameter to True. Default: False.

``--version``, ``-V`` - Print ``cazy_webscraper`` version number. When called and the version number is printed, ``cazy_webscraper`` is immediately terminated.

-----------
Basic Usage
-----------

The command-line options listed above can be used in combination to customise the scraping of CAZy. Some options (e.g. ``--families`` and ``--classes``) define the broad group of data that will be scraped, others (e.g. ``--species``) are used to filter and fine-tune the data that is scraped.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Defining CAZy families and classes to scrape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The 'definition' arguments (e.g. ``--classes`` and ``--families``) indicate which groups of data will be selected for scraping from CAZy, e.g.

.. code-block:: bash

  cazy_webscraper --families GH169 -o GH169.db
  cazy_webscraper --classes AA -o AA.db

will download all CAZymes from the GH169 family, and the AA class, respectively. More than one class or family can be specified, e.g.

.. code-block:: bash

  cazy_webscraper --families GH169,GH1,GH2,GH3 -o GH_families.db
  cazy_webscraper --classes AA,CBM -o other_classes.db

and members of distinct families and classes can be selected simultaneously, e.g.

.. code-block:: bash

  cazy_webscraper --families GH169,GH1,GH2,GH3 --classes AA,CBM -o complex_query.db

.. NOTE::
  CAZy families should be named using the standard CAZy syntax.
  GH1 is **accepted**.
  gh1 and Glycoside hydrolase 1 are **note** accepted.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifying output data location
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default ``cazy_webscraper`` writes out a SQL database file to the current working directory, with a 
filename with the following structure ``cazy_webscraper_<date>_<time>.db``, where the date and time mark 
the time ``cazy_webscraper`` was called.

To specify the location of the output database the ``--db_output``/``-o`` option can be used:

.. code-block:: bash

  cazy_webscraper --families GH169 -o GH169_output.db

will write an SQL database file to ``GH169_output.db``.

If the target output file already exists, ``cazy_webscraper`` by default will not overwrite the existing file and will terminate. To 
overwrite an existing file use the ``--force``/``-f`` options:

.. code-block:: bash

  cazy_webscraper --families GH169 -o GH169_output.db -f

A multi-layered path can be provided to ``cazy_webscraper``. If any of the parent directories for the target 
output path do not exist, ``cazy_webscraper`` will build the necessary output direcotires. In the following command if 
the ``cazy`` and ``families`` directories do not exist, ``cazy_webscraper`` will build these directories:

.. code-block:: bash

  cazy_webscraper --families GH169 -o cazy/families/GH169_output.db 

If any of the output directories exist, by default, ``cazy_webscraper`` will terminate. To write to an existing output 
directory use the ``--force``/``-f`` options:

.. code-block:: bash

  cazy_webscraper --families GH169 -o GH169_output.db -f

By default ``cazy_webscraper`` will delete the existing content in the existing output files. To not delete the content 
in the existing output directories use the ``--nodelete``/``-n``:

.. code-block:: bash

  cazy_webscraper --families GH169 -o GH169_output.db -f -n

If you already have an existing CAZy database, then specifying this database with the ``-d``/``--database`` option will cause the scraper to use the existing database rather than creating a new one:

.. code-block:: bash

  cazy_webscraper --families GH169 -d GH169/GH169_output.db

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Filtering CAZy families and classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Options that apply a *filter* to restrict which CAZymes from a class or familiy are scraped from CAZy (e.g.  ``--families`` and ``--species``) may be applied in combination. For example:

.. code-block:: bash

  cazy_webscraper --families GH169 \
      --species "Escherichia coli" \
      -o GH169_speciesEscherichia_coli.db

will download only the CAZymes in the GH169 family that are from the species *Escherichia coli*. The command:

.. code-block:: bash

  cazy_webscraper --families PL14,PL15,PL16 \
      -o PL14_ec1.2.3.4_kingdomBacteria

will download only CAZymes in the PL14, PL15 and PL16 families that are from the kingdom *Bacteria*.

.. NOTE::
  ``cazy_webscraper`` input options can also be specified in a **YAML configuration file**, to enable transparency and reproducibility.

-------------------------------
Configuration using a YAML file
-------------------------------

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


^^^^^^^^^^^^^^^^^^^^^^^^^
Synonyms for CAZy classes
^^^^^^^^^^^^^^^^^^^^^^^^^

A number of synonyms may be provided for CAZy classes, e.g. both "GH" and "Glycoside-Hydrolases" are accepted as synonyms for "Glycoside Hydrolases (GHs)" (the name recorded at CAZy). These alternatives are defined in the ``cazy_webscraper`` repository, in the file ``scraper/utilities/parse_configuration/cazy_dictionary.json``.

-------------------------
Scraping CAZy subfamilies
-------------------------

``cazy_webscraper`` can scrape CAZy subfamilies, using the standard CAZy notation for subfamilies 
(e.g. ``GH3_1``).

.. NOTE::
   If any subfamilies are specified for download/scraping in the YAML file, the command line argument ``--subfamilies`` must be used.

If a parent CAZy family is listed in the configuration file and ``--subfamilies`` is enabled at the command-line, all proteins catalogued under the named family and its subfamilies will be retrieved.
