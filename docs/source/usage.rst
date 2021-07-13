=========================
Using ``cazy_webscraper``
=========================

``cazy_webscraper`` can be used to retrieve user-specified data sets from the CAZy database. The ``cazy_webscraper`` application can be invoked *via* the command line

-------------
Example Usage
-------------

To download the single CAZy family GH169, use the command:

.. code-block:: bash

  cazy_webscraper --families GH169 -o GH169

This will create a new directory ``GH169`` in the current working directory, and will download all CAZy entries in the GH169 family to a new SQLite3 database in that directory.

This page provides a brief summary of command-line options for ``cazy_webscraper`` that control the retrieval of data sets from the CAZy database, including:

* Retrieve only specified CAZy classes and families or subfamilies
* Retrieve only CAZymes from specified taxonomic kingdoms, genera, species, or strains
* Recover only CAZymes with specified EC numbers
* Local SQLite3 database path
* Verbosity level and logging options

--------------------
Command line options
--------------------

.. list-table:: Command line options
   :header-rows: 1

   * - Short option
     - Long option
     - Action
     - Default
   * - ``-c``
     - ``--config``
     - path to a YAML configuration file
     - do not use YAML configuration file
   * -
     - ``--classes``
     - define CAZy classes to be retrieved (comma-separated list for multiple classes)
     - no specified value
   * - ``-d``
     - ``--database``
     - path to SQLite3 database
     - create a new database with default name
   * - 
     - ``--ec``
     - define EC numbers to filter CAZyme data (comma-separated list for multiple values)
     - no specified value
   * - ``-f``
     - ``--force``
     - force overwriting of existing output
     - do not force overwrite
   * -
     - ``--families``
     - define CAZy families to be retrieved (comma-separated list for multiple families)
     - no specified value
   * -
     - ``--genera``
     - filter CAZyme data on taxonomic genus (comma-separated list for multiple values)
     - do not filter on genus
   * - 
     - ``--get_pages``
     - retrieve HTML from CAZy for specified CAZy families and write to disk
     - do not retrieve HTML to disk
   * - ``-h``
     - ``--help``
     - display command line options
     -  
   * - 
     - ``--kingdoms``
     - filter CAZyme data on taxonomic kingdom (comma-separated list for multiple values)
     - do not filter on kingdom
   * - ``-l``
     - ``--log``
     - path to log file
     - do not write log file
   * - ``-n``
     - ``--nodelete``
     - do not delete ("clobber") existing output when overwriting
     - delete existing output when overwriting
   * - ``-o``
     - ``--output`` 
     - path to output database
     - write to STDOUT
   * - ``-r``
     - ``--retries``
     - number of times to retry CAZy web requests.
     - 10
   * -
     - ``--scrape_files``
     - path to local CAZy HTML files; data will be scraped from these files instead of CAZy website
     - do not use local CAZy HTML files
   * - ``-s``
     - ``--subfamilies``
     - define CAZy subfamilies to be retrieved (comma-separated list for multiple families)
     - no specified value
   * - 
     - ``--species``
     - filter CAZyme data on taxonomic species (comma-separated list for multiple values)
     - do not filter on species
   * - 
     - ``--strains``
     - filter CAZyme data on strain (comma-separated list for multiple values)
     - do not filter on strain
   * - 
     - ``--streamline``
     - override CAZy metadata for each recovered record
     - do not override CAZy metadata
   * - ``-t``
     - ``--timeout``
     - wait time before CAZy web connection is considered timed out
     - 45
   * - ``-v``
     - ``--verbose``
     - enable verbose logging and output
     - standard logging/output level

-----------
Basic Usage
-----------

The command-line options listed above can be used in combination to customise the scraping of CAZy. Some options (e.g. ``--families`` and ``--classes``) define the broad group of data that will be scraped, others (e.g. ``--species``) are used to filter and fine-tune the data that is scraped.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Defining CAZy families and classes to scrape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The 'definition' arguments (e.g. ``--classes`` and ``--families``) indicate which groups of data will be selected for scraping from CAZy, e.g.

.. code-block:: bash

  cazy_webscraper --families GH169 -o GH169
  cazy_webscraper --classes AA -o AA

will download all CAZymes from the GH169 family, and the AA class, respectively. More than one class or family can be specified, e.g.

.. code-block:: bash

  cazy_webscraper --families GH169,GH1,GH2,GH3 -o GH_families
  cazy_webscraper --classes AA,CBM -o other_classes

and members of distinct families and classes can be selected simultaneously, e.g.

.. code-block:: bash

  cazy_webscraper --families GH169,GH1,GH2,GH3 --classes AA,CBM -o complex_query

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Filtering CAZy families and classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Options that apply a *filter* to restrict which CAZymes from a class or familiy are scraped from CAZy (e.g. ``--species`` and ``--ec``) may be applied in combination. For example:

.. code-block:: bash

  cazy_webscraper --families GH169 \
      --ec 1.1.1.1 --species "Escherichia coli" \
      -o GH169_ec1.1.1.1_speciesEscherichia_coli

will download only the CAZymes in the GH169 family that have EC number 1.1.1.1 *and* are from the species *Escherichia coli*. The command:

.. code-block:: bash

  cazy_webscraper --families PL14 \
      --ec 1.2.3.4 --kingdoms bacteria \
      -o PL14_ec1.2.3.4_kingdomBacteria

will download only CAZymes in the PL14 familiy that have EC number 1.2.3.4 *and* are from the kingdom *Bacteria*.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifying output data location
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To specify the location of the output database and log files, the ``--output``/``-o`` option can be used:

.. code-block:: bash

  cazy_webscraper --families GH169 -o GH169_output

will write output to the directory ``GH169_output``, and create a new CAZy database in that directory.

If you already have an existing CAZy output directory with a database, then specifying this database with the ``-d``/``--database`` option will cause the scraper to use the existing database rather than creating a new one:

.. code-block:: bash

  cazy_webscraper --families GH169 -d GH169_output/cazy.db

To write output to an existing directory without deleting the content already present, use the ``--force``/``-f`` and ``--nodelete``/``-n`` options:

.. code-block:: bash

  cazy_webscraper --families GH169 -d GH169_output -f -n


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



A ``cazy_dictionary.json`` has been created and packaged within the ``cazy_webscraper`` 
(the specific location is ``./scraper/file_io/cazy_dictionary.json``, where '.' is the directory 
where the webscraper is installed). This allows users to use a variety of synonoms for the CAZy 
classes, for example both "GH" and "Glycoside-Hydrolases" are accepted as synonoms for 
"Glycoside Hydrolases (GHs)". Additionally, the retrieval of CAZy classes from the configuration 
file is **not** case sensitive, therefore, both "gh" and "GH" are excepted. The excepted class 
synonoms have beeen written out in a json file to enale easy editing of this file if additional 
accepted synonoms are to be added, of it a new CAZy class is defined then this class only needs 
to be added to the json file, without needing to modify the entire webscraper. 

If you having issues with the scraper retrieving the list of CAZy classes that are written under 
'classes' in the configuration file, please check the dictionary first to see the full list of 
accepted synonoms. If you are comfortable modifying json files then feel free to add your own 
synonoms to the dictionary.

Each class must be listed on a separate line, indented by 4 spaces, and the class name encapsulated 
with single or double quotation marks. For example:

.. code-block:: yaml

    classes:
        - "GH"
        - "pl"


Specifying specific families to scrape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under the each of the class names listed in the configuration file, list the names of specific 
**families** to be scraped from that class. The respective classes of the specificed families do 
**not** need to be added to the 'classes' list.

Write the true name of the family not only it's number, for example **GH1** is excepted by **1** is 
not. Name families using the standard CAZy nomenclature, such as **"GT2"** and 
**NOT "GlycosylTransferases_2"**. Additionally, use the standard CAZy notation for subfamilies 
(**GH3_1**).

.. warning::
   If any subfamilies are listed within the configuration file, the retrieval of subfamilies 
   **must** be enabled at the command line uisng ``--subfamilies``.

Each family must be listed on a separate line and the name surrounded by double or single quotation 
marks. For example:

.. code-block:: yaml

    Glycoside Hydrolases (GHs):
        - "GH1"
        - "GH2"


Configuration when scraping subfamilies
---------------------------------------

If any subfamilies are listed within the configuration file, the retrieval of subfamilies **must** 
be enabled at the command line uisng ``--subfamilies``.

If the parent family, e.g GH3, is listed in the configuration file and ``--subfamilies`` is enabled, 
all proteins catalogued under GH3 and its subfamilies will be retrieved. This is to save time 
having to write out all the subfamilies for a given CAZy family. The scraper will remove any 
duplicate proteins automatically.


An example configuration file
-----------------------------

A blank configuration file is packaged within ``cazy_webscraper``, within the ``scraper`` directory, 
called ``scraper_config.yaml``. This configuration file contains comments to assit filling in the 
file correctly. A new configuration file with any given name can be created and used. However, 
it **must** be a Yaml file and it **must** use the same headings/tags as used in the configuration 
file ``scraper_config.yaml``.Please find more information on writing lists in Yaml files 
[here](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html).

Below is an example of how the configuration file may look.

.. code-block:: yaml

    classes:
        - "AA"
    Glycoside Hydrolases (GHs):
        - "GH1"
        - "GH3"
    GlycosylTransferases (GTs):
    Polysaccharide Lyases (PLs):
        - "PL9"
    Carbohydrate Esterases (CEs):
    Auxiliary Activities (AAs):
    Carbohydrate-Binding Modules (CBMs):


..note::
    Indentations consist of 4 spaces.
