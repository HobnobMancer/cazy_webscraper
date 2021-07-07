=========================
Using ``cazy_webscraper``
=========================

``cazy_webscraper`` can be used to retrieve user-specified data sets from the CAZy database. The ``cazy_webscraper`` application can be invoked _via_ the command line

.. NOTE::
  ``cazy_webscraper`` options can be specified in a **YAML configuration file**, to enable transparency and reproducibility.

This page provides a brief summary of command-line options for ``cazy_webscraper`` that control the retrieval of data sets from the CAZy database, including:

* Retrieve only specified CAZy classes and families or subfamilies
* Retrieve only CAZymes from taxonomic kingdoms, genera, species, or strains
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

The command-line options listed above can be used in any combination to customise the scraping of CAZy. The options that apply a 'filter' 
to restrict which CAZymes are scraped from CAZy are applied in combination. For example, if the ``--families`` option and ``--ec`` option are called then 
only CAZymes from the specified families **and** annotated with the listed EC numbers will be retrieved.

Below are some example commands for invoking the ``cazy_webscraper`` to help demonstrate how to configure the webscraper at the command line.

1. Writing the output to the directory 'my_output' and enabling retrieval of subfamilies:  
``python3 cazy_webscraper.py -o my_output -s``

2. Retrieving all CAZymes derived from bacteria and annotated with the EC numbers EC1.2.3.4 or EC1.5.3.4
``python3 cazy_webscraper.py --kingdoms bacteria --ec EC1.2.3.4,EC``

3. Writing the output to an existing directory but do not delete the content already present in the directory:  
``python3 cazy_webscraper.py --output docs/my_output -f -n``

4. Write out the data retrieved from CAZy to an existing database, and only retrieve data for CAZymes derived from Aspergiulls species from families GH13, GH15 and PL9, and all CE familes:  
``python3 cazy_webscraper.py -d docs/my_cazy_database/cazy_scrape_2021-04-27--11-54-58.db --genera Aspergillus --families GH13,GH15,PL9 --classes CE``


Configuration via a YAML file
------------------------------

Using a configuration files produces reproducible documentation of how you used ``cazy_webscraper`` -- which is an essential part of all bioinformatic research.

An example/template YAML file is provided within the repository of the webscraper, located at: 
``./scraper/scraper_config.yaml``. A configuration YAML file must contain the same tags/headings as 
the example configuration file found in the repository. The headings are:

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


Specifying specific classes to scrape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under the **classes** heading list any classes to be scrapped. For classes listed under 'classes', 
all proteins catalogued under that class will be retrieved, **unless** specific families have been 
listed under the respective classes heading in the configuration file. Then scraping only the 
specific families takes precident and the entire class is not scraped. _If you believe this should 
be changed please raise an issue. It is invisioned that very few users would want to simultanious 
scrape an entire class and also scrape only specific families from that same class._

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
