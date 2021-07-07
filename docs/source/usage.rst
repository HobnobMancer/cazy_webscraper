=========================
Using ``cazy_webscraper``
=========================

``cazy_webscraper`` can be configured to retrieve user specified data sets from CAZy. The configuration 
applies to the retrieval of protein sequences from GenBank and protein structure files from PDB.

``cazy_webscraper`` can be configured via the **command line** and/or via a **YAML configuration file**.

This page provides a brief summary of each of the options that can be invoked to configure the scraping of CAZy 
using ``cazy_webscraper``. We are currently working on writing additional pages that work through examples of how to invoke 
and combine these options to fully customise the scraping of CAZy.


Configuration via the command line
-----------------------------------

There are no required/positional arguments for the webscraper, therefore the scraper can be enabled 
by simply calling the scraper at the command line in the terminal: 

.. code-block:: bash
  python3 <path_to_cazy_webscraper.py_file>

The ``cazy_webscraper.py`` file is located within the directory ``scraper``. Therefore, if the terminal 
is already pointing at the ``scraper`` directory, the command to invoke ``cazy_webscraper`` is:

.. code-block:: bash
  python3 cazy_webscraper.py

When optional arguments are provided the default behaviour of the scraper will be performed. 
The default behaviour is to:

* Scrape the entire CAZy databases
* Write the resulting database to standard out (STDOUT)
* Not to retrieve subfamilies (members of subfamilies will be retrieved but only their parent family will be listed)


Options configurable at the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following behaviours of the ``cazy_webscraper`` can be configured at the command-line in the terminal:  

* Limit the scraping of CAZy to specific CAZy classes, CAZy families, kingdoms, genuera, species, strains and/or EC numbers.
* Force writing out the database to a a new or existing directory
* Write out a log file of the packages operation
* Not delete content already present in the output directory
* Enable retrieving subfamilies
* Enable verbose logging during the operation of the webscraper


Command line options
^^^^^^^^^^^^^^^^^^^^

The below are listed the flags to be called when invoking ``cazy_webscraper`` in order to configure its operations. 
For example: ``cazy_webscraper -s -v``

..note::
    If an optional command is not given when calling the webscraper the relative default behaviour 
    will be performed.

    The short option is the shortened version of the command, and the long option is the long hand 
    method of writing the command, only one of the other needs be called not both.


.. list-table:: Command line options
   :header-rows: 1

   * - Short option
     - Long option
     - How to use
     - Default behaviour
   * - ``-c``
     - ``--config``
     - Pass a path to a YAML configuration file, to specify CAZy classes, families, kingdoms, genera, species, strains and EC numbers to scrape.
     - No path given, therefore, no configuration provided from a YAML file. Scrape all of CAZy unless other command-line options are provided.
   * -
     - ``--classes``
     - Define classes from which all families are to be scraped. Separate classes with a single comma (e.g. ``GH,PL``)
     - No CAZy classes specified for scraping, therefore, scrape all CAZy classes (unless specific families are provided)
   * - ``-d``
     - ``--database``
     - Path to an existing SQL database produced by ``cazy_webscraper``, this enables adding more data retrieved from CAZy to an existing database.
     - Build a new database.
   * - 
     - ``--ec``
     - Specify specific EC numbers so that only CAZymes annotated with these EC numbers are scraped. Separate EC numbers with a single comma (e.g. ``EC1.2.3.4,EC2.3.4.5``)
     - No EC numbers specified, therefore, retrieve data for CAZymes regardless of their EC number annotations
   * - ``-f``
     - ``--force``
     - Forcing writing out to directory that already exists, if specificed as the output directory. Add the option to the command and do not add anything else, for example: ``python3 cazy_webscraper -o my_dir/ -f``.
     - Force is false, thus if the output directory already exists the output will not be written to it. ``cazy_webscraper`` will raise an error and close.
   * -
     - ``--families``
     - Define CAZy families to be scraped. Separate families with a single comma (e.g. ``GH1,GH2``). Use the proper CAZy nomenclature for family names (e.g. GH1 not gh1)
     - No families specified to be scraped at the command line
   * -
     - ``--genera``
     - Specify specific genera of source organisms of CAZyme to scrape. Separate multiple genera with a comma (e.g. ``Aspergillus,Trichoderma``)
     - No genera specified, so scrapes CAZymes from all genera unless other configuration filters are applied.
   * - 
     - ``--get_pages``
     - Retrieve the HTML code from CAZy for the families that match the configuration criteria provided. The HTML pages are written out to the disk as HTML files, to create a local library of CAZy HTML pages.
     - Do not write out HTML files to disk and instead scrape the protein data from the retrieved HTML pages.
   * - ``-h``
     - ``--help``
     - Displays all command-line options for the webscraper in the terminal, including defining their default behaviour and required additional information.
     - Not to display the help information.
   * - 
     - ``--kingdoms``
     - Specify specific taxonomic kingdoms to retrieve CAZymes only sourced from organisms from these kingdoms. The available kingdoms are: archaea, bacteria, eukaryota, viruses, unclassified (not case sensitive).
     - No taxonomic kingdoms provided, therefore scrape from all kingdoms.
   * - ``-l``
     - ``--log``
     - Enable writing out a log file, logging the operation of the webscraper. Add the option to command followed by desired path of the resulting log file. This the path to the file not the directory to which the log file is to be written.
     - Not to write out a log file of the webscrapers operation.
   * - ``-n``
     - ``--nodelete``
     - Do not delete content in the already existing output directory, this applies for the dataframe, FASTA and protein structure output directories. Simply add this option to the command.
     - Nodelete is false, delete the content already presented in the already existent output directories.
   * - ``-o``
     - ``--output`` 
     - Specify the directory to write the resulting database of CAZymes to. This directory does not already have to exist, if it does not exist ``cazy_webscraper``` will make the directory. Add the option to the command followed by the path to the desired output directory.
     - Write the output database to standard out.
   * - ``-r``
     - ``--retries``
     - Number of times to retry scraping a family or class page if error encountered.
     - 10
   * -
     - ``--scrape_files``
     - Scrape CAZyme data from local HTML files containing HTML code from CAZy. Pass the path to the directory containing the CAZy HTML files.
     - Do not scrape from local HTML pages but call directly to CAZy and scrape the HTML code as it is retrieved.
   * - ``-s``
     - ``--subfamilies``
     - Enable retrieval of subfamilies. If not enabled then the parent CAZy family will be listed for the relevant CAZymes. Simply add the option to the command.
     - Do not retrieve the subfamily annotation. Only the parent CAZy family annotation will be added to the database for applicable CAZymes.
   * - 
     - ``--species``
     - Specify specific species to retrieve CAZymes from. Specifying the species will result in CAZymes from all strains of this species being retrieved. To list multiple species, separate them with a comma (e.g. ``Aspergillus niger,Aspergillus fumigatus``).
     - No specific species specified.
   * - 
     - ``--strains``
     - Specify specific strains of species to retrieve CAZymes from. To list multiple species, separate them with a comma (e.g. ``Aspergillus niger CBS 513.88,Aspergillus fumigatus Af293``).
     - No specific species specified.
   * - 
     - ``--streamline``
     - Specify attributes that are presumed to be the same each time the same CAZyme is parsed from multiple families. The options are: genbank, ec, uniprot and pdb. Any combination can be provided. GenBank refers to non-primary GenBank accessions.
     - Streamline mode not enabled, therefore, for every every CAZyme record, check all its provided data is catalogued into the database.
   * - ``-t``
     - ``--timeout``
     - Specify how long (in seconds) a connection is tried before it is called as timed out.
     - 45
   * - ``-v``
     - ``--verbose``
     - Enable verbose logging of the webscraper. This provides more detailed logging of the progress of the webscrapers operation. Simply add the option to the command.
     - Do not perform verbose logging. Only log if a warning or error is raised.


Basic examples of configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
