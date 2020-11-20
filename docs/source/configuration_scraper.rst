===========================================
Configuring the CAZy webscraper
===========================================

The CAZy webscraper is configurable to restrict the scraping to specific classes and/or families, 
and enable retrieval of subfamilies in order to meet your unique needs.

The webscraper is configured via the command line and via a YAML configuration file.


Configuration via the command line
-----------------------------------

There are no required/positional arguments for the webscraper. If no optional arguments are provided 
the default behaviour of the scraper will be performed. The default behaviour is to:

* scrape the entire CAZy database
* retrieval of subfamilies is not enabled so only the parent CAZy family will be listed those CAZymes classified under subfamilies
* the resulting dataframe of CAZyme data will not be split, creating a single dataframe containing all scraped CAZymes

Command line options
^^^^^^^^^^^^^^^^^^^^

* ``-c``, ``--config``: Path to the configuration file. Default: None, scrapes entire CAZy database
* ``-d``, ``--data_split``: Choices: None, class, family. Default: None, not to split data when written out; ``class`` will create a single dataframe per class; ``family`` will create a single dataframe per family
* ``-f``, ``--force``: (True/False). Force over writing in output directory if specified output directory already exists. Default: False, does not over write in already exising output dictory. 
* ``-l``, ``--log``: Path to write out a logger file. Default: None. Logger messages will be written out to the terminal and out to the specified file.
* ``-n``, ``--nodelete``: (True or False). Do not delete content in already exisiting output directory. Default: False, will delete content in already existing output directory. If set to true then the content in the output directory will not be deleted first before writing out output from the scrape.
* ``-o``, ``--output``: Path to output DIRECTORY for all output to be written to. Default: STDOUT. The dataframe names are pre-formated by the scraper so only pass the path to the directory into which the output data is to be written. If the directory does not already exist the ``cazy_webscraper`` will create the output dataframe.
* ``-s``, ``--subfamilies``: (True or False) Enable retrieval of subfamilies from CAZy. Default: false. If subfamilies are specified in the configuration file ensure ``-s`` is enabled.
* ``-v``, ``--verbose``: (True or False) Change the logger level from Warning to Info, resulting in logging of the scrapers progress. Default: false.


Configuration via a YAML file
------------------------------

The configuration file is for specifying specific CAZy classes and families to be scraped.


Specifying specific classes to scrape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under the **classes** heading list, any classes to be scrapped. For classes listed under 'classes', all proteins catalogued under that class will be retrieved, **unless** specific families have been listed under the respective classes heading in the configuration file. Then scraping only the specific classes takes precident and the entire class is not scraped. _If you believe this should be changed please raise an issue. It is invisioned that very few users would want to simultanious scrape an entire class and also scrape only specific families from that same class._

A ``cazy_dictionary.json`` has been created and packaged within the ``cazy_webscraper`` (the specific location is ./scraper/file_io/cazy_dictionary.json, where '.' is the directory where the webscraper is installed). 
This allows users to use a variety of synonoms for the CAZy classes, for example both "GH" and "Glycoside-Hydrolases" are accepted as synonoms for "Glycoside Hydrolases (GHs)". 
Additionally, the retrieval of CAZy classes from the configuration file is **not** case sensitive, therefore, both "gh" and "GH" are excepted.
If you having issues with the scraper retrieving the list of CAZy classes that are written under 'classes' in the configuration file, please check the dictionary first to see the full list of accepted synonoms. If you are comfortable modifying json files then feel free to add your own synonoms to the dictionary.

Each class must be listed on a separate line, indented by 4 spaces, and the class name encapsulated with single or double quotation marks. For example:

.. code-block:: yaml

    classes:
        - "GH"
        - "pl"


Specifying specific families to scrape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under the each of the class names listed in the configuration file list the names of specific **families** to be scraped from that class. You do not have to list the class of the families to be scraped under 'classes' as well, this is handled by the webscraper.

Write the true name of the family not only it's number, for example **GH1** is excepted by **1** is not. 
Name families using the standard CAZy nomenclature, such as "GT2" and NOT "GlycosylTransferases_2". 
Additionally, use the standard CAZy notation for subfamilies (**GH3_1**).

**Note:**
If any subfamilies are listed within the configuration file, the retrieval of subfamilies **must** 
be enabled at the command line uisng ``--subfamilies``.

Each family must be listed on a separate line and the name surrounded by double or single quotation marks. For example:

.. code-block:: yaml

    Glycoside Hydrolases (GHs):
        - "GH1"
        - "GH2"


Configuration when scraping subfamilies
---------------------------------------

If any subfamilies are listed within the configuration file, the retrieval of subfamilies **must** 
be enabled at the command line uisng ``--subfamilies``.

If the parent family, e.g GH3, is listed in the configuration file and `--subfamilies` is enabled, all proteins catalogued under GH3 and its subfamilies will be retrieved. This is to
save time having to write out all the subfamilies for a given CAZy family. The scraper will remove any duplicate proteins automatically.


An example configuration file
-----------------------------

A blank configuration file is packaged within `cazy_webscraper`, within the `scraper` directory, called `scraper_config.yaml`. 
This configuration file contains comments to assit filling in the file correctly. 
A new configuration file with any given name can be created and used. However, it **must** be a Yaml file and it **must** use the same headings/tags as used in the configuration file `scraper_config.yaml`.
Please find more information on writing lists in Yaml files [here](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html).

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


**Note:** indentations consist of 4 spaces.
