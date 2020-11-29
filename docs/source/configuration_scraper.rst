===========================================
Configuring the CAZy webscraper
===========================================

The CAZy webscraper is configurable to restrict the scraping of CAZy to specific classes and/or 
families, if scraping the entire CAZy database is not desired. Additionally, the retrieval of 
subfamilies, GenBank FASTA files and PDB protein structures can be enabled if desired.

The webscraper is configured via the **command line** and via a **YAML configuration file**.


Configuration via the command line
-----------------------------------

There are no required/positional arguments for the webscraper, therefore the scraper can be enabled 
by simply calling the scraper at the command line in the terminal: ``cazy_webscaper``

If no optional arguments are provided the default behaviour of the scraper will be performed. 
The default behaviour is to:

* Scrape the entire CAZy databases
* Not to split/separate the data, producing a single dataframe
* Write the resulting dataframe to standard out (STDOUT)
* Not to retrieve subfamilies (members of subfamilies will be retrieved but their parent family be listed)
* Not to retrieve FASTA files from GenBank
* Not to retrieve protein sequences from PDB


Options configurable at the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following behaviours of the ``cazy_webscraper`` can be configured at the command-line in the terminal:  

* Limit the scraping of CAZy to specific classes and/or families
* Split the scraped data by class and/or family
* Force writing out the data to an existing directory
* Write out a log file of the program being invoked
* Not delete content already present in the output directory
* Enable retrieving subfamilies
* Enable verbose logging during the operation of the webscraper


Command line options
^^^^^^^^^^^^^^^^^^^^

The below table lists the commands to be included with calling the webscraper in order to configure 
the webscraper's operation. For example:
``cazy_webscraper -s -v``

If an optional command is not given when calling the webscraper the relative default behaviour will 
be performed

The short option is the shortened version of the command, and the long option is the long hand 
method of writing the command, only one of the other needs be called not both.

==============  ======================  =============================================================================================================================================================================================================================================================================================  =============================================================================================================================================================================
 Short option    Long option             Behaviour and invoking/configuring the option                                                                                                                                                                                                                                                  Default behaviour
==============  ======================  =============================================================================================================================================================================================================================================================================================  =============================================================================================================================================================================
 ``-c``          ``--config``            Configure which classes and/or families that are to be scraped. Add the option to the command followed by the path to configuration YAML file.                                                                                                                                                 Write the generated dataframes to standard out.
 ``-d``          ``--data_split``        Split/separate the CAZymes into separate dataframes. Spliting the data by class or family. Call the option to the command then specify 'None' for not spliting the data, ' family' for separing the data by family, or 'class' for separating the data by CAZy family.                         Not to split the data. Produce a single dataframe containing all scraped CAZymes.
 ``-f``          ``--force``             Forcing writing out to directory that already exists, if specificed as the output directory. Simply add the option to the command.                                                                                                                                                             Force is false, thus if the output directory already exists the output will not be written to it.
 ``-g``          ``--genbank``           Enable retrieval of FASTA files from GenBank containing the protein sequence of the respective CAZyme. Add the option to the command then provide your email address (this is requirement of NCBI/GenBank).                                                                                    GenBank is false, FASTA files are not retrieved from GenBank.
                 ``--genbank_output``    Specify the output directory for FASTA files retrieved from GenBank to be written tothis applies for the dataframe, FASTA and protein structure output directories. Add the option to the command followed by the path to the desired directory to which the FASTA files are to be written.    Write the FASTA files to standard out.
 ``-h``          ``--help``              Displays all command-line options for the webscraper in the terminal, including defining their default behaviour and required additional information.                                                                                                                                          Not to display the help information
 ``-l``          ``--log``               Enable writing out a log file, logging the operation of the webscraper. Add the option to command followed by desired path of the resulting log file. This the path to the file not the directory to which the log file is to be written.                                                      Not to write out a log file of the webscrapers operation.
 ``-n``          ``--nodelete``          Do not delete content in the already existing output directory, this applies for the dataframe, FASTA and protein structure output directories. Simply add this option to the command.                                                                                                         Nodelete is false, delete the content already presented in the already existent output directories.
 ``-o``          ``--output``            Specify the directory to write the resulting dataframe of CAZymes to. This dataframe does not already have to exist, if it does not exist the webscraper will produce the directory. Add the option to the command followed by the path to the desired output directory.                       Write the output dataframes to standard out.
 ``-p``          ``--pdb``               Enable retrieval of protein structures from PDB. Call the option then specify the format of the resulting structure files. The available file formats are: mmCif, pdb, xml, mmtf, and bundle.                                                                                                  Not retrieve protein structures from PDB.
                 ``--pdb_output``        Specify the directory to which the protein structures are to be written. Add the option to the command followed by the path to desired output directory.                                                                                                                                       Write the structure files to the directory specified by ``--output``. If ``--output`` is standard out then the structure files are written to the current working directory
 ``-s``          ``--subfamilies``       Enable retrieval of subfamilies. If not enabled then the parent CAZy family will be listed for the relevant CAZymes. Simply add the option to the command.                                                                                                                                     Do not retrieve subfamilies from CAZy.
 ``-v``          ``--verbose``           Enable verbose logging of the webscraper. This provides more detailed logging of the progress of the webscrapers operation. Simply add the option to the command.                                                                                                                              Do not perform verbose logging. Only log if a warning or error is raised.
==============  ======================  =============================================================================================================================================================================================================================================================================================  =============================================================================================================================================================================


Example for configuring the webscraper
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below are some example commands for invoking the ``cazy_webscraper`` to help demonstrate how to configure the webscraper at the command line.

1. Writing the output to the directory 'my_output' and enabling retrieval of subfamilies: ``cazy_webscraper -o my_output -s``

2. Retrieving GenBank FASTA sequences and writing all output to standard out, not retrieve subfamilies, and verbose logging: ``cazy_webscraper -g example_email@domain.com -v``

3. Writing the output to an existing directory but not deleting the content already present in the directory: ``cazy_webscraper --output docs/my_output -f -n``

4. Retrieve protein structures, in the pdb format: ``cazy_webscraper -o my_output -p pdb --pdb_output my_output/cazyme_structures``


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
