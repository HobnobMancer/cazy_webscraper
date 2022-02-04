================================================================
Tutorials on configuring ``cazy_webscraper`` to scrape CAZy
================================================================

``cazy_webscraper`` can be configured to retrieve user specified data sets from CAZy. The same configuration 
applies to the retrieval of protein data from UniProt, GenBank and PDB.

``cazy_webscraper`` can be configured via the **command line** and/or via a **YAML configuration file**.

This page runs through examples of how to combine the various 'filters' that can be applied, to fully customised 
the scraping of CAZy. These tutorials are designed for those with less experience using command-line tools.

.. NOTE::
  If you installed ``cazy_webscraper`` using ``bioconda`` or ``pip`` to invoke ``cazy_webscraper`` call the application using the command ``cazy_webscraper`` - this is the method used in this tutorial.  

  If you installed ``cazy_webscraper`` from source then you will need to invoke ``cazy_webscraper`` using the command ``python3 <path to cazy_scraper.py>``. For example, if you were located in the root of the repository, you would use: ``python3 cazy_webscraper/cazy_scraper.py``.


----------------------------------
Configuration via the command line
----------------------------------

``cazy_webscraper`` has only one required argument, the user email address. Therefore, 
 the scraper can be enabled to scrape all of CAZy using the following command:

.. code-block:: bash

  cazy_webscraper myemail@domain.com

When no optional arguments are provided the default behaviour of the scraper will be performed. 
The default behaviour is to:

* Scrape the entirety of CAZy databases
* Write the resulting database to the current working directory
* Not retrieve subfamilies (members of subfamilies will be retrieved but only their parent family will be listed)

.. NOTE::
   **For those new to using command-line tools:**  
   Arguments are additional pieces of information we add onto the end of the command. They are used to configure the specific behaviour 
   performed by yjr computer when we tell it to perform a specific command. In the examples above the command is ``cazy_webscraper myemail@domain.com``, 
   where we have told to computer to run the Python program ``cazy_webscraper`` and submit the user email 
   address to NCBI for the retrieval of source orgnaism data. No additional information was provided, the computer 
   will invoke ``cazy_webscraper`` using its default behaviour. If you do not want the default behaviour of ``cazy_webscraper`` then we need to 
   pass additionally arguments to the computer when telling it to run ``cazy_webscraper``, which we cover in the section below.


-----------------------------------------
Options configurable at the command line 
-----------------------------------------

The following behaviours of the ``cazy_webscraper`` can be configured at the command-line in the terminal:  

* Limit the scraping of CAZy to specific CAZy classes, CAZy families, kingdoms, genuera, species, and/or strains
* Force writing out the database to a a new or existing directory
* Write out a log file of the packages operation
* Not delete content already present in the output directory
* Enable retrieving subfamilies
* Enable verbose logging during the operation of the webscraper

`Here <https://cazy-webscraper.readthedocs.io/en/latest/usage.html>`_ you can find a full list of the command-line flags and options.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
How to use the command-line options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The command-line options listed above can be used in any combination to customise the scraping of CAZy. The options that apply a 'filter' 
to restrict which CAZymes are scraped from CAZy are applied in combination. For example, if the ``--families`` option and ``--genera`` option are called then 
only CAZymes from the specified families **and** belonging to source organisms within the defined genera will be retrieved.

We will now walk through some examples of how to use ``cazy_webscraper``. All example code presumes ``cazy_webscraper`` was installled using 
``Bioconda`` or ``pip`` and therefore, be simply called using the command ``cazy_webscraper``.

.. NOTE::
   **For those new to using command-line tools: flags**
   Command-line flags are used to tell the computer specifically which option(s) to change. Flags **always** come after the command.
   
   The abbreivated 
   version of a flag is given the prefixed of a single dash, followed by a single letter. For example, ``-s``, ``-o`` and ``-l`` are all examples of short 
   hand flags.
   
   The long version of a flag is prefixed by two dashes, followed by complete words. For example, ``--output`` is the long version of the ``-o``. 

   The flags used by a program are defined within the program. This means the flag ``-o`` may represent different options for different programs. Always make 
   sure to check the documentation to see what flags are provided with the program, and what they do!
   
   You can use the command-line to list all flags for a program/tool by typing in the command to invoke that tool, followed by the flag ``--help`` or ``-h``. For example: 
   ``cazy_webscraper --help``.


-------------------------------------
Configuring were the output is saved
-------------------------------------

^^^^^^^^^^^^^^^^^^^^^^^
Creating a new database
^^^^^^^^^^^^^^^^^^^^^^^

Instead of writing out the database to the current working directory using the default database name 
(``cazy_webscraper_<date>_<time>.db``), we can name the database and directory that the database 
created by ``cazy_webscraper`` is written to by calling the ``--output`` flag. 

We add the flag to the command that invokes ``cazy_webscraper``. For example, to write the output to the directory 'cazyme_database' with the file 
name 'cazyme_database.db' we can use:

.. code-block:: bash

   cazy_webscraper --output cazyme_database/cazyme_database.db

OR we can use the short hand version of the ``--output`` flag, ``-o``:

.. code-block:: bash

   cazy_webscraper -o cazyme_database/cazyme_database.db

.. NOTE::
   The final element of the path provided after the ``--output`` / ``-o`` flag is the name of database compiled by 
   ``cazy_webscraper``.

The output directory does not have to exist when ``cazy_webscraper`` is invoked. ``cazy_webscraper`` can make 
the output directory, including all necessary parent directories. 

The ``--output`` flag can take an infinetly long path. For example, we could use:

.. code-block:: bash

   cazy_webscraper -o data/cazyme_research/cazyme_database/cazyme_database.db

If the directories ``cazymes_research`` and ``cazyme_database`` did not exist, then ``cazy_webscraper`` will build 
these for you.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Overwriting an existing database or directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to write the output CAZyme database to a directory and/or file that already exists, you will need to add the force (``--force`` *or* ``-f``) flag 
anywhere to the ``cazy_webscraper`` command. For example:

.. code-block:: bash

   cazy_webscraper -o data/cazyme_research/cazyme_database/cazyme_database.db -f

By default ``cazy_webscraper`` will delete all content in an already existing output directory. Therefore, in the above example, 
if the directory ``cazyme_database`` already existed, ``cazy_webscraper`` would delete all content in the directory then proceed. 

You may wish to retain the data already in that directory. To do this add the 'no delete' (``--nodelete`` *or* ``-n``) flag anywhere 
to the ``cazy_webscraper`` command. For example:

.. code-block:: bash

   cazy_webscraper -o data/cazyme_research/cazyme_database/cazyme_database.db -f -n

The order you invoke *any* of the optional flags does not matter, for example the following three examples perform the 
exact same operation as the code given above:

.. code-block:: bash

   cazy_webscraper --force -o data/cazyme_research/cazyme_database/cazyme_database.db -f

.. code-block:: bash

   cazy_webscraper -n -o data/cazyme_research/cazyme_database/cazyme_database.db -f

.. code-block:: bash

   cazy_webscraper --nodelete --force --output data/cazyme_research/cazyme_database/cazyme_database.db

The above examples also highlight that it does not matter if you use the long or short versions of each of the flags.

.. NOTE::
   If you elect to write the database to a file in the current working directory, you do not need to worry 
   about ``cazy_webscraper`` deleting content in the current working directory. This only applies if you chose to
   write the database to a directory over than the current working directory.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Add the scraped data to an existing CAZyme database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may wish to scrape CAZy in multiple stages; maybe your internet dropped out while scraping CAZy 
and you don't want to start again, or maybe you scraped CAZy but missed out a species of interest. No matter 
the reason ``cazy_webscraper`` allows you to add more CAZyme data to an existing database previously created by 
``cazy_webscraper``.

To do this add the database (``--database`` or ``-d``) flag to the ``cazy_webscraper`` command, followed by the path 
to the CAZyme database you want to add your scraped CAZy data to. For example, to add data to an existing 
database in ``cazy/cazyme_db.db`` use the command:

.. code-block:: bash

   cazy_webscraper -- database cazy/cazyme_db.db

.. note::
   Don't forget the .db file extension at the end of the path!

All the paths we pass to ``cazy_webscraper`` are a *relative* path. This means ``cazy_webscraper`` will start in the directory 
the terminal is currently pointed out, and follow the path from there. For example, if we used the command:

.. code-block:: bash

   cazy_webscraper -d my_cazyme_databases/my_cazyme_database.db

Then the computer will look for a directory called ``my_cazyme_databases`` in the directory the terminal is looking at, then within the 
``my_cazyme_databases`` directory the computer will look for the file ``my_cazyme_database.db``.


----------------------------------------------
Specifying CAZy classes and families to scrape
----------------------------------------------

^^^^^^^^^^^^^^^^^^^^^^^^^
Scraping specific classes
^^^^^^^^^^^^^^^^^^^^^^^^^

If instead of scraping all of CAZy, you want to only scrape CAZymes from specific CAZy classes then add the 
``--classes`` flag followed by the classes you want to scrape. If you want to list multiple classes, separate the classes 
with a single comma. When you specify a CAZy class to scrape, *all* CAZy families within that class will be scraped.

For example, if you want to scrape all CAZymes from Glycoside Hydrolase and Carbohydrate Esterases then use the command:

.. code-block:: bash

   cazy_webscraper --classes Glycoside Hydrolases,Carbohydrate Esterases

``cazy_webscraper`` excepts multiple synonyms for each CAZy class:

* **Glycoside Hydrolases (GHs):** Glycoside-Hydrolases, Glycoside-Hydrolases, Glycoside_Hydrolases, GlycosideHydrolases, GLYCOSIDE-HYDROLASES, GLYCOSIDE-HYDROLASES, GLYCOSIDE_HYDROLASES, GLYCOSIDEHYDROLASES, glycoside-hydrolases, glycoside-hydrolases, glycoside_hydrolases, glycosidehydrolases, GH, gh

* **GlycosylTransferases (GTs):** Glycosyl-Transferases, GlycosylTransferases, Glycosyl_Transferases, Glycosyl Transferases, GLYCOSYL-TRANSFERASES, GLYCOSYLTRANSFERASES, GLYCOSYL_TRANSFERASES, GLYCOSYL TRANSFERASES, glycosyl-transferases, glycosyltransferases, glycosyl_transferases, glycosyl transferases, GT, gt

* **Polysaccharide Lyases (PLs):** Polysaccharide Lyases, Polysaccharide-Lyases, Polysaccharide_Lyases, PolysaccharideLyases, POLYSACCHARIDE LYASES, POLYSACCHARIDE-LYASES, POLYSACCHARIDE_LYASES, POLYSACCHARIDELYASES, polysaccharide lyases, polysaccharide-lyases, polysaccharide_lyases, polysaccharidelyases, PL, pl

* **Carbohydrate Esterases (CEs):** Carbohydrate Esterases, Carbohydrate-Esterases, Carbohydrate_Esterases, CarbohydrateEsterases, CARBOHYDRATE ESTERASES, CARBOHYDRATE-ESTERASES, CARBOHYDRATE_ESTERASES, CARBOHYDRATEESTERASES, carbohydrate esterases, carbohydrate-esterases, carbohydrate_esterases, carbohydrateesterases, CE, ce

* **Auxiliary Activities (AAs):** Auxiliary Activities, Auxiliary-Activities, Auxiliary_Activities, AuxiliaryActivities, AUXILIARY ACTIVITIES, AUXILIARY-ACTIVITIES, AUXILIARY_ACTIVITIES, AUXILIARYACTIVITIES, auxiliary activities, auxiliary-activities, auxiliary_activities, auxiliaryactivities, AA, aa

* **Carbohydrate-Binding Modules (CBMs):** Carbohydrate-Binding-Modules, Carbohydrate_Binding_Modules, Carbohydrate_Binding Modules, CarbohydrateBindingModules, CARBOHYDRATE-BINDING-MODULES, CARBOHYDRATE_BINDING_MODULES, CARBOHYDRATE_BINDING MODULES, CARBOHYDRATEBINDINGMODULES, carbohydrate-binding-modules, carbohydrate_binding_modules, carbohydrate_binding modules, carbohydratebindingmodules, CBMs, CBM, cbms, cbm

.. TIP::
   These synonyms are stored in a JSON found at ``scraper/utilities/parse_configuration/cazy_dictionary.json``. 
   Storing these synonyms allows you to modify this file if you wish to add your own synonoms for each CAZy class.

If you have your own synonyms these can be used by using the ``--class_synonyms`` flag, followed by the path to your JSON file. This JSON file **must** have the same 
architecture as the JSON filed used by ``cazy_webscraper``.

^^^^^^^^^^^^^^^^^^^^^^^^^^
Scraping specific families
^^^^^^^^^^^^^^^^^^^^^^^^^^

To specify specific CAZy families to scrape, add the ``--families`` flag followed by the families you want 
to scrape. If you want to scrape multiple families, list all the families you wish to scrape, with each family 
separated with a single comma.

For example, if you want to scrape all CAZymes from GH2, PL5, CE1, CE2 and AA10 use:

.. code-block:: bash

   cazy_webscraper --families GH2,PL5,CE1,CE2,AA10

.. WARNING::
   Make sure to use the accepted CAZy nomenclature; 'GH2' is accepted but 'gh2' is not.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Scraping specific classes AND families
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to specify specific CAZy classes *and* families to scrape then add *both* the ``--classess`` *and* ``-families`` 
flags, because you can combine, mix-and-match, any combination of optional flags when invoking ``cazy_webscraper``.

For example, if we wanted to scrape all CAZymes from GH1, PL9 and *all* of CE we would use the command:

.. code-block:: bash

   cazy_webscraper --families GH1,PL9 --classes CE

It does **not** matter what order you add the optional flags to your command. Therefore, if we wanted to 
scrape all CAZymes from PL1, PL2, PL3 and *all* of GH and CE, both:

.. code-block:: bash

   cazy_webscraper --families PL1,PL2,PL3 --classes GH,CE

**AND**

.. code-block:: bash

   cazy_webscraper --classes GH,CE --families PL1,PL2,PL3

are accepted.

.. note::
   In the example ``cazy_webscraper --classes GH,CE --families PL1,PL2,PL3`` all CAZymes from PL1, 
   PL2 and PL3 would be retrieved, but no CAZymes from the other PL families, in addition all CAZymes from all GH and CE 
   families would be retrieved, but no CAZymes from AA, GT or CBM families would be retrieved.


------------------
Applying taxonomic
------------------

^^^^^^^^^^^^^^^^^^^
Specifying kingdoms
^^^^^^^^^^^^^^^^^^^

You may only be interest in CAZymes that are derived from species from a specific taxonomic kingdom. 
CAZy classifies source organisms under one of 5 kingdoms:

* Archaea
* Bacteria
* Eukaryota
* Viruses
* Unclassified

To restrict the scraping of CAZy to retrieve CAZymes only derived from species from specific taxonomic kingdoms 
add the ``--kingdoms`` flag to the ``cazy_webscraper`` command followed by the kingdoms to limit the retrieval 
of CAZymes to. To list multiple kingdoms you need only add the ``--kingdoms`` flag *once*, then list all the kingdoms 
you are interested in, separated by a single comma.

For example, if you want to retrieve CAZymes only from bacterial and eukaryotic species then use the command 

.. code-block:: bash

   cazy_webscraper --kingdoms bacteria,eukaryota


.. warning::
   The kingdoms must be spelt the same way CAZy spells them, for example use 'eukaryot**a**' instead of 'eukaryot**e**'.
   
.. NOTE:: 
   The kingdoms are **not** case sensitive, therefore, both ``bacteria`` *and* ``Bacteria`` are accepted. 

.. NOTE::
   You can list the kingdoms in *any* order. Thus, both ``bacteria,eukaryota`` *and* ``eukaryota,bacteria`` are accepted.


^^^^^^^^^^^^^^^^^^^^^^^^^^
Speciying Genera to scrape
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can customise the scraping of CAZy to retrieve only CAZymes from *all* species from specific 
genera. To do this add the ``--genera`` flag to the ``cazy_webscraper`` command followed by your
genera of interes.

To list multiple genera, you need to only add the ``--genera`` flag *once* followed 
by a list of all your genera, with each genera separated with a single comma and *no* spaces.

For example, if we wanted to retrieve all CAZymes from *all* Aspergillus, Trichoderma and Streptomyces species 
we would use the command:

.. code-block:: bash

   cazy_webscraper --genera Aspergillus,Trichoderma,Streptomyces


.. note::
   The order that the genera are listed does **not** matter. 


.. warning::
   Make sure to use the expect practise for writing genera names, each genus starts with a **captial** letter and 
   all other letters are lower case.

   Aspergillus is **correct**

   asepergillus is **incorrect**

   ASPERGILLUS is **incorrect**


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifying species of organisms to scrape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can specify to retrieve only CAZymes derived from specific species. To do this add the ``--species`` 
flag to the ``cazy_webscraper`` command, followed by a list of all species you wish to retrist the retrieval of 
CAZymes to. Separate each species with a single comma. Also for each species use the full scientific name for the species.

For example, if we wanted to retrieve all CAZymes from *Aspergillus niger* and *Aspergillus fumigatus* we would use the command:  

.. code-block:: bash

   cazy_webscraper --species Aspergillus niger,Asepergillus fumigatus


.. note::
   The order that the species are listed does **not** matter, and separate multiple species names with a single comma 
   with **no** spaces.

.. warning::
   Use the standard scientific name formating. Captialise the first letter of *genus* and write a lower 
   case letter for the first letter of the species.

   Aspergillus niger is **correct**

   asepergillus niger is **incorrect**

   ASPERGILLUS NIGER is **incorrect**


.. warning::
   When you specify a species ``cazy_webscraper`` will retrieval CAZymes from *all* strains of the species.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specify specific strains of species to scrape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may only be interested in specific strains of a species. Instead of scraping CAZymes for all strains 
of a given speices, add the ``--strains`` flag followed by the specific species strains you wish to restrict 
the retrieval of CAZymes to.

List the full scientific name followed by the strain name. To specify multiple strains, list all 
strains of interest and separate with a single comma with **no** space.

For example, if we wanted to retrieve all CAZymes from Aspergillus niger ATCC 1015 and Aspergillus uvarum CBS 121591  we would use the command:

.. code-block:: bash

   cazy_webscraper --strains Aspergillus niger ATCC 1015,Aspergillus uvarum CBS 121591

he order that the strains are listed does **not** matter.

.. NOTE::
   If you use the ``--species``, ``--genera`` and ``--strains`` flags in any combination and a source organism matches 
   multiple of the taxonomy critera, the CAZymes derived from that species will only be retrieved **once**.
   
   For example, using the command ``cazy_webscraper --genera Aspergillus --species Aspergillus niger --strains Aspergillus niger ATCC 1015`` 
   will retrieve all CAZymes from *all* Aspergillus species *once*.
   
The higher taxonomy levels take president, and the command will not retrieve all CAZymes from all Aspergillus species once AND all CAZymes from Aspergillus niger strains as well, and then retrieve another copy of all CAZymes from Aspergillus niger ATCC 1015.


^^^^^^^^^^^^^^^^^^^^^^^^^^^
Combining taxonomic filters
^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can combine any combination of ``cazy_webscraper`` optional flags, including combining the taxonomic filtersFor example,
you may wish to retrieve all CAZyme derived from all viral species, Aspergillus species, Layia carnosa, Layia chrysanthemoides, Trichoderma reesei QM6a and 
Trichoderma reesei QM9414. To do this we would combine the respective flags for a single ``cazy_webscraper`` command. The command 
we would use would be:

.. code-block:: bash

   cazy_webscraper --kingdoms viruses --genera Aspergillus --species Layia carnosa,Layia chrysanthemoides --strains Trichoderma reesei QM6a,Trichoderma reesei QM9414

.. note::
   This is a single command written on a single line. When typing the command into the terminal do not fit enter until you have finished the command. 

.. warning::
   If you use the ``--species``, ``--genera`` and ``--strains`` flags in any combination and a source organism matches 
   multiple of the taxonomy critera, the CAZymes derived from that species will only be retrieved **once**. For example, 
   using the command ``cazy_webscraper --genera Aspergillus --species Aspergillus niger --strains Aspergillus niger ATCC 1015`` 
   will retrieve all CAZymes from *all* Aspergillus species *once*.

When combining taxonomy filters, the higher taxonomy levels take president. For example, the :command:
   
.. code-block:: bash

   cazy_webscraper --genera Aspergillus --species Aspergillus niger --strains Aspergillus niger ATCC 1015

will not retrieve all CAZymes from all Aspergillus species once AND all CAZymes from Aspergillus niger strains as well. 
``cazy_webscraper`` will retrieval all CAZymes for all strains of *Aspergillus niger*.

-----------------------------------------
Enabling retrieving subfamily annotations
-----------------------------------------

By default ``cazy_webscraper`` only retrieves the CAZy family annotation for each protein, it does not 
retrieve the CAZy subfamily annotation. For example, a CAZyme within the CAZy subfamily GH3_1, will be 
stored in the local CAZyme database as only a GH3 CAZyme.

To retrieve the CAZy family **and** CAZy subfamily annotations, add the ``-subfamilies`` / ``-s`` flag, anywhere in the 
``cazy_webscraper`` command. For example:

.. code-block:: bash

   cazy_webscraper --families GH3 --subfamilies

This command will retrieve all CAZymes from GH3, and will retrieve the CAZy family **and** CAZy subfamily 
annotations. For example, a CAZyme in CAZy subfamily GH3_1 will be stored in the local database under the 
CAZy family GH3 and the CAZy subfamily GH3_1.

------------------------------------------------------
Combining CAZy class, CAZy family and taxonomy filters
------------------------------------------------------

You can use any combination of the CAZy class, CAZy family and taxonomy filters to fully customise the scrape of 
CAZy.

Below are some examples:

**Example 1**
To retrieve all CAZymes from all CBM families, GH1, GH2 and PL9, and that are derived from any Aspergillus species:

.. code-block:: bash

   cazy_webscraper --classes CBM --families GH1,GH2,PL9 --genera Aspergillus

**Example 2**  
To retrieve all CAZymes from GH1, and GH2 that are derived from any bacterial species:

.. code-block:: bash

   cazy_webscraper --families GH1,GH2 --kingdoms bacteria 

**Example 3**  
To retrieve CAZymes from all viral species, and all Aspergillus niger strains which are catalogued within GH3_1 and GH3_2

.. code-block:: bash

   cazy_webscraper --families GH3_1,GH3_2 --subfamilies --species Aspergillus niger --kingdoms Bacteria

------------------
Configuration file
------------------

Whenever ``cazy_webscraper`` is invoked and adds data to a database, the configuration of ``cazy_webscraper`` 
(this is the kingdoms, genera, species, strains, CAZy classes and CAZy family filters which were applied) 
and the data and time the scrape was initiated is logged in the database. However, for optimal reproduction of 
how ``cazy_webscraper`` was used in your research, you can create shareable documentation that others can use to 
reproduce your CAZy dateset. This is achieved by creating a configuration file 
rather than configuring the performance of ``cazy_webscraper`` at the command line.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Creating a configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``cazy_webscraper`` uses the YAML file type for its configuraiton file; 
if you are new to YAML files please find more detailed information on YAML files [here](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html).

A template and example configuration file for scrapping CAZy using ``cazy_webscraper`` can be found in 
the repo, in the ``configuration_files`` directory.

The configuration YAML **must** contain the following tags/headings (identical to how they are presented below):

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
* kingdoms

.. NOTE::
   The order of the tags/headings does not matter.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Scraping specific CAZy classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under the **classes** heading list any classes to be scrapped. For each CAZy class listed under 'classes', CAZymes 
will be retrieved for every CAZy family within the CAZy class.

Each class must be listed on a separate line, indented by 4 spaces, and the class name encapsulated 
with single or double quotation marks. For example:

.. code-block:: yaml

    classes:
        - "GH"
        - "PL"

The same CAZy class name synonyms used for the command line are accepted for the configuration file.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Scraping specific CAZy families
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under the each of the class names listed in the configuration file, list the names of specific 
**families** to be scraped from that class. The respective classes of the specificed families do 
**not** need to be added to the 'classes' list.

Write the true name of the family not only it's number, for example **GH1** is excepted by **1** is 
not.

Name families using the standard CAZy nomenclature, such as **"GT2"** and 
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
        - "GH3_1"


^^^^^^^^^^^^^^^^^^^^^^^^^^
Example configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^

Below is an example of the content you may wish to put in a configuration file. Using this file 
will retrieve all CAZymes in CAZy class AA, CAZy families GH1, GH3 and PL9 that are either derived from 
a bacterial or *Trichoderma* species.

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
   genera:
      - "Trichoderma"
   species:
   strains:
   kingdoms:
      - "Bacteria"


.. note::
    Indentations consist of 4 spaces.


You can add 'comments' to configuration file. Comments are section of text that are not read by ``cazy_webscraper`` and 
allow you to add notes to your configuration file. For example:


.. code-block:: yaml

   # This is a comment, text following a hashtag '#' on the same line is not read by cazy_webscraper
   # https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html 
   classes:  # classes from which all proteins will be retrieved
   Glycoside Hydrolases (GHs):  # include two spaces between the end of the code and the hashtag
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
   ECs:  # only CAZymes with at least one of these EC numbers will be scrapped


Example configuration files and tempate files can be found `hre <https://github.com/HobnobMancer/cazy_webscraper/tree/master/configuration_files>`_.


^^^^^^^^^^^^^^^^^^^^^^^^^^
Using a configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have created a configuration file (we recommend modifying the template one provided with ``cazy_webscraper`` 
you then need to invoke ``cazy_webscraper`` and tell it you are using a configuration file. To do this we add the 
``--config`` / ``-c`` flag to the ``cazy_webscraper`` command, followed by the path to the configuration file.

The path we pass to ``cazy_webscraper`` is a *relative* path. This means ``cazy_webscraper`` will start in the directory 
the terminal is currently pointed out, and follow the path from there. For example, if we used the command:

.. code-block:: bash

   cazy_webscraper -c scraper/scraper_config.yaml

Then the computer will look for a directory called ``scraper`` in the directory the terminal is looking at, then look within the 
``scraper`` directory for a yaml file called ``scraper_config.yaml``.

.. note::
   To check which directory ``cazy_webscraper`` is pointed at type ``pwd`` into the terminal and hit enter. This is the 
   'Present Working Directory' command, which will print the path to the directory the terminal is presently looking at.

.. warning::
   Your path must point directly to the YAML file. Don't forget the '.yaml' file extension!

------------------------------------------
Using a configuration and the command-line
------------------------------------------

You can configure ``cazy_webscraper`` using a combination of command line arguments and a configuration file. 

If a CAZyme matches at least one of the configuration data (whether if be from the terminal of the configuration file),  
one copy of the CAZyme record will be added to the SQL database, and only **one copy**, no matter how many of the 
configuration data the CAZyme meets.

To use a configuration file and a the command-line to configure ``cazy_webscraper``, use the configuration file 
``--config`` flag followed by the path to the configuration file and any of the additional optional flags you wish to use.

.. note::
   The order you invoke the optional flags **does not** matter.


-------------------------------------------------------------------
Additional operations to fine tune how ``cazy_webscraper`` operates
-------------------------------------------------------------------


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Scraping data from a previously downloaded CAZy txt file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CAZy provides access to data within its database via text files. ``cazy_webscraper`` downloads the CAZy 
text file containing all data within the CAZy database, providing a database dump. This file is then written to the cache directory 
(by default, called ``.cazy_webscraper_<date>_<time>``).

For consistency in the dataset, you may wish to perform multiple scrapes of CAZyme data from the same CAZy text file. 
This could be a CAZy text file you have downloaded from CAZy or a text file downloaded by ``cazy_webscrapper``.

To direct ``cazy_webscraper`` to retrieve CAZyme data from a previously downloaded CAZy text file, using the 
``--cazy_data`` flag, followed by the path to the text file. For example:

.. code-block:: bash
   
   cazy_webscraper --cazy_data cazy_db/cazy_data.txt

.. WARNING::
   ``--cazy_data`` must be pointed directly at the text file, **not** a zipped file containing the CAZy 
   data text file.


^^^^^^^^^^^^^^^^^^^^^^
Writing out a log file
^^^^^^^^^^^^^^^^^^^^^^

If you want to have a log file of all terminal output produced by ``cazy_webscraper`` then add the log 
``--log`` / ``-l`` anywhere to the ``cazy_webscraper`` command, followed by a 
path to write the log file to. This path is a *relative* path and must include target a log file specifically. 
For example:

.. code-block:: bash

   cazy_webscraper --subfamilies --genera Aspergillus --log log_dir/cazy_webscraper_log.log

.. warning::
   The log file does not already have to exist for ``cazy_webscraper`` to write to it; however, all 
   directories included in the path must already exist.

^^^^^^^^^^^^^^^
Verbose logging
^^^^^^^^^^^^^^^

For more detailed logging (which includes not only error and warning messages (the default) but also 
configuration setup, number of proteins retrieved etc.), add the verbose logging flag (``--verbose`` or ``-v``) anywhere to the ``cazy_webscraper`` 
command. For example:

.. code-block:: bash

   cazy_webscraper --subfamilies --genera Aspergillus -v

The verbose flag can be used in combination with the log flag to write all terminal output to a log file.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changing the connection timeout limit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes the connection to the CAZy server times out. By default if a connection is attempted to made to CAZy 
and no response is recieved within 45 seconds, then ``cazy_webscraper`` interprets this as the connection 
timing out. ``cazy_webscraper`` then waits 10 seconds and retries the connection.

You can change how long the computer waits for a 
response from the CAZy server before classifying the connection as timed out by adding the ``--timeout`` flag to the 
``cazy_webscraper`` command, followed by the number of seconds you want the computer to wait for a response from CAZy 
before classifying the connection as timing out.

For example, to set the connection timeout limit to 30 seconds use the command:

.. code-block:: bash

   cazy_webscraper --timeout 30

The timeout flag can be used in combination with other flags, for example:

.. code-block:: bash

   cazy_webscraper --subfamilies --genera Aspergillus -v --timeout 30
