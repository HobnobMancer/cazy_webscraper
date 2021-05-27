================================================================
Tutorials on installing, using and configuring cazy_webscraper
================================================================

``cazy_webscraper`` can be configured to retrieve user specified data sets from CAZy. The configuration 
applies to the retrieval of protein sequences from GenBank and protein structure files from PDB.

``cazy_webscraper`` can be configured via the **command line** and/or via a **YAML configuration file**.

This page runs through examples of how to combine the various 'filters' that can be applied, to fully customised 
the scraping of CAZy. These tutorials are designed for those with less experience using command-line tools.


Install cazy_webscraper
--------------------------

Checking the requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use ``cazy_webscraper`` you need to be using a POISx operating system (OS), Mac OS or a Linux emulator. This is because 
you need access to a Unix shell (a form of command-line) that writes in a language called Bash.

**Linux**
POISx includes operating systems such as Unix and it derivaties, Linux, and its derivatives such as Ubuntu. On Linux 
OS the default Unix shell is Bash.

**Mac**
If you are working on a Mac then you already have a Unix shell program installed. If you are running a 
Mac computer running Mac OS Mojave or earlier releases, the default Unix shell is using Bash. To access 
the Unix shell (or terminal) try one of the following:

* In Finder, select the Go menu, then select Utilities. Locate Terminal in the Utilities folder and open it
* Use the Mac ‘Spotlight’ computer search function. Search for: Terminal and press Return

To check if your machine is set up to use something other than Bash, type ``echo $SHELL`` in your terminal window.

For a Mac computer running macOS Catalina or later releases, or if you computer is set up to use something other 
than Bash, you can run Bash by opening a terminal, typing the command ``bash`` and hitting return.

**Windows**
If you use Windows you can use a Linux emulator, which allows your Windows computer to behave like a windows computer. Several 
emulators exist. Here we'll run through two different emulators you could use, and how to install them. These are not the only two 
and feel free to use which ever Linux emulator you wish!

**Git**
One Linux emulator you can use if included in *Git for Windows*. The following installation instructions are adapted from 
Software Carpentry `lessons <https://carpentries.github.io/workshop-template/#shell>`_. You can find a video tutorial with a step-by-step guide `here <https://youtu.be/339AEqk9c-8>`_.

Installation instrucitions:
* Download the Git for Windows install from `here <https://gitforwindows.org/>`_.
* Run the installer
* Click the "Next" button four times (only two times if you already have Git installed). *You do not need to change anything the Information, location, components, and start menu screens.
* From the dropdown menu select "Use the Nano editor by default" (NOTE: you will need to scroll up to find it) and click "Next"
* On the page that says "Adjusting the name of the initial branch in new repositories", ensure that "Let Git decide" is selected.
* Ensure that "Git from the command line and also from 3rd-party software" is selected and click on "Next"
* Ensure that "Use the native Windows Secure Channel Library" is selected and click on "Next".
* Ensure that "Checkout Windows-style, commit Unix-style line endings" is selected and click on "Next".
* Ensure that "Use Windows' default console window" is selected and click on "Next".
* Ensure that "Default (fast-forward or merge) is selected and click "Next"
* Ensure that "Git Credential Manager Core" is selected and click on "Next".
* Ensure that "Enable file system caching" is selected and click on "Next".
* Click on "Install".
* Click on "Finish" or "Next".
* If your "HOME" environment variable is not set (or you don't know what this is):
* Open command prompt (Open the Windows Start Menu, then type 'cmd', and press Enter). Wait for the command promt window to open, and type the following line into the command prompt window exactly as shown:
``setx HOME "%USERPROFILE%"``. Press Enter, you should see ``SUCCESS: Specified value was saved``. Quit command prompt by typing `exit` then pressing Enter

**Ubuntu**
Another Linux emulator that you can run on Windows is provided by Ubuntu. Ubuntu is a version of Linux with a graphical 
user interface (GUI), and they also provide a Windows emulator version that is fully supported, free and available via the `Microsoft app store <https://www.microsoft.com/en-gb/p/ubuntu-2004-lts/9n6svws3rx71#activetab=pivot:overviewtab>`_.

To install the Ubuntu emulator, navigate the `Microsoft store page <https://www.microsoft.com/en-gb/p/ubuntu-2004-lts/9n6svws3rx71#activetab=pivot:overviewtab>`_, and click install. 
This will handel all installation for you.

To start the terminal, open the Windows Start Menu and type 'Ubuntu', this will find the Ubuntu program. Click on Ubuntu and this will open a Ubuntu terminal. 

Ubuntu sets the Home as its own set of directories. To navigate to any hardrive in your system, open the Ubuntu terminal then type:  
``cd /mnt/<letter_of_harddrive (in lower case)>``, for example to access the C drive you would use ``cd /mnt/c``. The ``/mnt`` prefix 
accesses the Windows 'mounting' system, which is how it accesses harddrives within the computer.

**The shell terminal**
For more information on using the Unix shell checkout the content and lessons hosted at `SoftWare Carpentry <https://swcarpentry.github.io/shell-novice/01-intro/index.html>`_, which 
will walk through an introduction to the Unix shell and how the Unix shell can be used to optimise computational work.

Installing cazy_webscraper
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Preparation for installing**

The easiest way to install ``cazy_webscraper`` is to use `pip`. We can check if `pip` installed by opening a Bash terminal and typing:

.. code-block:: bash

   pip help

If `pip` is installed, text will appear in the terminal telling you how to use `pip`. If `pip` is not install an error message will be displayed, stating that the computer could not find a command called `pip`. If this happens install pip using:  

.. code-block:: bash
   
   python get-pip.py

To update the version of `pip` install on your computer use:  


.. code-block:: bash
   
   python -m pip install --upgrade pip
 
 
**Quick and easy installation**

The easiest way to install ``cazy_webscraper` is via Pypi distribtion, but using the command:

.. code-block:: bash

   pip3 install cazy-webscraper

Now ``cazy_webscraper`` is fully installed and you can skip to the part of the tutorial that explains how to use it!

.. note::

   ``cazy_webscraper`` is spelt with a DASH (-) when installing using `pip` **not** an underscore.

**Preparation for installing using the files stored in the GitHub repository**

Instead of installing ``cazy_webscraper`` via Pypi, you can install ``cazy_webscraper`` directly from the GitHub repository. Open up your Bash terminal. Then we need to navigate to directory where you want to install ``cazy_webscraper``. To do this we will use the ``cd`` command.

Just like how the windows explorer points at a single directory at any time, and shows you the content of the directory, the terminal acts the same way. 
To check at what directory your terminal is pointed at (or looking at, type the command ``pwd`` and press enter. The terminal will then 
return the path of the directory at which it is currently looking at. For example, if the terminal is pointed at a directory called 'my_data' within another directory called 'Documents', located on the C drive, 
``pwd`` will return ``c/Documents/my_dir/``.

We can change directory using the 'change directory' command (``cd``). Continuing on from the example above, 
if we wanted to move from the directory 'my_dir' into the directory 'cazyme_research' located within it, and then into 
the directory 'cazy_dir' within that, we would use the command ``cd cazyme_research/cazy_dir``.

The change directory (``cd``) command is called then provided a path to the directory that we wish to 
have the terminal pointed at. The ``cd`` command starts at the directory the terminal is currently looking at, then 
follows the path we provide it. This is why to move from 'my_dir' > 'cazy_research' > 'cazy_dir' we can type 
``cd cazy_research/cazy_dir``, becuase the terminal will looking within the 'my_dir' directory for the 'cazy_research' directory.

Using the ``cd`` command navigate to the directory you wish to install ``cazy_webscraper``. 
**If the directory where you wish to install ``cazy_webscraper`` does not exist we can use the terminal to make it**. 
To use the terminal, first use the ``cd`` comamand to navigate to the parent directory of where you wish to house the 
directory that you will install ``cazy_webscraper``. Then call the 'make directory command' ``mkdir`` followed by the name 
you wish to give the directory. For example, once we have navigated to the 'cazy_dir', we can make the directory 
'cazyme_databases' by using ``mkdir cazyme_database``. We can then navigate into the 'cazyme_database' directory we have justed made 
by typing ``cd cazy_database`` into the terminal and hitting Return.

**Installing ``cazy_webscraper``**

First we clone the GitHub repository, by using the code:

.. code-block:: bash

   git clone https://github.com/HobnobMancer/cazy_webscraper 

This creates a new directory into the directory that the terminal is currently pointed at, called 
'cazy_webscraper'. The command also downloads all files in the GitHub repository, and writes them into 
the new 'cazy_webscraper' directory.

We then need to move into the 'cazy_webscraper' directory:

.. code-block:: bash

   cd cazy_webscraper

We then use the Python package manage ``pip`` to install ``cazy_webscraper``.

.. code-block:: bash

   pip3 install -e .

Do not forget the **-e** from this command, otherwise ``cazy_webscraper`` will not be installed correctly 
and you will run into constant issues when trying to use ``cazy_webscraper``.

**If you ever invoke ``cazy_webscraper`` and want to cancle the command, simple press the ``Ctrl`` and ``c`` keys at the same time.**


Configuration via the command line
-----------------------------------

There are no required/positional arguments for the webscraper, therefore the scraper can be enabled 
by simply calling the scraper at the command line in the terminal: 

.. code-block:: bash
  python3 <path_to_cazy_webscraper.py_file>

The `cazy_webscraper.py` file is located within the directory `scraper`. Therefore, if the terminal 
is already pointing at the `scraper` directory, the command to invoke ``cazy_webscraper`` is:

.. code-block:: bash
  python3 cazy_webscraper.py

If you installed ``cazy_webscraper`` using the instructions above, the terminal will be pointing at the 
'cazy_webscraper' directory. We can move to the 'scraper' directory using the command ``cd scraper``.

When NO optional arguments are provided the default behaviour of the scraper will be performed. 
The default behaviour is to:

* Scrape the entire CAZy databases
* Write the resulting database to standard out (STDOUT)
* Not to retrieve subfamilies (members of subfamilies will be retrieved but only their parent family will be listed)


Options configurable at the command line
------------------------------------------

The following behaviours of the ``cazy_webscraper`` can be configured at the command-line in the terminal:  

* Limit the scraping of CAZy to specific CAZy classes, CAZy families, kingdoms, genuera, species, strains and/or EC numbers.
* Force writing out the database to a a new or existing directory
* Write out a log file of the packages operation
* Not delete content already present in the output directory
* Enable retrieving subfamilies
* Enable verbose logging during the operation of the webscraper


How to use the command-line options
--------------------------------------

The command-line options listed above can be used in any combination to customise the scraping of CAZy. The options that apply a 'filter' 
to restrict which CAZymes are scraped from CAZy are applied in combination. For example, if the ``--families`` option and ``--ec`` option are called then 
only CAZymes from the specified families **and** annotated with the listed EC numbers will be retrieved.

We will now walk through some examples of how to use ``cazy_webscraper``. All example code provided in this section will presume that the terminal is 
pointed at the `scraper` directory, which contains the `cazy_webscraper.py` file.


Configuring were the output is saved
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can name the directory that the database created by ``cazy_webscraper`` is written to by calling the ``--output`` flag. 
We add the flag to the command that invokes ``cazy_webscraper``. For example, to write the output to the directory 'cazyme_database' we can use:

.. code-block:: bash

   python3 cazy_webscraper.py --output cazyme_database

OR we can use the short hand version of the ``--output`` flag, ``-o``:

.. code-block:: bash

   python3 cazy_webscraper.py -o cazyme_database

The output directory does not have to exist when ``cazy_webscraper`` is invoked. ``cazy_webscraper`` can make 
a new directorty within the directory the terminal is currently pointing at. For example, if we are in the directory 
'scraper', and the directory 'cazyme_database' did not exist, then ``cazy_webscraper`` would make the directory 
'cazyme_database'.

The ``--output`` flag can take an infinetly long path. For example, we could use:

.. code-block:: bash

   python3 cazy_webscraper.py -o data/cazyme_research/cazyme_database

In the Bash terminal paths are **relative**, meaning that the terminal starts in the directory it is currently 
looking at and follows the path from there. The installation section of this tutorial covers this when 
discussing how to change directory.

.. warning::
   When requesting ``cazy_webscraper`` make an output directory, the parent of the directory we wish to make 
   **must already exist**. For examlple, if we asked ``cazy_webscraper`` to write its output to the directory 
   'data/cazyme_research/cazyme_database' and the directory 'cazyme_database' did not exist, *if* the directory 
   'cazyme_research' did exist ``cazy_webscraper`` would build the directory 'cazyme_database' within 'cazyme_research'. 
   However, if 'cazyme_research' *and* 'cazyme_database' did not exist, then ``cazy_webscraper`` would raise an error saying 
   the path 'data/cazyme_research' does not exist.


**Writing the output to an existing database**
If you want to write the output CAZyme database to a directory that already exists, you will need to add the force (``--force`` *or* ``-f``) flag 
anywhere to the ``cazy_webscraper`` command. For example:

.. code-block:: bash

   python3 cazy_webscraper.py -o data/cazyme_research/cazyme_database -f

By default ``cazy_webscraper`` will delete or content in an already existing output directory. Therefore, in the above example, 
if the directory ``cazyme_database`` already existed, ``cazy_webscraper`` would delete all content in the directory then proceed. 

You may wish to retain the data already in that directory. To do this add the 'no delete' (``--nodelete`` *or* ``-n``) flag anywhere 
to the ``cazy_webscraper`` command. For example:

.. code-block:: bash

   python3 cazy_webscraper.py -o data/cazyme_research/cazyme_database -f -n

The order you invoke *any* of the optional flags does not matter, for example the following three examples perform the 
exact same operation as the code given above:

.. code-block:: bash

   python3 cazy_webscraper.py --force -o data/cazyme_research/cazyme_database -f

.. code-block:: bash

   python3 cazy_webscraper.py -n -o data/cazyme_research/cazyme_database -f

.. code-block:: bash

   python3 cazy_webscraper.py --nodelete --force --output data/cazyme_research/cazyme_database

The above examples also highlight that it does not matter if you use the long or short versions of each of the flags.



Configuring CAZy classes and families
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Scraping specific classes**

If instead of scraping all of CAZy, you want to only scrape CAZymes from specific CAZy classes then add the 
``--classes`` flag followed by the classes you want to scrape. If you want to list multiple families, separate the families 
with a single comma. When you specify a CAZy class to scrape, *all* CAZy families within that class will be scraped.

For example, if you want to scrape all CAZymes from Glycoside Hydrolase and Carbohydrate Esterases then use the command:

.. code-block:: bash

   python3 cazy_webscraper.py --classes Glycoside Hydrolases,Carbohydrate Esterases

``cazy_webscraper`` excepts multiple synonyms for each CAZy class:

* **Glycoside Hydrolases (GHs):** Glycoside-Hydrolases, Glycoside-Hydrolases, Glycoside_Hydrolases, GlycosideHydrolases, GLYCOSIDE-HYDROLASES, GLYCOSIDE-HYDROLASES, GLYCOSIDE_HYDROLASES, GLYCOSIDEHYDROLASES, glycoside-hydrolases, glycoside-hydrolases, glycoside_hydrolases, glycosidehydrolases, GH, gh
* **GlycosylTransferases (GTs):** Glycosyl-Transferases, GlycosylTransferases, Glycosyl_Transferases, Glycosyl Transferases, GLYCOSYL-TRANSFERASES, GLYCOSYLTRANSFERASES, GLYCOSYL_TRANSFERASES, GLYCOSYL TRANSFERASES, glycosyl-transferases, glycosyltransferases, glycosyl_transferases, glycosyl transferases, GT, gt
* **Polysaccharide Lyases (PLs):** Polysaccharide Lyases, Polysaccharide-Lyases, Polysaccharide_Lyases, PolysaccharideLyases, POLYSACCHARIDE LYASES, POLYSACCHARIDE-LYASES, POLYSACCHARIDE_LYASES, POLYSACCHARIDELYASES, polysaccharide lyases, polysaccharide-lyases, polysaccharide_lyases, polysaccharidelyases, PL, pl
* **Carbohydrate Esterases (CEs):** Carbohydrate Esterases, Carbohydrate-Esterases, Carbohydrate_Esterases, CarbohydrateEsterases, CARBOHYDRATE ESTERASES, CARBOHYDRATE-ESTERASES, CARBOHYDRATE_ESTERASES, CARBOHYDRATEESTERASES, carbohydrate esterases, carbohydrate-esterases, carbohydrate_esterases, carbohydrateesterases, CE, ce
* **Auxiliary Activities (AAs):** Auxiliary Activities, Auxiliary-Activities, Auxiliary_Activities, AuxiliaryActivities, AUXILIARY ACTIVITIES, AUXILIARY-ACTIVITIES, AUXILIARY_ACTIVITIES, AUXILIARYACTIVITIES, auxiliary activities, auxiliary-activities, auxiliary_activities, auxiliaryactivities, AA, aa
* **Carbohydrate-Binding Modules (CBMs):** Carbohydrate-Binding-Modules, Carbohydrate_Binding_Modules, Carbohydrate_Binding Modules, CarbohydrateBindingModules, CARBOHYDRATE-BINDING-MODULES, CARBOHYDRATE_BINDING_MODULES, CARBOHYDRATE_BINDING MODULES, CARBOHYDRATEBINDINGMODULES, carbohydrate-binding-modules, carbohydrate_binding_modules, carbohydrate_binding modules, carbohydratebindingmodules, CBMs, CBM, cbms, cbm

.. note::
   These synonyms are stored in a JSON found at ``scraper/utilities/parse_configuration/cazy_dictionary.json``. 
   Storing these synonyms allows you to modify this file if you wish to add your own synonoms for each CAZy class.


**Scraping specific families**


To specify specific CAZy families to scrape, add the ``--families`` flag followed by the families you want 
to scrape. If you want to scrape multiple families, add the ``--families`` flag *once* followed by a list of *all* 
the CAZy families you want to scrape, separated by a single comma.

For example, if you want to scrape all CAZymes from GH2, PL5, CE1, CE2 and AA10 use:

.. code-block:: bash

   python3 cazy_webscraper.py --families GH2,PL5,CE1,CE2,AA10

Make sure to use the accepted CAZy nomenclature; 'GH2' is accepted but 'gh2' is not.

.. note::
   When ``--families`` is invoked any CAZy classes that do **not** include an of the CAZy families specified will 
   **not** be scraped. Therefore, using the example above, CAZymes from the families GH2, PL5, CE1, CE2 and AA10 
   **will** be retrieved; however, CAZymes from any other families from those classes **will not** be retrieved, and CAZymes 
   from the Carbohydrate Binding Modules (CBM) and GlycoslyTransferases classes will **not** be retrieved.


**Scraping specific classes AND families**


If you want to specify specific CAZy classes *and* families to scrape then add *both* the ``--classess`` *and* ``-families`` 
flags, because you can combine, mix-and-match, any combination of optional flags when invoking ``cazy_webscraper``.

For example, if we wanted to scrape all CAZymes from GH1, PL9 and *all* of CE we would use the command:

.. code-block:: bash

   python3 cazy_webscraper.py --families GH1,PL9 --classes CE

It does **not** matter what order you add the optional flags to your command. Therefore, if we wanted to 
scrape all CAZymes from PL1, PL2, PL3 and *all* of GH and CE we both:

.. code-block:: bash

   python3 cazy_webscraper.py --families PL1,PL2,PL3 --classes GH,CE

**AND**

.. code-block:: bash

   python3 cazy_webscraper.py --classes GH,CE --families PL1,PL2,PL3

are accepted.

.. note::
   In the example ``python3 cazy_webscraper.py --classes GH,CE --families PL1,PL2,PL3`` all CAZymes from PL1, 
   PL2 and PL3 would be retrieved, but no CAZymes from the other PL families, in addition all CAZymes from all GH and CE 
   families would be retrieved, but no CAZymes from AA, GT or CBM families would be retrieved.


Applying taxonomic and EC number filters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Specifying kingdoms**


You may only be interest in CAZymes that are derived from species from a specific taxonomic kingdom. 
CAZy classifies source organisms under one of 5 kingdoms:

* Archaea
* Bacteria
* Eukaryota
* Viruses
* Unclassified

To restrict the scraping of CAZy to retrieve CAZymes only derived from species from specific taxonomic kingdoms 
then add the ``--kingdoms`` flag to the ``cazy_webscraper`` command followed by the kingdoms to limit the retrieval 
of CAZymes to. To list multiple kingdoms you need only add the ``--kingdoms`` flag *once*, then list all the kingdoms 
you want to restrict the restrival of CAZymes to, separated by a single comma.

For example, if you want to retrieve CAZymes only from bacterial and eukaryotic species then use the command 

.. code-block:: bash

   python3 cazy_webscraper.py --kingdoms bacteria,eukaryota


.. warning::
   The kingomds must be spelt the same way CAZy spells them, for example use 'eukaryot**a**' instead of 'eukaryot**e**'. The kingdoms 
   are **not** case sensitive, therefore, both ``bacteria`` *and* ``Bacteria`` are accepted. You can also list the kingdoms in 
   *any* order. Thus, both ``bacteria,eukaryota`` *and* ``eukaryota,bacteria`` are accepted.


**Genera**

You can customise the scraping of CAZy to retrieve only CAZymes from *all* species from specific 
genera. To do this add the ``--genera`` flag to the ``cazy_webscraper`` command followed by all 
the genera you want to retrieve CAZymes from. CAZymes from any genera that you do not list will 
**not** be retrieved. To list multiple genera, you need to only add the ``--genera`` flag once followed 
by a list of all genera, with each genera separated with a single comma and *no* spaces.

For example, if we wanted to retrieve all CAZymes from *all* Aspergillus, Trichoderma and Streptomyces species 
we would use the command:

.. code-block:: bash

   python3 cazy_webscraper.py --genera Aspergillus,Trichoderma,Streptomyces


.. note::
   The order that the genera are listed does **not** matter. 


.. warning::
   Make sure to use the expect practise for writing genera names, each genus starts with a **captial** letter and 
   all other letters are lower case.

   Aspergillus is **correct**

   asepergillus is **incorrect**

   ASPERGILLUS is **incorrect**


**Species**


You can specify to retrieve CAZymes only derived from specific species. To do this add the ``--species`` 
flag to the ``cazy_webscraper`` command, followed by a list of all species you wish to retrist the retrieval of 
CAZymes to. Separate each species with a single comma. Also for each species use the full scientific name for the species.

For example, if we wanted to retrieve all CAZymes from Aspergillus niger and Aspergillus fumigatus we would use the command:  

.. code-block:: bash

   python3 cazy_webscraper.py --species Aspergillus niger,Asepergillus fumigatus


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

**Strains**

You may only be interested in specific strains of a species. Therefore, ``cazy_webscraper`` allows you to 
restrict the retrieval of CAZymes to only those derived from specific strains of species. To do this 
add the ``--strains`` flag to the ``cazy_webscraper`` command, followed by a list of all the strains 
of interest. Separate each strain with a single command and no spaces.

For example, if we wanted to retrieve all CAZymes from Aspergillus niger ATCC 1015 and Aspergillus uvarum CBS 121591  we would use the command:

.. code-block:: bash

   python3 cazy_webscraper.py --strains Aspergillus niger ATCC 1015,Aspergillus uvarum CBS 121591

.. note::
   The order that the strains are listed does **not** matter, and separate multiple species names with a single comma 
   with **no** spaces.

.. note::
   Sometimes in CAZy only the species name is given and no specific strain identifer. To retrieve CAZymes from these 
   species then you can list the species name and it will only retrieve CAZymes that are listed with the exact species 
   and with no strain identifers. For example, listing 'Aspergillus niger' will only retrieve CAZymes with their source 
   organism specifically listed as 'Aspergillus niger' and will not retrieve CAZymes from ''.

.. warning::
   If you use the ``--species``, ``--genera`` and ``--strains`` flags in any combination and a source organism matches 
   multiple of the taxonomy critera, the CAZymes derived from that species will only be retrieved **once**. For example, 
   using the command ``python3 cazy_webscraper --genera Aspergillus --species Aspergillus niger --strains Aspergillus niger ATCC 1015`` 
   will retrieve all CAZymes from *all* Aspergillus species *once*. The higher taxonomy levels take president, and the command 
   will not retrieve all CAZymes from all Aspergillus species once AND all CAZymes from Aspergillus niger strains as well, and then 
   retrieve another copy of all CAZymes from Aspergillus niger ATCC 1015.


**Combining taxonomic filters**

You can combine any combination of ``cazy_webscraper`` optional flags, including combining the taxonomic filters. For example,
you may wish to retrieve all CAZyme derived from all viral and Aspergillus species, Layia carnosa, Layia chrysanthemoides, Trichoderma reesei QM6a and 
Trichoderma reesei QM9414, we would combine the respective flags for a single ``cazy_webscraper`` command. The command 
we would use would be:

.. code-block:: bash

   python3 cazy_webscraper.py --kingdoms viruses --genera Aspergillus --species Layia carnosa,Layia chrysanthemoides --strains Trichoderma reesei QM6a,Trichoderma reesei QM9414

.. note::
   This is a single command written on a single line. When typing the command into the terminal do not fit enter until you have finished the command. 
   Visually the command may spread over multiple lines but it is a *single* command.

.. warning::
   If you use the ``--species``, ``--genera`` and ``--strains`` flags in any combination and a source organism matches 
   multiple of the taxonomy critera, the CAZymes derived from that species will only be retrieved **once**. For example, 
   using the command ``python3 cazy_webscraper --genera Aspergillus --species Aspergillus niger --strains Aspergillus niger ATCC 1015`` 
   will retrieve all CAZymes from *all* Aspergillus species *once*. The higher taxonomy levels take president, and the command 
   will not retrieve all CAZymes from all Aspergillus species once AND all CAZymes from Aspergillus niger strains as well, and then 
   retrieve another copy of all CAZymes from Aspergillus niger ATCC 1015.


**EC numbers**

If you are interested in CAZymes with specific activities you can limit the retrieval of CAZymes from CAZy to only those 
annotated with *at least one* EC number from a set of EC numbers you specify. To specify a set of EC numbers 
add the ``--ec`` flag to the ``cazy_webscraper`` command, followed by a list of EC numbers. Separate each EC number with a single 
comma and *no* spaces. Do **not** forget to include the 'EC' prefix from your EC numbers. 

.. note::
   Use the international accepted '-' (dash) to indicate missing identifiers (numbers) in the EC number.
   EC1.2.3.- is accepted but EC1.2.3. and EC1.2.3.* are not.

To limit the scraping of CAZy to only retrieve CAZymes that are annotated with *at least one* of the EC numbers 
EC4.2.2.-, EC1.3.2.- and EC5.4.-.-, use the command:

.. code-block:: bash

   python3 cazy_webscraper.py --ec "EC4.2.2.-,EC1.3.2.-,EC5.4.-.-"

.. warning::
   Some terminals may misinterpret ``EC1.2.-.-`` as trying to invoke the options ``.``, therefore, it is 
   recommend practise to encase the entire EC number list in single or double quotation marks if any of the EC numbers 
   include missing identifiers. ``"EC4.2.2.-,EC1.3.2.-,EC5.4.-.-"`` or ``'EC4.2.2.-,EC1.3.2.-,EC5.4.-.-'`` are recommended, 
   ``EC4.2.2.-,EC1.3.2.-,EC5.4.-.-`` is not recommended, and ``"EC4.2.2.-,EC1.3.2.-,EC5.4.-.-'`` (mismatching double and single 
   quotation marks) will raise errors.


**Taxonomy and EC numbers**

You can use any combination of ``cazy_webscraper`` optional flags to fully customise the scraping of CAZy. 
For example, you may which to retrieve all CAZymes annotated with the EC number EC4.2.2.- which are only from bacterial 
species. To do that you would add the ``--kingdoms`` and ``--ec`` flags:

.. code-block:: bash

   python3 cazy_webscraper.py --ec "EC4.2.2.-" --kingdoms bacteria

The order you add the optional flags **does not** matter, and you can specify multiple EC numbers, kingdoms, strains etc.


Combining Taxonomy, EC numbers, CAZy classes and CAZy families filters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The optional flags for ``cazy_webscraper`` can be used in any combination and any order. For example, 
you can combine the EC number, taxonomy, CAZy class and CAZy family configurations. Below are some examples:

**Example 1**  
To retrieve all CAZymes from all CBM families, GH1, GH2 and PL9, and that are derived from any Aspergillus species:

.. code-block:: bash

   python3 cazy_webscraper.py --classes CBM --families GH1,GH2,PL9 --genera Aspergillus

**Example 2**  
To retrieve all CAZymes from GH1, and GH2, if they are annotated with EC1.2.-.-, and are derived from any bacterial species:

.. code-block:: bash

   python3 cazy_webscraper.py --families GH1,GH2 --ec "EC1.2.-.-" --kingdoms bacteria 

.. warning::
   Some terminals may misinterpret ``EC1.2.-.-`` as trying to invoke the options ``.``, therefore, it is 
   recommend practise to encase the entire EC number list in single or double quotation marks if any of the EC numbers 
   include missing identifiers. ``"EC4.2.2.-,EC1.3.2.-,EC5.4.-.-"`` or ``'EC4.2.2.-,EC1.3.2.-,EC5.4.-.-'`` are recommended, 
   ``EC4.2.2.-,EC1.3.2.-,EC5.4.-.-`` is not recommended, and ``"EC4.2.2.-,EC1.3.2.-,EC5.4.-.-'`` (mismatching double and single 
   quotation marks) will raise errors.

**Example 3**  
To retrieve CAZymes from all viral species, and all Aspergillus niger strains which are catalogued within GH3_1 and GH3_2

.. code-block:: bash

   python3 cazy_webscraper.py --families GH3_1,GH3_2 --subfamilies --species Aspergillus niger --kingdoms Bacteria


Configuration file
------------------------------------

Whenever ``cazy_webscraper`` is invoked and adds data to a database, the configuration of ``cazy_webscraper`` 
(this is the kingdoms, genera, species, strains, EC numbers, CAZy classes and CAZy family filters which were applied) 
and the data and time the scrape was initiated is logged in the database. However, for optimal reproduction of 
how ``cazy_webscraper`` was used in your research, you can create shareable documentation that others can use to 
invoke ``cazy_webscraper`` and apply the exact sample filters as yourself. This is achieved by creating a configuration file 
rather than configuring the performance of ``cazy_webscraper`` at the command line.


Creating a configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An example and template configuration file is included in ``cazy_webscraper``, it can be found at ``scraper/scraper_config.yaml``. 
This is a YAML file; if you are new to YAML files please find more detailed information on YAML files [here](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html).

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
* kingoms
* ECs


**Specifying CAZy classes to scrape**

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

**Specifying CAZy families to scrape**

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

**Example configuration file**

Below is an example of the content you may wish to put in a configuration file.

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
   ECs:
      - EC4.2.2.-
      - EC5.4.-.-


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


Using a configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have created a configuration file (we recommend modifying the template one provided with ``cazy_webscraper`` 
you then need to invoke ``cazy_webscraper`` and tell it you are using a configuration file. To do this we add the 
``--config`` flag to the ``cazy_webscraper`` command, followed by the path to the configuration file.

.. note::
   You can use the long form of the configuration file flag (``--config``) *or* the short hand (``-c``).

The path we pass to ``cazy_webscraper`` is a *relative* path. This means ``cazy_webscraper`` will start in the directory 
the terminal is currently pointed out, and follow the path from there. For example, if we used the command:

.. code-block:: bash

   python3 cazy_webscraper.py -c scraper/scraper_config.yaml

Then the computer will look for a directory called ``scraper`` in the directory the terminal is looking at, then within the 
``scraper`` directory it will look for a yaml file called ``scraper_config.yaml``.

.. note::
   To check which directory ``cazy_webscraper`` is pointed at type ``pwd`` into the terminal and hit enter. This is the 
   'Present Working Directory' command, which will print the path to the directory the terminal is presently looking at.

.. warning::
   Your path must point directly to the YAML file. Don't forget the '.yaml' file extension!


Using a configuration and the command-line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you so wished, you can use a configuration file *and* the command line to configure ``cazy_webscraper``. If you do this 
``cazy_webscraper`` will **not** retrieve duplicates of the data. If a CAZyme matches at least one of the configuration data then 
one copy of the CAZyme record will be added to the SQL database, and only one copy will be added to the database no matter how many of the 
configuration data the CAZyme meets.

To use a configuration file and a the command-line to configure ``cazy_webscraper``, use the configuration file 
``--config`` flag followed by the path to the configuration file and any of the additional optional flags you wish to use.

.. note::
   The order you invoke the optional flags **does not** matter.


Additional operations to fine tune how ``cazy_webscraper`` operates
-----------------------------------------------------------------------


**Retrieving CAZy family and CAZy subfamily annotations**

The default behaviour of ``cazy_webscraper`` retrieves only the CAZy family annotations of CAZymes, and 
does **not** catalogue the child CAZy subfamily annotations as well. If you want to retrieve the CAZy subfamily 
annotations then add the ``--subfamilies`` flag anywhere to the ``cazy_webscraper`` command. For example:

.. code-block:: bash

   python3 cazy_webscraper.py --subfamilies


**Writing out a log file**

If you want to have a log file of all terminal output produced by ``cazy_webscraper`` then add the log 
``--log`` flag (or the shorthand version ``-l``) anywhere to the ``cazy_webscraper`` command, followed by a 
path to write the log file to. This path is a *relative* path and must include target a log file specifically. 
For example:

.. code-block:: bash

   python3 cazy_webscraper.py --subfamilies --genera Aspergillus --log log_dir/cazy_webscraper_log.log

.. warning::
   The log file does not already have to exist for ``cazy_webscraper`` to write to it; however, all 
   directories included in the path must already exist.


**Verbose logging**

For more detailed logging (logging more detail and not only when warnings and errors are raised by 
``cazy_webscraper``), add the verbose logging flag (``--verbose`` or ``-v``) anywhere to the ``cazy_webscraper`` 
command. You need only add the verbose flag and nothing else, for example:

.. code-block:: bash

   python3 cazy_webscraper.py --subfamilies --genera Aspergillus -v

The verbose flag can be used in combination with the log flag to write all terminal output to a log file.


**Changing connection timeout**

Sometimes the connection to the CAZy server times out. By default if a connection is attempted to made to CAZy 
and no response is recieved within 45 seconds, then ``cazy_webscraper`` interprets this as the connection 
timing out, waits 10 seconds and retries the connection.  You can change how long the computer waits for a 
response from the CAZy server before classifying the connection as timed out by adding the timeout flag to the 
``cazy_webscraper`` command, followed by the number of seconds you want the computer to wait for a response from CAZy 
before classifying the connection as timing out.

For example, to set the connection timeout limit to 30 seconds use the command:

.. code-block:: bash

   python3 cazy_webscraper.py --timeout 30

The timeout flag can be used in combination with other flags, for example:

.. code-block:: bash

   python3 cazy_webscraper.py --subfamilies --genera Aspergillus -v --timeout 30

You can use the long version ``--timeout`` or short version ``-t`` of the timeout flag.

.. code-block:: bash

   python3 cazy_webscraper.py --subfamilies --genera Aspergillus -v -t 60


**Add the scraped data to an existing CAZyme database**

You may wish to scrape CAZy in multiple stages, maybe your internet dropped out while scraping CAZy 
and you don't want to start again, or maybe you scraped CAZy but forget missed out a species of interest. No matter 
the reason ``cazy_webscraper`` allows you to add more CAZyme data to an existing database previously created by 
``cazy_webscraper``.

To do this add the database (``--database`` or ``-d``) flag to the ``cazy_webscraper`` command, followed by the path 
to the SQL database you want to add your scraped CAZy data to.

.. note::
   Don't forget the .db file extension at the end of the path!

All the paths we pass to ``cazy_webscraper`` are a *relative* path. This means ``cazy_webscraper`` will start in the directory 
the terminal is currently pointed out, and follow the path from there. For example, if we used the command:

.. code-block:: bash

   python3 cazy_webscraper.py -d my_cazyme_databases/my_cazyme_database.db

Then the computer will look for a directory called ``my_cazyme_databases`` in the directory the terminal is looking at, then within the 
``my_cazyme_databases`` directory the computer will look for the file ``my_cazyme_database.db``.



Configuration when scraping subfamilies
---------------------------------------

The default behaviour of ``cazy_webscraper`` retrieves only the CAZy family annotations of CAZymes, and 
does **not** catalogue the child CAZy subfamily annotations as well. If you want to retrieve the CAZy subfamily 
annotations then add the ``--subfamilies`` flag anywhere to the ``cazy_webscraper`` command. For example:

.. code-block:: bash

   python3 cazy_webscraper.py --subfamilies

This will retrieve both the parent CAZy family annotations and the child CAZy subfamily annotations for all applicable CAZymes. 
If a CAZyme is not part of a subfamily only its CAZy family annotations will be catagloued.

If any subfamilies are listed within the configuration file, the retrieval of subfamilies **must** 
be enabled at the command line uisng ``--subfamilies``.

If the parent family, e.g GH3, is listed in the configuration file and `--subfamilies` is enabled, 
all proteins catalogued under GH3 and its subfamilies will be retrieved. This is to save time 
having to write out all the subfamilies for a given CAZy family. The scraper will remove any 
duplicate proteins automatically.
