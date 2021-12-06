==========================================================================
Tutorials on configuring ``cazy_webscraper`` to retrieve data from PDB
==========================================================================

``cazy_webscraper`` can be configured to retrieve user specified data sets from UniProt for specific date 
sets of CAZymes in a local CAZyme database. Many of the configuration options 
apply to the retrieval of protein data from CAZy, UniProt, GenBank and PDB.

.. NOTE::
   ``cazy_webscraper`` retrieves protein data from UniProt for CAZymes in a local CAZyme database. It 
   will not add new proteins to the database.

``cazy_webscraper`` can be configured via the **command line** and/or via a **YAML configuration file**.

This page runs through examples of how to combine the various 'filters' that can be applied, to fully customised 
the retrieval of data from UniProt. These tutorials are designed for those with less experience using command-line tools.

.. NOTE::
  If you installed ``cazy_webscraper`` using ``bioconda`` or ``pip`` to invoke ``cazy_webscraper`` to retrieve UniProt data call it using ``cw_get_uniprot_data`` - this is the method used in this tutorial.  
  If you installed ``cazy_webscraper`` from source then you will need to invoke ``cazy_webscraper`` from the root of the repo using the command ``python3 cazy_webscraper/expand/uniprot/get_uniprot_data.py``.

From this point on, we will be discussed the ``cw_get_uniprot_data``, which is the entry point for 
retrieving data from UniProt. We also presume you are comfortable configuring ``cazy_webscraper`` for the 
scraping of data from CAZy.

All data retrieved from UniProt by ``cw_get_uniprot_data`` is added to the local CAZyme database.

----------------------------------
Configuration via the command line
----------------------------------

``cw_get_uniprot_data`` has only one required argument: the path to the local CAZyme database created 
using ``cazy_webscraper``. Therefore, ``cw_get_uniprot_data`` can be enabled using a simple command structure:

.. code-block:: bash

  cazy_webscraper <path to the local CAZyme db>

For example, if our database was stored in ``cazy/cazyme.db``, we would used:

.. code-block:: bash
   
  cazy_webscraper cazy/cazyme.db

When no optional arguments are provided, the default behaviour is invoked. The default behaviour is to:

* Retrieve protein data for **all** CAZymes in the local CAZyme db
* Retrieve UniProt protein accessions and protein names

-----------------------------------------
Options configurable at the command line 
-----------------------------------------

The following behaviours of the ``cw_get_uniprot_data`` can be configured at the command-line in the terminal:  

* Limit the retrieve of UniProt protein data to CAZymes in the local databaes from specific CAZy classes, CAZy families, kingdoms, genuera, species, strains and/or EC numbers
* Enable retrieving EC number annotations
* Enable retrieving PDB accessions
* Enable retrieving protein amino acid sequences
* Enable updating protein sequences in the local CAZyme database if newer versions are retrieved from UniProt
* Enable verbose logging during the operation of the webscraper

`Here <https://cazy-webscraper.readthedocs.io/en/latest/configuration_scraper.html>`_ you can find a full list of the command-line flags and options.


---------------------------------------------------------------
Retrieving protein data for CAZy classes and families to scrape
---------------------------------------------------------------

The ``--classes`` and ``--families`` flags from scraping data from CAZy are applied in the extact same way 
for retrieve data from UniProt.

For instance, if instead of retrieving protein data for all CAZymes in your local CAZyme database, you want to 
retrieve protein data for CAZymes in specific CAZy classes then add the 
``--classes`` flag followed by the classes you want to retrieve protein data for.

.. TIP::
   To list multiple classes, separate the classes with a single comma. 

For example, if you want to retrieve protein data for all CAZymes from Glycoside Hydrolase and Carbohydrate Esterases then use the command:

.. code-block:: bash

   cw_get_uniprot_data cazy/cazyme.db --classes GH,CE

OR

.. code-block:: bash

   cw_get_uniprot_data cazy/cazyme.db --classes Glycoside Hydrolases,Carbohydrate Esterases

Retrieving protein data for proteins from specific specific CAZy families is achieved using the ``--families`` flag. For 
example, to retrieve protein data for all proteins in PL1, PL2 and PL3 in the local CAZyme database use the 
following command:

.. code-block:: bash

   cw_get_uniprot_data cazy/cazyme.db --families PL1,PL2,PL3

.. WARNING::
   ``cw_get_uniprot_data`` only accpets families written in the proper CAZy family syntax.
   GH1 is accepted.
   gh1 and GlycosideHydrolases1 are not accepted.

As with scraping data from CAZy, the ``--classes`` and ``--families`` flags can be combined. To retrieve 
protein data for all CAZymes in PL1, PL2, PL3 and *all* of GH and CE both:

.. code-block:: bash

   cw_get_uniprot_data cazy/cazyme.db --families PL1,PL2,PL3 --classes GH,CE

**AND**

.. code-block:: bash

   cw_get_uniprot_data cazy/cazyme.db --classes GH,CE --families PL1,PL2,PL3

are accepted.


------------------
Applying taxonomic
------------------

The ``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` flags can be used to refine the dataset 
of proteins to retrieve protein data by taxonomy. These flags are applied in the exact same way as they 
are used for the scraping of data from CAZy. Only proteins in the local CAZyme database and matching at least on of the provided taxonomy 
criteria will have protein data retrieved from UniProt and added to the local CAZyme datbase.

For example, if you want to retrieve protein data for all CAZymes in a local CAZyme database from bacterial and eukaryotic species then use the command 

.. code-block:: bash

   cw_get_uniprot_data cazy/cazyme.db --kingdoms bacteria,eukaryota

.. warning::
   The kingdoms must be spelt the same way CAZy spells them, for example use 'eukaryot**a**' instead of 'eukaryot**e**'.
   
.. NOTE:: 
   The kingdoms are **not** case sensitive, therefore, both ``bacteria`` *and* ``Bacteria`` are accepted. 

.. NOTE::
   You can list the kingdoms in *any* order. Thus, both ``bacteria,eukaryota`` *and* ``eukaryota,bacteria`` are accepted.

You can combine any combination of the optional flags, including combining the taxonomic filters. For example,
you may wish to retrieve protein data for all CAZymes in a local CAZyme database that are derived from all viral species, Aspergillus species, Layia carnosa, Layia chrysanthemoides, Trichoderma reesei QM6a and 
Trichoderma reesei QM9414. To do this we would combine the respective flags for a single ``cw_get_uniprot_data`` command. The command 
we would use would be:

.. code-block:: bash

   cw_get_uniprot_data cazy/cazyme.db --kingdoms viruses --genera Aspergillus --species Layia carnosa,Layia chrysanthemoides --strains Trichoderma reesei QM6a,Trichoderma reesei QM9414

.. note::
   The order that the flags are used and the order taxa  are listed does **not** matter, and separate multiple taxa names with a single comma 
   with **no** spaces.

.. warning::
   Use the standard scientific name formating. Captialise the first letter of *genus* and write a lower 
   case letter for the first letter of the species.

   Aspergillus niger is **correct**

   asepergillus niger is **incorrect**

   ASPERGILLUS NIGER is **incorrect**

.. warning::
   When you specify a species ``cw_get_uniprot_data`` will retrieval CAZymes from *all* strains of the species.


-------------------------
Applying EC number filter
-------------------------

The retrieval of protein data from UniProt can also be limited to proteins in a local CAZyme database that are
annotated with specific EC numbers.

Having previously retrieved EC number annotations and added them to the local CAZyme database, you  may 
wish to retrieve protein data for CAZymes annotated with specific EC numbers. To do this add the 
``--ec_filter`` flag to the command, follwed by a list of EC numbers.

.. code-block:: bash
   
   cw_get_uniprot_data cazy/cazyme.db --ec_filter "EC1.2.3.4,EC2.3.4.5"


.. NOTE::
    Provide complete EC numbers. 
    Both dashes ('-') and asterixes ('*') are accepted for missing digits in EC numbers.

    EC1.2.3.- and EC1.2.3.* are accepted.
    EC1.2.3. and EC 1.2.3 are **not** accepted.

.. NOTE::
   The 'EC' prefix is not necessary.
   EC1.2.3.4 and 1.2.3.4 are accepted.

.. WARNING::
    If using dashes to represent missing digits in EC numbers, it is recommended to bookend the entire 
    EC number list in single or double quotation marks. Some terminals may misinterpret EC1.2.-.- as trying to invoke the options '.'

.. NOTE::
    ``cazy_webscraper`` will retrieve the specified UniProt data for all proteins in the local CAZyme 
    database that are annotated with **at least one** of the given EC numbers. Therefore, if multiple 
    EC numbers are given this **does not mean** only CAZymes will all provided EC numbers will have data retrieved
    from UniProt for them.

``--ec_filter`` is based upon EC number annotations stored within the local CAZyme database. For 
example, if protein A is annotated with the EC1.2.3.4, but this annotation is not stored in the 
local CAZyme database, using ``--ec_filter EC1.2.3.4`` will **not** cause ``cazy_webscraper`` to retrieve
data for protein A. This is because ``cazy_webscraper`` does not know protein A is annotated with 
EC1.2.3.4, because this annotation is not within its database.

.. WARNING::
    If ``--ec_filter`` is used along side ``--ec``, ``cazy_webscraper`` will retrieve **all** EC number 
    annotations from UniProt for all proteins in the local CAZyme database that are associated with 
    at least one of the EC numbers provided via ``--ec_filter`` within the CAZyme database.


---------------------
Combining all filters
---------------------

The ``--classes``, ``--families``, ``--ec_filter``, ``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` flags can 
be used in any combination to define a specific subset of proteins in the local CAZyme database for whom
protein data from UniProt will be retrieved. These flags can be used with any combination of 
``--ec``, ``--pdb``, ``--sequence``, ``--update_seq`` to customise what data is retrieved from UniProt and 
added to the local CAZyme database.

Below we run through 3 example commands of combining these flags, and the resulting behaviour.

**Example 1:**
To retrieve PDB accessions for all CAZymes in GH, GT, CE1, CE5 and CE8, and which are derived from baceterial species we use the command:

.. code-block:: bash

   cw_get_uniprot_data cazy/cazyme.db --pdb --classes GH,CE --families CE1,CE5,CE8 --kingdoms bacteria


**Example 2:**
To retrieve EC numbers and PDB accessions for all CAZymes in GH and which are derived from *Aspegillus* and *Trichoderma* species we use the command:

.. code-block:: bash

   cw_get_uniprot_data cazy/cazyme.db --pdb --classes GH --genera Aspegillus,Trichoderma


**Example 3:**
To retrieve EC numbers and sequences for all CAZymes in GH,CE and CBM which are derived from baceterial species and are annotated with at least one of 
EC3.2.1.23, EC3.2.1.37 and EC3.2.1.85, we use the command:

.. code-block:: bash

   cw_get_uniprot_data cazy/cazyme.db --pdb --sequences --classes GH,CE,CBM --kingdoms bacteria --ec_filter "3.2.1.23,3.2.1.37,3.2.1.85"

