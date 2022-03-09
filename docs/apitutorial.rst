================================================
Tutorial on interrogating the data using the API
================================================

``cazy_webscraper`` includes an API that can be used to interrogate the data in the local CAZyme database and write out the retrieved data 
in ``JSON`` and/or ``CSV`` format.

By default ``cazy_webscraper`` only includes the GenBank accessions of proteins matching the provided 
criteria, but the inclusion of additional data (such as protein squences, UniProt accessions, EC numbers, etc) 
is fully customisable. 

For interrogating the data and retrieving proteins matching the provided criteria, many of the same configuration options 
apply to the retrieval of protein data from CAZy, UniProt, GenBank and PDB.

``cazy_webscraper`` can be configured via the **command line** and/or via a **YAML configuration file**.

This page runs through examples of how to combine the various 'filters' that can be applied, to fully customised 
the retrieval of protein data from the local CAZyme database. These tutorials are designed for those with less experience using command-line tools.

.. NOTE::
  If you installed ``cazy_webscraper`` using ``bioconda`` or ``pip`` to invoke ``cazy_webscraper`` to retrieve UniProt data call it using ``cw_query_database`` - this is the method used in this tutorial.  
  If you installed ``cazy_webscraper`` from source then you will need to invoke ``cazy_webscraper`` from the root of the repo using the command ``python3 cazy_webscraper/api/cw_query_database.py``.

From this point on, we will be discusseing the ``cw_query_database`` command, which is used by ``cazy_webscraper`` for 
retrieving protein data from the local CAZyme database. We also presume you are comfortable configuring ``cazy_webscraper`` for the 
scraping of data from CAZy.


----------------------------------
Configuration via the command line
----------------------------------

``cw_query_database`` requires two arguments:
* The path to the local CAZyme database created using ``cazy_webscraper``
* The file formats to write out the output

Therefore, ``cw_query_database`` can be enabled using a simple command structure:

.. code-block:: bash

  cazy_webscraper <path to the local CAZyme db> <file formats>


For example, if our database was stored in ``cazy/cazyme.db`` and we want to write the output to a csv file, we would used:

.. code-block:: bash
   
  cazy_webscraper cazy/cazyme.db csv

.. NOTE::
   Make sure ``cw_query_database`` is pointed directly at the database file.

When no optional arguments are provided, the default behaviour is invoked. The default behaviour is to: 
retrieve only the GenBank accessions of **all** CAZymes in the local CAZyme CAZyme db


---------------------
Accepted file formats
---------------------

``cw_query_database`` can write the output to a csv or json file.

These are provided as the second arguments to ``cw_query_database``. To write out both a csv and json 
file use both ``csv`` and ``json`` after the path to the local CAZyme database, separted with a single space.

.. code-block:: bash

    cw_query_database <path to local CAZyme db> csv json

.. NOTE::
    The order ``csv`` and ``json`` are written does not matter.

.. WARNING::
    Both ``csv`` and ``json`` are case sensitive.


-----------------------------------------
Options configurable at the command line 
-----------------------------------------

The following behaviours of the ``cw_query_database`` can be configured at the command-line in the terminal:  

* Limit the retrieval of protein data to CAZymes in the local databaes from specific CAZy classes, CAZy families, kingdoms, genuera, species, strains and/or EC numbers
* Including any combination of the following in the output, along side the GenBank accessions:
    * CAZy class
    * CAZy family
    * CAZy subfamily
    * Kingdom
    * Genus
    * Scientific name of the source organsim
    * Protein sequence retrieved from GenBank
    * UniProt accession
    * Protein name retrieved from UniProt
    * EC numbers
    * PDB accessions
    * Protein sequence retrieved from UniProt
* Write the output in JSON and/or CSV format
* Choose an output directory
* Force overwriting existing files
* Enable verbose logging during the operation of the webscraper

`Here <https://cazy-webscraper.readthedocs.io/en/latest/api.html>`_ you can find a full list of the command-line flags and options.

----------------------------
Choosing an output directory
----------------------------

By default, ``cw_query_database`` writes all output files to the current working directory.

To specify an alternative output directory, using the ``--output_dir``/``-o`` flag, followed by the path to the target output directory. ``cw_query_database`` will build all necessary parent and child output directories.

For example, to write the output to the directory ``my_cazy_data`` use the following command:

.. code-block:: bash
  cw_query_database cazy/cazyme.db json csv --output_dir my_cazy_data
  
If the output directory already exists, ``cw_query_database`` will raise an error warning the output directory already exists and close. This is to prevent accidently writing data to the wrong location.

To force ``cw_query_database`` to write the data to an existing output directory, add the ``--force``/``-f`` flag.

.. code-block:: bash
  cw_query_database cazy/cazyme.db json csv --output_dir my_cazy_data --force

By default ``cw_query_database`` will delete all content already present in the existing output directory. To retain the data in the existing output directory, add the ``--nodelete``/``-n`` flag.

.. code-block:: bash
  cw_query_database cazy/cazyme.db json csv --output_dir my_cazy_data -- force --nodelete
  
.. note::
  The ``--force`` and ``--nodelete`` flags are only applied when the ``--output_dir`` flag is used. ``cw_query_database`` will **not** delete content in the current working directory when writing to the current working directory when the ``--output_dir`` flag is **not** used.
  
------------------------
Overwrite existing files
------------------------

``cw_query_database`` automatically compiles the names of the output files. 

The file names of all output files are composed of the name of the local CAZyme database, followed by the names of the data retrieved from the local CAZyme database. For example, retrieving the following data from the local CAZyme database called ``cazy_database.db``:
* CAZy family annotation
* CAZy subfamily annotations
* EC numbers
* PDB accessions
Will produce the following file name: ``cazy_database_gbkAcc_fams_subfams_ec_pdb``.  

.. note::
  ``_gbkAcc`` is always included in the file name because GenBank accessions are always retrieved and written to the output by ``cw_query_database``.

Both the `json` and `csv` files are given the same name, the files only differ in their file extension.

An optional prefix can be applied to all output file names using the ``-p``/``--prefix`` flag, followed by the desired prefix. For example, using the same example as above, the prefix 'engineering_candidates_` can be applied to every output file by adding the following to command:

.. code-block:: bash
  --prefix engineering_candidates_

This will produce output files with the file name ``engineering_candidates_cazy_database_fams_subfams_ec_pdb``.

If files matching the file names compiled by ``cw_query_database`` already existing at the target output location, ``cw_query_database`` will raise a warning that output files already existing and terminate. This is to prevent accidently overwriting data files.

To overwrite existing datafiles add the ``--overwrite`` flag to the command. For example, the following command will retrieve all GenBank accessions stored in the local CAZyme database located at ``cazy/cazyme.db`` and write out the GenBank accessions to a file called ``all_gbk_acc_cazyme_gbkAcc.csv`` to ``my_cazy_data``, and will not delete content in the existing output directory and will overwrite the existing output file ``my_cazy_data/all_gbk_acc_cazyme_gbkAcc.csv``.

.. code-block:: bash
  cw_query_database cazy/cazyme.db csv \
  --output_dir my_cazy_data \
  --prefix all_gbk_accs_
  -- force \
  --nodelete \
  --overwrite
 
--------------------------------------------------------------------
Retrieving protein data for CAZy classes and families to scrape
--------------------------------------------------------------------

The ``--classes`` and ``--families`` flags from scraping data from CAZy are applied in the extact same way 
for retrieving protein data from the local CAZyme databases.

For instance, if instead of retrieving protein data for all CAZymes in your local CAZyme database, you want to 
retrieve protein data for CAZymes in specific CAZy classes then add the 
``--classes`` flag followed by the classes you want to retrieve protein data for.

.. TIP::
   To list multiple classes, separate the classes with a single comma. 

For example, if you want to retrieve protein data for all CAZymes from Glycoside Hydrolase and Carbohydrate Esterases, and write the data to a csv file, then use the command:

.. code-block:: bash

   cw_query_database cazy/cazyme.db csv --classes GH,CE

OR

.. code-block:: bash

   cw_query_database cazy/cazyme.db csv --classes Glycoside Hydrolases,Carbohydrate Esterases

Retrieving protein data for proteins from specific specific CAZy families is achieved using the ``--families`` flag. For 
example, to retrieve protein data for all proteins in PL1, PL2 and PL3 in the local CAZyme database, and write the 
data to csv and json files, use the following command:

.. code-block:: bash

   cw_query_database cazy/cazyme.db json csv --families PL1,PL2,PL3

.. WARNING::
   ``cw_query_database`` only accpets families written in the proper CAZy family syntax.
   GH1 is accepted.
   gh1 and GlycosideHydrolases1 are not accepted.

As with scraping data from CAZy, the ``--classes`` and ``--families`` flags can be combined. To retrieve 
protein data for all CAZymes in PL1, PL2, PL3 and *all* of GH and CE both, and write the data to a json file:

.. code-block:: bash

   cw_query_database cazy/cazyme.db json --families PL1,PL2,PL3 --classes GH,CE

**AND**

.. code-block:: bash

   cw_query_database cazy/cazyme.db json --classes GH,CE --families PL1,PL2,PL3

are accepted.


------------------
Applying taxonomic
------------------

The ``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` flags can be used to refine the dataset 
of proteins to retrieve protein data by taxonomy. These flags are applied in the exact same way as they 
are used for the scraping of data from CAZy. Only proteins in the local CAZyme database and matching at least on of the provided taxonomy 
criteria will have protein data retrieved from GenBank and added to the local CAZyme datbase.

For example, if you want to retrieve protein data for all CAZymes in a local CAZyme database from bacterial and eukaryotic species then use the command 

.. code-block:: bash

   cw_query_database cazy/cazyme.db csv --kingdoms bacteria,eukaryota

.. warning::
   The kingdoms must be spelt the same way CAZy spells them, for example use 'eukaryot**a**' instead of 'eukaryot**e**'.
   
.. NOTE:: 
   The kingdoms are **not** case sensitive, therefore, both ``bacteria`` *and* ``Bacteria`` are accepted. 

.. NOTE::
   You can list the kingdoms in *any* order. Thus, both ``bacteria,eukaryota`` *and* ``eukaryota,bacteria`` are accepted.

You can combine any combination of the optional flags, including combining the taxonomic filters. For example,
you may wish to retrieve protein data for all CAZymes in a local CAZyme database that are derived from all viral species, Aspergillus species, Layia carnosa, Layia chrysanthemoides, Trichoderma reesei QM6a and 
Trichoderma reesei QM9414. To do this we would combine the respective flags for a single ``cw_query_database`` command. The command 
we would use would be:

.. code-block:: bash

   cw_query_database cazy/cazyme.db csv --kingdoms viruses --genera Aspergillus --species Layia carnosa,Layia chrysanthemoides --strains Trichoderma reesei QM6a,Trichoderma reesei QM9414

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
   When you specify a species ``cw_query_database`` will retrieval CAZymes from *all* strains of the species.


-------------------------
Applying EC number filter
-------------------------

The retrieval of protein data from the local CAZyme database can also be limited to proteins in a local CAZyme database that are
annotated with specific EC numbers.

Having previously retrieved EC number annotations and added them to the local CAZyme database, you  may 
wish to retrieve protein data for CAZymes annotated with specific EC numbers. To do this add the 
``--ec_filter`` flag to the command, follwed by a list of EC numbers.

.. code-block:: bash
   
   cw_query_database cazy/cazyme.db csv --ec_filter "EC1.2.3.4,EC2.3.4.5"


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
protein data from GenBank will be retrieved.

Below we run through 3 example commands of combining these flags, writing the output to a csv file, and the resulting behaviour.

**Example 1:**
To retrieve protein data for all CAZymes in GH, GT, CE1, CE5 and CE8, and which are derived from baceterial species we use the command:

.. code-block:: bash

   cw_query_database cazy/cazyme.db csv --classes GH,CE --families CE1,CE5,CE8 --kingdoms bacteria


**Example 2:**
To protein data for all CAZymes in GH and which are derived from *Aspegillus* and *Trichoderma* species we use the command:

.. code-block:: bash

   cw_query_database cazy/cazyme.db csv -classes GH --genera Aspegillus,Trichoderma


**Example 3:**
To retrieve protein data for all CAZymes in GH,CE and CBM which are derived from baceterial species and are annotated with at least one of 
EC3.2.1.23, EC3.2.1.37 and EC3.2.1.85, we use the command:

.. code-block:: bash

   cw_query_database cazy/cazyme.db csv --ec --sequences --classes GH,CE,CBM --kingdoms bacteria --ec_filter "3.2.1.23,3.2.1.37,3.2.1.85"


----------------------
Customising the output
----------------------

By defauly ``cw_query_database`` only includes the GenBank accessions of the CAZymes matching the provided 
criteria in the final output. Any combination of the following can also be included in the output:
* CAZy class
* CAZy family
* CAZy subfamily
* Kingdom
* Genus
* Scientific name of the source organsim
* Protein sequence retrieved from GenBank
* UniProt accession
* Protein name retrieved from UniProt
* EC numbers
* PDB accessions
* Protein sequence retrieved from UniProt

To include additional data in the output use the ``--include`` flag followed by any combination (and any order) of the following options:
``--include`` - List additional data to include in the output. Multiple fields can be named, separating each with a single space (' '). The accepted fields are:
* 'class' - Include the CAZy class annotations
* 'family' - Include the CAZy family annotations
* 'subfamily' - Include the subfamily class annotations
* 'kingdom' - Include the taxonomic kingdom of the source organism
* 'genus' - Include the genus of the source organism
* 'organism' - Include the scientific name of the source organism
* 'uniprot_acc' - Include the UniProt accession
* 'uniprot_name' - Include the protein name retrieved from UniProt
* 'ec' - Include the EC number annotations
* 'pdb' - Include the PDB accessions
* 'genbank_seq' - Include the GenBank protein sequence
* 'uniprot_seq' - Include the Uniprot protein sequence

.. NOTE::
   The quotation marks around the terms do not need to be included.

.. NOTE::
   No matter what additional data is included in the output, the data will be presented in the same 
   order as presented above. For example, 'Kingdom' will always come before all fields listed below it. 
   Changing the order the fields are listed in the command will not change the order data is presented in the 
   output. For example, using ``--include kingdom genus organism`` and ``--include organism genus kingdom`` will 
   both result in the respective columns being placed in the order of: 'Kingdom', 'Genus', 'Organism'

To list multiple fields to include in the final output, separate each field with a singel space (' ').

**Example 1:**
To retrieve protein data for all CAZymes in GH, GT, CE1, CE5 and CE8, and which are derived from baceterial species, and include the CAZy family annotations and 
scientific names of the source organisms we use the command:

.. code-block:: bash

   cw_query_database cazy/cazyme.db csv --classes GH,CE --families CE1,CE5,CE8 --kingdoms bacteria --include family organism


**Example 2:**
To protein data for all CAZymes in GH and which are derived from *Aspegillus* and *Trichoderma* species, and include the CAZy class, EC number and PDB accessions 
in the output we use the command:

.. code-block:: bash

   cw_query_database cazy/cazyme.db csv --include class ec pdb --classes GH --genera Aspegillus,Trichoderma


**Example 3:**
To retrieve protein data for all CAZymes in GH,CE and CBM which are derived from baceterial species and are annotated with at least one of 
EC3.2.1.23, EC3.2.1.37 and EC3.2.1.85, and include the EC number annotations, CAZy family and CAZy subfamily annotations we use the command:

.. code-block:: bash

   cw_query_database cazy/cazyme.db csv \
      --ec --sequences \
      --classes GH,CE,CBM \
      --kingdoms bacteria \
      --ec_filter "3.2.1.23,3.2.1.37,3.2.1.85" \
      --include family subfamily ec
