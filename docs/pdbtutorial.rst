==========================================================================
Tutorials on configuring ``cazy_webscraper`` to retrieve data from PDB
==========================================================================

``cazy_webscraper`` can be configured to retrieve protein structures files user specified sets of 
CAZymes in a local CAZyme database. Many of the same configuration options 
apply to the retrieval of protein data from CAZy, UniProt, GenBank and PDB.

``BioPython`` is used to perform the retrieval of protein structure files from PDB. The retrieved structure 
files are written to the disk, *not* into the local CAZyme database.

> Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., … others. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423.

``cazy_webscraper`` can be configured via the **command line** and/or via a **YAML configuration file**.

This page runs through examples of how to combine the various 'filters' that can be applied, to fully customised 
the retrieval of protein structured from PDB. These tutorials are designed for those with less experience using command-line tools.

.. NOTE::
  If you installed ``cazy_webscraper`` using ``bioconda`` or ``pip`` to invoke ``cazy_webscraper`` to retrieve UniProt data call it using ``cw_get_pdb_structures`` - this is the method used in this tutorial.  
  If you installed ``cazy_webscraper`` from source then you will need to invoke ``cazy_webscraper`` from the root of the repo using the command ``python3 cazy_webscraper/expand/uniprot/get_uniprot_data.py``.

From this point on, we will be discussed the ``cw_get_pdb_structures``, which is the entry point for 
retrieving data from UniProt. We also presume you are comfortable configuring ``cazy_webscraper`` for the 
scraping of data from CAZy.

----------------------------------
Configuration via the command line
----------------------------------

``cw_get_pdb_structures`` has two required argument:
* The path to the local CAZyme database created using ``cazy_webscraper``
* The structure file formats to retrieve the data from PDB in

The accepted structure file formats are:
* mmCif (default, PDBx/mmCif file),
* pdb (format PDB),
* xml (PDBML/XML format),
* mmtf (highly compressed),
* bundle (PDB formatted archive for large structure}

Any combination of file formats can be provided to ``cw_get_pdb_structures`` to download every file type for each PDB accession in the local CAZyme database. To list multiple file formats, separate each file format with a single space (' '). For example, to download the mmCif and xml files for every PDB accession in a local CAZyme database (located at ``cazy/cazyme_db.db``), use the following command:

.. code-block:: bash
    
    cw_get_pdb_structures cazy/cazyme_db.db mmCif xml

.. WARNING::
    The file formats are case sensitive. For example, make sure to use 'mmCif' not 'mmcif'.

When no optional arguments are provided, the default behaviour is invoked. The default behaviour is to: 
Retrieve protein data for **all** CAZymes in the local CAZyme db

-----------------------------------------
Options configurable at the command line 
-----------------------------------------

CAZymes of interest can be defined via providing:

* A set of GenBank accessions
* A set of UniProt accessions
* CAZy classes
* CAZy families
* Taxonomic kingdoms
* Genera
* Species
* Strains
* EC numbers

`Here <https://cazy-webscraper.readthedocs.io/en/latest/pdb.html>`_ you can find a full list of the command-line flags and options.


--------------------------------------------------------------------
Retrieving protein structures for specific CAZy classes and families
--------------------------------------------------------------------

The ``--classes`` and ``--families`` flags from scraping data from CAZy are applied in the extact same way 
for retrieving protein structure files from PDB.

For instance, if instead of retrieving protein data for all CAZymes in your local CAZyme database, you want to 
retrieve protein data for CAZymes in specific CAZy classes then add the 
``--classes`` flag followed by the classes you want to retrieve protein data for.

.. TIP::
   To list multiple classes, separate the classes with a single comma. 

For example, if you want to retrieve protein data for all CAZymes from Glycoside Hydrolase and Carbohydrate Esterases then use the command:

.. code-block:: bash

   cw_get_pdb_structures cazy/cazyme.db mmCif --classes GH,CE

OR

.. code-block:: bash

   cw_get_pdb_structures cazy/cazyme.db mmCif --classes Glycoside Hydrolases,Carbohydrate Esterases

Retrieving protein data for proteins from specific specific CAZy families is achieved using the ``--families`` flag. For 
example, to retrieve protein data for all proteins in PL1, PL2 and PL3 in the local CAZyme database, in mmCif and PDB format, use the 
following command:

.. code-block:: bash

   cw_get_pdb_structures cazy/cazyme.db mmCif pdb --families PL1,PL2,PL3

.. WARNING::
   ``cw_get_pdb_structures`` only accpets families written in the proper CAZy family syntax.
   GH1 is accepted.
   gh1 and GlycosideHydrolases1 are not accepted.

As with scraping data from CAZy, the ``--classes`` and ``--families`` flags can be combined. To retrieve 
protein data for all CAZymes in PL1, PL2, PL3 and *all* of GH and CE both:

.. code-block:: bash

   cw_get_pdb_structures cazy/cazyme.db mmCif pdb --families PL1,PL2,PL3 --classes GH,CE

**AND**

.. code-block:: bash

   cw_get_pdb_structures cazy/cazyme.db mmCif pdb --classes GH,CE --families PL1,PL2,PL3

are accepted.


------------------
Applying taxonomic
------------------

The ``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` flags can be used to refine the dataset 
of proteins to retrieve protein data by taxonomy. These flags are applied in the exact same way as they 
are used for the scraping of data from CAZy. Only proteins in the local CAZyme database and 
matching at least on of the provided taxonomy criteria will have protein data retrieved from PDB.

For example, if you want to retrieve protein data for all CAZymes in a local CAZyme database from bacterial and eukaryotic species, in pdb and xml formats, then use the command 

.. code-block:: bash

   cw_get_pdb_structures cazy/cazyme.db pdb xml --kingdoms bacteria,eukaryota

.. warning::
   The kingdoms must be spelt the same way CAZy spells them, for example use 'eukaryot**a**' instead of 'eukaryot**e**'.
   
.. NOTE:: 
   The kingdoms are **not** case sensitive, therefore, both ``bacteria`` *and* ``Bacteria`` are accepted. 

.. NOTE::
   You can list the kingdoms in *any* order. Thus, both ``bacteria,eukaryota`` *and* ``eukaryota,bacteria`` are accepted.

You can combine any combination of the optional flags, including combining the taxonomic filters. For example,
you may wish to retrieve protein data for all CAZymes in a local CAZyme database that are derived from all viral species, Aspergillus species, Layia carnosa, Layia chrysanthemoides, Trichoderma reesei QM6a and 
Trichoderma reesei QM9414. To do this we would combine the respective flags for a single ``cw_get_pdb_structures`` command. The command 
we would use would be:

.. code-block:: bash

   cw_get_pdb_structures cazy/cazyme.db pdb xml --kingdoms viruses --genera Aspergillus --species Layia carnosa,Layia chrysanthemoides --strains Trichoderma reesei QM6a,Trichoderma reesei QM9414

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
   When you specify a species ``cw_get_pdb_structures`` will retrieval CAZymes from *all* strains of the species.


-------------------------
Applying EC number filter
-------------------------

The retrieval of protein data from PDB can also be limited to proteins in a local CAZyme database that are
annotated with specific EC numbers.

Having previously retrieved EC number annotations from UniProt and adding them to the local CAZyme database, you may 
wish to retrieve protein data for CAZymes annotated with specific EC numbers. To do this add the 
``--ec_filter`` flag to the command, follwed by a list of EC numbers.

.. code-block:: bash
   
   cw_get_pdb_structures cazy/cazyme.db pdb --ec_filter "EC1.2.3.4,EC2.3.4.5"


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
    ``cw_get_pdb_structures`` will retrieve the PDB structure files for all proteins in the local CAZyme 
    database that are annotated with **at least one** of the given EC numbers. Therefore, if multiple 
    EC numbers are given this **does not mean** PDB structure files will only be retrieved for 
    CAZymes annotated for all provided EC numbers.

``--ec_filter`` is based upon EC number annotations stored within the local CAZyme database. For 
example, if protein A is annotated with the EC1.2.3.4, but this annotation is not stored in the 
local CAZyme database, using ``--ec_filter EC1.2.3.4`` will **not** cause ``cw_get_pdb_structures`` to retrieve
data for protein A. This is because ``cw_get_pdb_structures`` does not know protein A is annotated with 
EC1.2.3.4, because this annotation is not within its database.

.. WARNING::
    If ``--ec_filter`` is used along side ``--ec``, ``cw_get_pdb_structures`` will retrieve **all** EC number 
    annotations from UniProt for all proteins in the local CAZyme database that are associated with 
    at least one of the EC numbers provided via ``--ec_filter`` within the CAZyme database.


---------------------
Combining all filters
---------------------

The ``--classes``, ``--families``, ``--ec_filter``, ``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` flags can 
be used in any combination to define a specific subset of proteins in the local CAZyme database for whom
protein structure files will be retrieved from PDB. These flags can be used with any combination of 
``--ec``, ``--pdb``, ``--sequence``, ``--update_seq`` to customise what data is retrieved from UniProt and 
added to the local CAZyme database.

Below we run through 3 example commands of combining these flags, and the resulting behaviour.

**Example 1:**
To retrieve PDB structure file for all CAZymes in GH, GT, CE1, CE5 and CE8, and which are derived from baceterial species, and retrieve the files in mmCif and xml format we use the command:

.. code-block:: bash

   cw_get_pdb_structures cazy/cazyme.db mmCif xml --classes GH,CE --families CE1,CE5,CE8 --kingdoms bacteria


**Example 2:**
To retrieve PDB structure files for all CAZymes in GH and which are derived from *Aspegillus* and *Trichoderma* species in bundle format we use the command:

.. code-block:: bash

   cw_get_pdb_structures cazy/cazyme.db bundle --classes GH --genera Aspegillus,Trichoderma


**Example 3:**
To retrieve PDB structure files for all CAZymes in GH,CE and CBM which are derived from baceterial species and are annotated with at least one of 
EC3.2.1.23, EC3.2.1.37 and EC3.2.1.85, and retrieve the files in pdb and xml format we use the command:

.. code-block:: bash

   cw_get_pdb_structures cazy/cazyme.db xml pdb --classes GH,CE,CBM --kingdoms bacteria --ec_filter "3.2.1.23,3.2.1.37,3.2.1.85"

.. NOTE::
   The order the structure file formats are provided does **not** matter.

------------------------------
Providing a list of accessions
------------------------------

Instead of retrieving protein structure files for all CAZymes matching a defined set of criteria, 
``cw_get_pdb_structures`` can retrieve protein structure files for a set of CAZymes defined by their 
GenBank and/or UniProt accession.

The flag ``--genbank_accessions`` can be used to provide ``cw_get_pdb_structures`` a list of GenBank accessions 
to identify the specific set of CAZymes to retrieve protein structure files for.

The flag ``--uniprot_accessions`` can be used to provide ``cw_get_pdb_structures`` a list of UniProt accessions 
to identify the specific set of CAZymes to retrieve protein structure files for.

In both instances (for ``--genbank_accessions`` and ``--uniprot_accessions``) the list of respective accessions 
are provided via a plain text file, with a unique protein accession of each line. The path to this file is 
then passed to ``cw_get_pdb_structures`` via the respective ``--genbank_accessions`` and ``--uniprot_accessions`` flag.

``--genbank_accessions`` and ``--uniprot_accessions`` can be used at the same time to define all 
CAZymes of interest.

.. WARNING::
   ``--genbank_accessions`` and ``--uniprot_accessions`` take president over the filter flags.

   When either ``--genbank_accessions`` or ``--uniprot_accessions`` is used, ``cw_get_pdb_structures`` will 
   **not** retrieve any CAZymes from the local database matching a set of criteria.

   Therefore, if ``--genbank_accessions`` and ``--classes`` are used, ``cw_get_pdb_structures`` will ignore 
   the ``--classes`` flag and only retrieve PDB structure files for the proteins listed in the file provided via 
   the ``--genbank_accessions``.
