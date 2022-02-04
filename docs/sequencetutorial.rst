============================================================
Tutorials on configuring the extraction of protein sequences
============================================================

``cazy_webscraper`` can be configured to extract GenBank and/or UniProt protein sequences for user specified sets of proteins from 
a local CAZyme database. Many of the configuration options 
apply to the retrieval of protein data from CAZy, UniProt, GenBank and PDB.

``cazy_webscraper`` can be configured via the **command line** and/or via a **YAML configuration file**.

This page runs through examples of how to combine the various 'filters' that can be applied, to fully customised 
the extraction of protein sequences. These tutorials are designed for those with less experience using command-line tools.

.. NOTE::
  If you installed ``cazy_webscraper`` using ``bioconda`` or ``pip`` to invoke ``cazy_webscraper`` to retrieve UniProt data call it using ``cw_extract_db_sequences`` - this is the method used in this tutorial.  
  If you installed ``cazy_webscraper`` from source then you will need to invoke ``cazy_webscraper`` from the root of the repo using the command ``python3 cazy_webscraper/expand/extract/extract_sequences.py``.

From this point on, we will be discussed the ``cw_extract_db_sequences``, which is the entry point for 
extract protein sequences from the local CAZyme database. We also presume you are comfortable configuring ``cazy_webscraper`` for the 
scraping of data from CAZy.


----------------------------------
Configuration via the command line
----------------------------------

``cw_extract_db_sequences`` has at 2 required arguments:

1. The path to the local CAZyme databases created using ``cazy_webscraper``
2. The names of the database from which the proteins were sourced

When no optional arguments are provided, the default behaviour is invoked. The default behaviour is to: 
Extract (GenBank and/ot UniProt) protein sequences for **all** CAZymes in the local CAZyme db


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Defining the source of the proteins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Protein sequences previously retrieved from GenBank **and/or** UniProt can be extracted.

* To extract only protein sequences from GenBank, use ``genbank``.  
* To extract only protein sequences from UniProt, use ``uniprot``.
* To extract GenBank and UniProt protein sequences, use ``genbank uniprot`` or ``uniprot genbank``.

.. TIP::
    The order the databases (i.e. 'genbank' and 'uniprot') does not matter, and they are **not** casesensitive.

``cw_extract_db_sequences`` can be enabled using a simple command structure:

.. code-block:: bash

  cw_extract_db_sequences <path to the local CAZyme db> genbank uniprot

For example, if our database was stored in ``cazy/cazyme.db``, we would used:

.. code-block:: bash
   
  cw_extract_db_sequences cazy/cazyme.db genbank uniprot


^^^^^^^^^^^^^^^^^^^^^^^^^^
Defining the target output
^^^^^^^^^^^^^^^^^^^^^^^^^^

The extracted protein sequences can be written to any combination of the following:
* A single FASTA file containing all extracted protein sequences
* One FASTA file per extracted sequence
* A BLAST database

To define a **single FASTA file** to write all extracted sequences to use the ``--fasta_file`` flag, followed by the 
path to the target file. These file nor its parent directorties need to already exist, ``cazy_webscraper`` will build 
all necessary parent directories. For example:  

.. code-block:: bash

    cw_extract_db_sequences cazy/cazyme_database.db genbank --fasta_file cazy/protein_sequences/all_genbanks.FASTA

To write out each extracted sequence **to its own FASTA file** use the ``--fasta_dir`` flag, followed by the 
path to the target directory. ``cazy_webscraper`` will build 
all necessary parent directories. For example:  

.. code-block:: bash

    cw_extract_db_sequences cazy/cazyme_database.db genbank --fasta_dir cazy/protein_sequences

To build a BLAST database of the extracted protein squences, use the ``--blastdb`` or ``-b`` flag, followed by the path 
of where to write out the database, including the database name. ``cazy_webscraper`` will build 
all necessary parent directories. For example:  

.. code-block:: bash

    cw_extract_db_sequences cazy/cazyme_database.db genbank --blastdb cazy/protein_sequences/all_genbanks_blast_db.db

Any combination of ``--fasta_file``, ``--fasta_dir`` and ``--blastdb`` can be used to produce multiple outputs. For example, 
to generate a single FASTA file of all extracted UniProt protein sequences **and** write the extracted sequences to a BLAST database:  

.. code-block:: bash

    cw_extract_db_sequences cazy/cazyme_database.db genbank \
        --fasta_file cazy/protein_sequences/all_genbanks.FASTA \
        --blastdb cazy/protein_sequences/all_genbanks_blast_db.db

.. TIP::
    Backward slashes '\' can be used to break up a long command into multiple lines to make it easier to read.


------------------
FASTA file formats
------------------

The FASTA files generated by ``cw_extract_db_sequences`` have a very simple protein ID line. The line always and only contains:
* The GenBank or UniProt accession
* The name of the source database: 'GenBank' or 'UniProt'

For example, a protein sequence from GenBank which is extracted from a local CAZyme datbase will be presented as:

.. code-block:: bash
    > AIP21820.1 GenBank
    MPVALAVAAALGACSGDDDATLESRADAIVERMTTRQKVGQKLMMAFRYWCPDGQPACTT
    GMTEFPDAARDALRENGIGGVILFSNNLTGIEQTRRLIDGIRAAPAADSPLGLMIGIDEE
    GGNVFRLPRVEATAFAGNMALGAAYEATRDDRLAYDMGRVLAAEIAAVGFNVNFAPDVDV
    NSNPLNPVINVRAFGDDPATIGLLGRRMVQGMKSERVIGTFKHFPGHGDTDTDSHYGLPV
    VIKSRADAYAIDLAPYRQAIEAGEAPDMIMTAHIQYPSLDDTRVATRTGEQMIAPATMSR
    RIQHDILRGEFGYQGVTITDALDMKGIAGFFDEDDAVVKVFQADVDIALMPVEFRTAADA
    GRLAALVDRVAAAVDSGRIDRAEFDRSVRRIVLTKLRHGIVASDRGRPVDELASIGGPAH
    RAIERDIAQKSITVLRNENGALPLQAAGRRIFILTPWGEQAEAMRRRFVELGHPLVTGAK
    LSAITWAEQQQAIDAADVVIVGTLSTGVTPVEHNGDPNARVSPPAPSAVRMRQAAPANGE
    EEGSVIFDHVERADAAKDIGARPSVLAAIAAPSEAQQMRDAMDYAKARRKTVIHVTMRAP
    YDVISYDDVADATLATYAYYGYEGGLRGPSLPAAVDAMLGVGRPVGRLPVAIHALNADGS
    TGPLRYARGFGLQY


-----------------------------------------
Options configurable at the command line 
-----------------------------------------

The following behaviours of the ``cw_extract_db_sequences`` can be configured at the command-line in the terminal to 
limit the extraction of protein sequences to CAZymes in the local databaes from specific:

* CAZy classes
* CAZy families and subfamilies
* Taxonomic kingdoms
* Genuera
* Species
* Species strains
* Annotated with at least one of a set of specified EC numbers

`Here <https://cazy-webscraper.readthedocs.io/en/latest/sequence.html>`_ you can find a full list of the command-line flags and options.


----------------------------------------------------------------
Extract protein sequences for specific CAZy classes and families
----------------------------------------------------------------

The ``--classes`` and ``--families`` flags from scraping data from CAZy are applied in the extact same way 
for extracting protein sequences for proteins of interest.

For instance, if instead of extracting protein sequences for all CAZymes in your local CAZyme database, you want to 
extract protein sequebces for CAZymes in specific CAZy classes then add the 
``--classes`` flag followed by the classes you want to extract protein sequences for.

.. TIP::
   To list multiple classes, separate the classes with a single comma. 

For example, if you want to extract protein sequences for all CAZymes from Glycoside Hydrolase and Carbohydrate Esterases then use the command:

.. code-block:: bash

   cw_extract_db_sequences cazy/cazyme.db genbank --classes GH,CE

OR

.. code-block:: bash

   cw_extract_db_sequences cazy/cazyme.db genbank --classes 'Glycoside Hydrolases','Carbohydrate Esterases'

.. WARNING::
    When including spaces in a parameter value, such as 'Glycoside Hydrolases' single or double quotation marks must be written around the value.

Extracting protein sequences for proteins from specific specific CAZy families is achieved using the ``--families`` flag. For 
example, to extract GenBank protein sequences for all proteins in PL1, PL2 and PL3 in the local CAZyme database use the 
following command:

.. code-block:: bash

   cw_extract_db_sequences cazy/cazyme.db genbank --families PL1,PL2,PL3

.. WARNING::
   ``cw_extract_db_sequences`` only accpets families written in the proper CAZy family syntax.
   GH1 is accepted.
   gh1 and GlycosideHydrolases1 are not accepted.

As with scraping data from CAZy, the ``--classes`` and ``--families`` flags can be combined. To extract UniProt protein sequences 
for all CAZymes in PL1, PL2, PL3 and *all* of GH and CE both:

.. code-block:: bash

   cw_extract_db_sequences cazy/cazyme.db uniprot --families PL1,PL2,PL3 --classes GH,CE

**AND**

.. code-block:: bash

   cw_extract_db_sequences cazy/cazyme.db uniprot --classes GH,CE --families PL1,PL2,PL3

are accepted.


------------------
Applying taxonomic
------------------

The ``--kingdoms``, ``--genera``, ``--species`` and ``--strains`` flags can be used to refine the dataset 
of proteins to extract protein sequences by taxonomy. These flags are applied in the exact same way as they 
are used for the scraping of data from CAZy. Only proteins in the local CAZyme database and matching at least on of the provided taxonomy 
criteria will have protein data retrieved from UniProt and added to the local CAZyme datbase.

For example, if you want to extract GenBank protein sequences for all CAZymes in a local CAZyme database from bacterial and eukaryotic species then use the command 

.. code-block:: bash

   cw_extract_db_sequences cazy/cazyme.db genbank --kingdoms bacteria,eukaryota

.. warning::
   The kingdoms must be spelt the same way CAZy spells them, for example use 'eukaryot**a**' instead of 'eukaryot**e**'.
   
.. NOTE:: 
   The kingdoms are **not** case sensitive, therefore, both ``bacteria`` *and* ``Bacteria`` are accepted. 

.. NOTE::
   You can list the kingdoms in *any* order. Thus, both ``bacteria,eukaryota`` *and* ``eukaryota,bacteria`` are accepted.

You can combine any combination of the optional flags, including combining the taxonomic filters. For example,
you may wish to extract GenBank and UniProt protein sequences for all CAZymes in a local CAZyme database that are derived from all viral species, Aspergillus species, Layia carnosa, Layia chrysanthemoides, Trichoderma reesei QM6a and 
Trichoderma reesei QM9414. To do this we would combine the respective flags for a single ``cw_extract_db_sequences`` command. The command 
we would use would be:

.. code-block:: bash

   cw_extract_db_sequences cazy/cazyme.db genbank uniprot --kingdoms viruses --genera Aspergillus --species Layia carnosa,Layia chrysanthemoides --strains Trichoderma reesei QM6a,Trichoderma reesei QM9414

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
   When you specify a species ``cw_extract_db_sequences`` will retrieval CAZymes from *all* strains of the species.


-------------------------
Applying EC number filter
-------------------------

The extraction of protein sequences an also be limited to proteins in a local CAZyme database that are
annotated with specific EC numbers.

Having previously retrieved EC number annotations from UniProt and added them to the local CAZyme database, you  may 
wish to extract protein sequences for CAZymes annotated with specific EC numbers. To do this add the 
``--ec_filter`` flag to the command, follwed by a list of EC numbers.

.. code-block:: bash
   
   cw_extract_db_sequences cazy/cazyme.db genbank --ec_filter "EC1.2.3.4,EC2.3.4.5"


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
protein sequences will be extracted. These flags can be used with any combination of 
``--ec``, ``--pdb``, ``--sequence``, ``--update_seq`` to customise what data is retrieved from UniProt and 
added to the local CAZyme database.

Below we run through 3 example commands of combining these flags, and the resulting behaviour.

**Example 1:**
To extract GenBank protein sequences for CAZymes:
* In GH, GT, CE1, CE5 and CE8
* Derived from bacterial species

.. code-block:: bash

   cw_extract_db_sequences cazy/cazyme.db genbank --classes GH,CE --families CE1,CE5,CE8 --kingdoms bacteria


**Example 2:**
To extract GenBank protein sequences for CAZymes:
* In GH
* From *Aspegillus* and *Trichoderma* species
.. code-block:: bash

   cw_extract_db_sequences cazy/cazyme.db genbank --classes GH --genera Aspegillus,Trichoderma


**Example 3:**
To extract GenBank and UniProt protein sequences for CAZymes:
* In GH,CE and CBM
* Derived from baceterial species
* Annotated with at least one of EC3.2.1.23, EC3.2.1.37 and EC3.2.1.85

.. code-block:: bash

   cw_extract_db_sequences cazy/cazyme.db genbank uniprot --classes GH,CE,CBM --kingdoms bacteria --ec_filter "3.2.1.23,3.2.1.37,3.2.1.85"

------------------------------
Providing a list of accessions
------------------------------

Instead of extracting protein sequences for all CAZymes matching a defined set of criteria, 
``cw_extract_db_sequences`` can extract protein sequences a set of CAZymes defined by their 
GenBank and/or UniProt accession.

The flag ``--genbank_accessions`` can be used to provide ``cw_extract_db_sequences`` a list of GenBank accessions 
to identify the specific set of CAZymes to extract protein sequences for.

The flag ``--uniprot_accessions`` can be used to provide ``cw_extract_db_sequences`` a list of UniProt accessions 
to identify the specific set of CAZymes to extract protein sequences for.

In both instances (for ``--genbank_accessions`` and ``--uniprot_accessions``) the list of respective accessions 
are provided via a plain text file, with a unique protein accession of each line. The path to this file is 
then passed to ``cw_extract_db_sequences`` via the respective ``--genbank_accessions`` and ``--uniprot_accessions`` flag.

``--genbank_accessions`` and ``--uniprot_accessions`` can be used at the same time to define all 
CAZymes of interest.

The ``sources`` of the proteins operates independently of the ``--genbank_accessions`` and ``--uniprot_accessions`` 
flags. Therefore, the ``--uniprot_accessions`` flag can be used to identify a set of CAZymes of interest 
by their UniProt accession, and their protein sequence ``source`` can be defined as 'genbank', which will 
retrieve the GenBank protein sequence for the specified CAZymes of interest.

.. WARNING::
   ``--genbank_accessions`` and ``--uniprot_accessions`` take president over the filter flags.

   When either ``--genbank_accessions`` or ``--uniprot_accessions`` is used, ``cw_extract_db_sequences`` will 
   **not** retrieve any CAZymes from the local database matching a set of criteria.

   Therefore, if ``--genbank_accessions`` and ``--classes`` are used, ``cw_extract_db_sequences`` will ignore 
   the ``--classes`` flag and only extract protein squences for the proteins listed in the file provided via 
   the ``--genbank_accessions``.
