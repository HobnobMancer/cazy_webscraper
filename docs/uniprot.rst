============================
Retrieving data from UniProt
============================

``cazy_webscraper`` can be used to retrieve user-specified data sets from the UniProt database, for a given subset
of proteins in a local CAZyme database created using ``cazy_webscraper``. The ``cazy_webscraper`` application can be invoked *via* the command line

-----------
Quick Start
-----------

To download UniProt protein accessions and names from UniProt for all protein in the local CAZyme database, and save the data to
the local CAZyme database, use the following command structure:  

.. code-block:: bash
    
   cw_get_uniprot_data <path to local CAZyme db>

.. NOTE::
   The ``cw`` prefix on command is an abbreviation of ``cazy_webscraper``.

.. WARNING::
    Please do not download data from UniProt for the entire CAZy database unless absolute necessary. 
    Retrieving the data from any of these exteranl databases for the entire CAZy 
    dataset will take several hours and may unintentionally deny the service to others.

To download UniProt protein accessions and names from UniProt for all protein in the local CAZyme database, including 
EC number annotations, PDB accessions and protein sequences, and save the data to
the local CAZyme database, use the following command structure:  

.. code-block:: bash
    
   cw_get_uniprot_data <path to local CAZyme db> --ec --pdb --sequence

For example:

.. code-block:: bash
    cw_get_uniprot_data cazy/cazyme.db --ec --pdb --sequence


--------------------
Command line options
--------------------

Below are listed the required and optional command-line options for configuring the retrieval of data from UniProt.

``database`` - **REQUIRED** Path to a local CAZyme database to add UniProt data to.

``--cache_dir`` - Path to cache dir to be used instead of default cache dir path.

``--cazy_synonyms`` - Path to a JSON file containing accepted CAZy class synonsyms if the default are not sufficient.

``--config``, ``-c`` - Path to a configuration YAML file. Default: None.

``--classes`` - list of classes to retrieve UniProt data for.

``--ec``, ``-e`` - Enable retrieval of EC number annotations from UniProt. Default, EC number annotations are **not** retrieved.

``--ec_filter`` - List of EC numbers to limit the retrieval of protein data for proteins annotated with at least one of the given EC numbers **in the local CAZyme database**.

``--families`` - List of CAZy (sub)families to retrieve UniProt protein data for.

``--force`` - Force writing in existing cache directory.

``--genbank_accessions`` - Path to text file containing a list of GenBank accessions to retrieve protein data for. A unique accession per line.

``--genera`` - List of genera to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given genera.

``--kingdoms`` - List of taxonomic kingdoms to restrict the scrape to. Default: None, filter is not applied.

``--log``, ``-l`` - Target path to write out a log file. If not called, no log file is written. Default: None (no log file is written out).

``--name_update`` - Boolean, whether to overwrite the existing protein name (previously retrieved from UniProt). Default: do not update.

``--nodelete_cache`` - When called, content in the existing cache dir will **not** be deleted. Default: False (existing content is deleted).

``--nodelete_log`` - When called, content in the existing log dir will **not** be deleted. Default: False (existing content is deleted).

``--pdb``, ``-p`` - Enable retrieval of PDB accessions. Default, PDB accessions not retrieved.

``--retries``, ``-r`` - Define the number of times to retry making a connection to CAZy if the connection should fail. Default: 10.

``--sequence``, ``-s`` - Enable retrieving protein amino acid sequences. Default, sequences are **not** retrieved.

``--update_name`` - If a newer version of the protein name is available, overwrite the existing name for the protein in the database. Default is false, the protein name is **not** overwritten and updated.

``--update_seq`` - If a newer version of the protein sequence is available, overwrite the existing sequence for the protein in the database. Default is false, the protein sequence is **not** overwritten and updated.

``--delete_old_ec_relationships`` - Boolean, delete old Genbanks-EC number relationships - For those proteins in the local db for whom data is downloaded from UniProt, compare the current links between the proteins in the Genbanks table and EC numbers in the Ecs table. Delete Genbanks-Ecs relationships that are not longer listed in the respective protein records in UniProt.

``--delete_old_ecs`` - Boolean, delete EC number - Delete EC numbers that are not linked to any proteins listed in the Genbanks table. These can arise from multiple retrievals of data from the UniProt data over a period of time during UniProt records have been updated.

``--delete_old_pdb_relationships`` - Boolean, delete old Genbanks-PDB relationships - For those proteins in the local db for whom data is downloaded from UniProt, compare the current links between the proteins in the Genbanks table and PDB accessions in the Pdbs table. Delete Genbanks-Pdbs relationships that are not longer listed in the respective protein records in UniProt.

``--delete_old_pdbs`` - Boolean, delete PDB accessions - Protein relationships that are no longer listed in UniProt, i.e. an PDB accessions that are no longer included in UniProt but is in the local database. If set to TRUE these relationships will be DELETED from the database.

``--use_uniprot_cache`` - Path to a JSON file, keyed by UniProt accessions/IDs and valued by dicts containing `{'gbk_acc': str, 'db_id': int}`. This file part of the cache created by `cw_get_uniprot_data`. This is option to skip retrieving the UniProt IDs for a set of GenBank accessions, if retrieving data for the same dataset (this save a lot of time!)

``skip_download`` - Bool, default False. If set to True, only uses data from UniProt cache and will not download new data from UniProt.

``--sql_echo`` - Set SQLite engine echo parameter to True, causing SQLite to print log messages. Default: False.

``--species`` - List of species (organsim scientific names) to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given species.

``--strains`` - List of species strains to restrict the retrieval of protein to data from UniProt to proteins belonging to one of the given strains.

``--taxonomy``, ``-t`` - Retrieve taxonomic classifications (genus species) and add to the local CAZyme db.

``--timeout`` - Connection timout limit (seconds). Default: 45.

``--use_uniprot_cache`` - Path to JSON file containing data previosuly retrieved from UniProt by ``cazy_webscraper``, use if an error occurred while adding the data to the local CAZyme database. This will skip the retrieval of data from UniProt, and the cached data will be added to the local CAZyme database. This can also be shared with others to add the same data to their local CAZyme database.

``--bioservices_batch_size`` - Size of an individual batch query submitted to the `UniProt REST API <https://www.uniprot.org/help/programmatic_access>_` to retrieve protein data from UniProt. Default is 1000. 

``--verbose``, ``-v`` - Enable verbose logging. This does **not** set the SQLite engine ``echo`` parameter to True. Default: False.

-----------
Batch sizes
-----------

Note that according to Uniprot (June 2022), there are various limits on ID Mapping Job Submission:

========= =====================================================================================
Limit	  Details
========= =====================================================================================
100,000	  Total number of ids allowed in comma separated param ids in /idmapping/run api
500,000	  Total number of "mapped to" ids allowed
100,000	  Total number of "mapped to" ids allowed to be enriched by UniProt data
10,000	  Total number of "mapped to" ids allowed with filtering
========= =====================================================================================

-----------
Basic Usage
-----------

The command-line options listed above can be used in combination to customise the retrieval of protein data from UniProt. Some options (e.g. ``--families`` and ``--classes``) define the broad group of proteins for which data will be retrieved from UniProt, others (e.g. ``--species``) are used to filter and fine-tune the protein dataset for which protein data will be retrieved.

The ``--classes``, ``--families``, ``--kingdoms``, ``--genera``, ``--species``, and ``--strains`` filteres are applied 
in the exactly same for retrieving data from CAZy as retrieving data from UniProt. Examples of using these flags 
can be found in the ``cazy_webscraper`` tutorial in this documentation.

Here we discuss using the new flags ``--ec``, ``--pdb``, ``--sequence``, ``--seq_update``, and ``--ec_filter``.

.. NOTE::
    To retrieve data for members of specific CAZy subfamilies, list the subfamilies after the ``--families`` 
    flag.

.. NOTE::
    The command for retrieving protein data from UniProt for proteins in a local CAZyme database is ``cw_get_uniprot_data``.

-----------------------------
Data retrievable from UniProt
-----------------------------

By default ``cw_get_uniprot_data`` retrieves the UniProt protein accession and protein name from UniProt, for proteins in a 
local CAZyme database. ``cw_get_uniprot_data`` can also retrieve from UniProt:

* EC number annotations
* PDB accessions
* Protein amino acid sequences

.. WARNING:: 
    When performing retrievals of large datasets from UniProt, please perform these retrievals during quiet periods (e.g. at the weekend). 

.. WARNING::
    It is strongly advised to only download data from UniProt for the necessary protein datasets. Retrieval of very large 
    datasets from UniProt (e.g. for +1,000,000 proteins or all proteins in the GH class) can result in UniProt 
    terminating the connection early due to a high bandwidth command over an extended period of time.

    To retrieve data from large datasets, it is recommend to break down the dataset into subgroups, and periodically 
    retrieve the data for one subgroup. This reduces the burdeon on the UniProt server and will not break 
    the expected use and practises of the UniProt database.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Retrieving EC number annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To retrieve EC number annotations from UniProt add the ``--ec`` /  ``-e`` flag to the command:

.. code-block:: bash

    cw_get_uniprot_data cazy_db.db --ec

OR

.. code-block:: bash

    cw_get_uniprot_data cazy_db.db -e

.. NOTE::
    **All** EC number annotations are retrieved for all CAZymes matching the given filter criteria. In the example 
    command above, no filters were provided therefore, all EC number annotations will be retrieved for all
    CAZymes in the local CAZyme database (in this case called ``cazy_db.db``).


^^^^^^^^^^^^^^^^^^^^^^^^^
Retrieving PDB accessions
^^^^^^^^^^^^^^^^^^^^^^^^^

To retrieve all PDB accessions for all CAZymes in the local CAZyme database matching the given filter criteria,
add the ``--pdb`` /  ``-p`` flag to the command:

.. code-block:: bash

    cw_get_uniprot_data cazy_db.db --pdb

OR

.. code-block:: bash

    cw_get_uniprot_data cazy_db.db -p


^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Retrieving protein sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To retrieve all protein amino acid sequences for all CAZymes in the local CAZyme database matching the given filter criteria,
add the ``--sequence`` /  ``-s`` flag to the command:

.. code-block:: bash

    cw_get_uniprot_data cazy_db.db --sequence

OR

.. code-block:: bash

    cw_get_uniprot_data cazy_db.db -s

``cw_get_uniprot_data`` stores the protein amino acids sequence within the local CAZyme database, as well 
as the 'last modified date' retrieved from UniProt.


^^^^^^^^^^^^^^^^^^^^^^^^
Updating local sequences
^^^^^^^^^^^^^^^^^^^^^^^^

When using ``--sequence`` flag, ``cw_get_uniprot_data`` will only add *new* protein sequences to the database, i.e.
it will only add protein sequences to records that do not have a sequence. Therefore, if a protein
already has a sequence in the local database, this sequence is **not** overwritten.

You may wish to update the protein sequences in your local CAZyme database. To do this use the ``--sequence`` /  ``-s`` 
flag to tell ``cw_get_uniprot_data`` to retrieve protein sequences, **and** use the ``--seq_update`` flag.

.. code-block:: bash

    cw_get_uniprot_data cazy_db.db -s --seq_update

This instructs ``cw_get_uniprot_data`` to overwriting existing protein sequences in the local database *if* a newer version 
of the sequence is retrieved from UniProt. This is checked by comparing the 'last modified date' of the 
protein sequence in the local database against the sequence retrieved from UniProt.


--------------------------
Using the EC number filter
--------------------------

Having previously retrieved EC number annotations and added them to the local CAZyme database, you  may 
wish to retrieve protein data for CAZymes annotated with specific EC numbers. To do this add the 
``--ec_filter`` flag to the command, follwed by a list of EC numbers.

.. NOTE::
    Provide complete EC numbers. 
    Both dashes ('-') and asterixes ('*') are accepted for missing digits in EC numbers.

    EC1.2.3.- and EC1.2.3.* are accepted.
    EC1.2.3. and EC 1.2.3 are **not** accepted.

.. WARNING::
    If using dashes to represent missing digits in EC numbers, it is recommended to bookend the entire 
    EC number list in single or double quotation marks. Some terminals may misinterpret EC1.2.-.- as trying to invoke the options '.'

.. NOTE::
    ``cw_get_uniprot_data`` will retrieve the specified UniProt data for all proteins in the local CAZyme 
    database that are annotated with **at least one** of the given EC numbers. Therefore, if multiple 
    EC numbers are given this **does not mean** only CAZymes will all provided EC numbers will have data retrieved
    from UniProt for them.

``--ec_filter`` is based upon EC number annotations stored within the local CAZyme database. For 
example, if protein A is annotated with the EC1.2.3.4, but this annotation is not stored in the 
local CAZyme database, using ``--ec_filter EC1.2.3.4`` will **not** cause ``cw_get_uniprot_data`` to retrieve
data for protein A. This is because ``cw_get_uniprot_data`` does not know protein A is annotated with 
EC1.2.3.4, because this annotation is not within its database.

.. WARNING::
    If ``--ec_filter`` is used along side ``--ec``, ``cw_get_uniprot_data`` will retrieve **all** EC number 
    annotations from UniProt for all proteins in the local CAZyme database that are associated with 
    at least one of the EC numbers provided via ``--ec_filter`` within the CAZyme database.

-------------------------------
Configuration using a YAML file
-------------------------------

As with scraping CAZy, a YAML file can be provided to define the filters for retrieving data from UniProt. 
The same YAML file can be used both for scraping CAZy and UniProt. However, the configuration file for
retrieving data from UniProt can contain the additionl ``ec`` tag.

Using a config file supports reproducible documentation of ``cazy_webscraper`` usage.

An template YAML file is provided in the ``cazy_webscraper`` repository (``configuration_files/template-get_data_config.yaml``):

.. code-block:: yaml

    # Under 'classes' list class from which all proteins will retrieved
    # Under each families respective name, list the specific families/subfamilies to be scraped
    # Write the FULL family name, e.g. 'GH1', NOT only its number, e.g. '1'
    # To list multiple families, each familiy must be on a new line starting indented once
    # relative to the parent class name, and the name written within quotation marks.
    # For more information on writing lists in Yaml please see:
    # https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html 
    classes:  # classes from which all proteins will be retrieved
    - "GH"
    - "CE"
    Glycoside Hydrolases (GHs):
    GlycosylTransferases (GTs):
    Polysaccharide Lyases (PLs):
    - "GT1"
    - "GT5"
    - "GT6"
    Carbohydrate Esterases (CEs):
    Auxiliary Activities (AAs):
    Carbohydrate-Binding Modules (CBMs):
    genera:  # list genera to be scraped
    - "Trichoderma"
    - "Aspergillus"
    species:  # list species, this will scrape all strains under the species
    - "Pythium ultimum"
    strains:  # list specific strains to be scraped
    kingdoms:  # Archaea, Bacteria, Eukaryota, Viruses, Unclassified
    ec:
    - "EC1.2.3.4"

.. ATTENTION::
  The YAML configuration file must contain all tags/headings indicated in the example configuration file found in the repository:

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
  * ec

Each value in the YAML mappings for these arguments must be listed on a separate line, indented by 4 spaces, and the class name encapsulated with single or double quotation marks. For example:

.. code-block:: yaml

    classes:
        - "GT"
        - "pl"
    Glycoside Hydrolases (GHs):
        - "GH1"
        - "GH2"


^^^^^^^^^^^^^^^^^^^^^^^^^
Synonyms for CAZy classes
^^^^^^^^^^^^^^^^^^^^^^^^^

A number of synonyms may be provided for CAZy classes, e.g. both "GH" and "Glycoside-Hydrolases" are accepted as synonyms for "Glycoside Hydrolases (GHs)" (the name recorded at CAZy). These alternatives are defined in the ``cazy_webscraper`` repository, in the file ``scraper/utilities/parse_configuration/cazy_dictionary.json``.
