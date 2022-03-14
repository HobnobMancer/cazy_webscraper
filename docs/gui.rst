============================
The Graphical User Interface
============================

``cazy_webscraper`` includes a graphical user interface (GUI) that individaully wraps each of the 
subcommands in ``cazy_webscraper``. Therefore, a separate GUI is built for:
* Scraping CAZy
* Retrieving protein data from UniProt
* Retrieving protien sequences from GenBank
* Extracting protein sequences from a local CAZyme database
* Retrieving structure files from PDB
* Interrogating and extracting information from a local CAZyme database

-----------
Basic Usage
-----------

The GUIs are called from the command-line. The table below lists the commands for invoking each 
subcommand/feature of ``cazy_webscraper`` when ``cazy_webscraper`` was installed using ``pip`` or ``bioconda`` 
or installed from source.

.. note::
    These commands provided for when ``cazy_webscraper`` was installed from source, presume the current 
    working directory is the root of the ``cazy_webscraper`` repository.

+-------------------------------------+--------------------------------+--------------------------------------------------+
| Subcommand                          | Installed with pip / bioconda  | Installed from source                            |
+=====================================+================================+==================================================+
| Scrape CAZy                         | cw_gui_cazy_webscraper         | python3 cazy_webscarper/gui/gui_cazy_scraper.py  |
| Get data from UniProt               | cw_gui_get_uniprot_data        | python3 cazy_webscarper/gui/gui_get_uniprot_data.py  |
| Get protein sequences from GenBank  | cw_gui_get_genbank_seqs        | python3 cazy_webscarper/gui/gui_genbank_seqs.py  |
| Extract protein sequences from db   | cw_gui_extract_db_seqs         | python3 cazy_webscarper/gui/gui_extract_db_seqs.py  |
| Get structure files from PDB        | cw_gui_get_pdb_structures      | python3 cazy_webscarper/gui/gui_get_pdb_structures.py  |
| Extract information from the db     | cw_gui_query_database          | python3 cazy_webscarper/gui/gui_query_database.py  |
+-------------------------------------+--------------------------------+--------------------------------------------------+

------------------------------
Installing ``gooey`` on Ubuntu
------------------------------

There are well documented issues with specifically installing ``gooey`` and its requirements on Unbuntu (>=v18.0.1), 
when using ``pip``. If you are having issues we recommend the following method as addapted from ...:

1. Install ``gooey`` requirements using ``conda``

.. code-block::bash
    conda install -c anaconda wxpython
    conda install -c conda-forge pygtrie
    conda install -c conda-forge typing-extensions

2. Install ``gooey`` using ``conda``

.. code-block::bash
    conda install -c conda-forge gooey

3. Install ``gooey`` using ``pip``

.. code-block::bash
    pip3 install gooey
