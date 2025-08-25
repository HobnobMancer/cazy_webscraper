========================================================
Quick Start: Configuring ``cazy_webscraper`` to Scrape CAZy
========================================================

``cazy_webscraper`` retrieves user-specified data sets from CAZy and related databases. Configuration is possible via the **command line** or a **YAML configuration file**.

-------------------------------
Basic Usage and Output Settings
-------------------------------

The only required argument is your email address:

.. code-block:: bash

   cazy_webscraper download myemail@domain.com

By default, this scrapes all of CAZy and writes the database to the current directory.

To specify an output file:

.. code-block:: bash

   cazy_webscraper download myemail@domain.com --output cazyme_database/cazyme_database.db

To overwrite an existing database, add ``--force`` (or ``-f``):

.. code-block:: bash

   cazy_webscraper -o data/cazyme_database.db -f

To prevent deleting existing files in the output directory, add ``--nodelete`` (or ``-n``):

.. code-block:: bash

   cazy_webscraper -o data/cazyme_database.db -f -n

-------------------------------
Filtering by CAZy Classes/Families
-------------------------------

Scrape specific CAZy classes:

.. code-block:: bash

   cazy_webscraper --classes GH,CE

This retrieves all CAZymes from Glycoside Hydrolases and Carbohydrate Esterases.

Scrape specific families:

.. code-block:: bash

   cazy_webscraper --families GH2,PL5,CE1

You can combine class and family filters:

.. code-block:: bash

   cazy_webscraper --classes CE --families GH1,PL9

-------------------------------
Taxonomic Filters
-------------------------------

Scrape by kingdom:

.. code-block:: bash

   cazy_webscraper --kingdoms bacteria,eukaryota

Scrape by genus:

.. code-block:: bash

   cazy_webscraper --genera Aspergillus,Trichoderma

Scrape by species (use quotes if there are spaces):

.. code-block:: bash

   cazy_webscraper --species "Aspergillus niger,Trichoderma reesei"

Scrape by strain:

.. code-block:: bash

   cazy_webscraper --strains "Aspergillus niger ATCC 1015"

You can combine these filters as needed. For example, to retrieve all CAZymes from viral species, Aspergillus genus, and specific Layia species:

.. code-block:: bash

   cazy_webscraper --kingdoms viruses --genera Aspergillus --species "Layia carnosa,Layia chrysanthemoides"

-------------------------------
Retrieving Subfamily Annotations
-------------------------------

To include CAZy subfamily annotations:

.. code-block:: bash

   cazy_webscraper --families GH3 --subfamilies

-------------------------------
Using a Configuration File
-------------------------------

You can specify filters in a YAML configuration file. Example:

.. code-block:: yaml

   classes:
      - "AA"
   Glycoside Hydrolases (GHs):
      - "GH1"
      - "GH3"
   Polysaccharide Lyases (PLs):
      - "PL9"
   genera:
      - "Trichoderma"
   kingdoms:
      - "Bacteria"

Invoke with:

.. code-block:: bash

   cazy_webscraper -c path/to/config.yaml

-------------------------------
Other Useful Options
-------------------------------

- Use ``--cazy_data <file>`` to scrape from a downloaded CAZy text file.
- Use ``--log <file>`` to write terminal output to a log file.
- Use ``--verbose`` for detailed logging.
- Use ``--timeout <seconds>`` to set the connection timeout.

-------------------------------
Further Help
-------------------------------

See the full documentation or run ``cazy_webscraper --help`` for more options.
