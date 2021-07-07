==============================
``cazy_webscraper`` Quickstart
==============================

------------
Installation
------------

The most recent version of ``cazy_webscraper`` can be installed on your local machine using ``conda`` or ``pip``. Both methods will install the ``cazy_webscraper`` command-line tool, and the Python package ``cazy_webscraper``.

^^^^^^^^^
``conda``
^^^^^^^^^

.. code-block:: bash

   conda install cazy_webscraper

^^^^^^^
``pip``
^^^^^^^

``pip`` should distribute the latest version of ``cazy_webscraper``, although there may be some minor lag between GitHub releases and ``pip``.

.. code-block:: bash

   pip3 install cazy_webscraper

.. TIP::
    ``cazy_webscraper`` can also be installed directly from source. More detailed, and alternative installation instructions can be found in the :ref:`installation` section.


----------------------
Getting Started Poster
----------------------

For a quick summary of how to get started, check out our poster:

    Hobbs, Emma E. M.; Pritchard, Leighton; Gloster, Tracey M.; Chapman, Sean (2021): cazy_webscraper - getting started. FigShare. Poster. `https://doi.org/10.6084/m9.figshare.14370869.v3 <https://doi.org/10.6084/m9.figshare.14370869.v3>`_ 

-------------
Default Usage
-------------

To invoke the webscraper with its default functionality call the webscraper at the command line:  

.. code-block:: bash

  cazy_webscraper 

The default behaviour of the scraper is:

* Scrape all entries in the CAZy database
* Do not split/separate the data: produce a single output
* Write the resulting data to STDOUT
* Do not retrieve subfamilies (subfamily members will be retrieved but only their parent family be listed)
* Do not retrieve FASTA files from GenBank
* Do not retrieve protein sequences from PDB