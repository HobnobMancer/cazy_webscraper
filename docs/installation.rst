.. _installation:

============
Installation
============

We support three ways to install ``cazy_webscraper``.

* Using `bioconda <https://bioconda.github.io/>`_ **[Recommended]**
* Using `pip/PyPI <https://pypi.python.org/pypi/cazy_webscraper>`_
* Installation from source (`git clone <https://github.com/cazy-project/cazy_webscraper>`_)

------------
Requirements
------------

* A POSIX-compliant operating system, e.g. Linux or MacOS.
* Python 3.8 or later
* An internet connection (to access CAZy and download data)

----------------------------
Installing with ``Bioconda``
----------------------------

.. TIP::
   The most recent stable release of ``cazy_webscraper`` should always be avaiable from the ``bioconda`` channel.

To install ``cazy_webscraper`` using `bioconda <https://bioconda.github.io/>`_ and all required packages, you can use the following command:

.. code-block:: bash

   conda install -c bioconda cazy_webscraper

.. NOTE::
   If ``conda`` is not installed on your system, please see the ``conda`` website for `instructions <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

----------------------------
Installing with ``pip``/PyPI
----------------------------

.. TIP::
   The most recent stable release of ``cazy_webscraper`` should always be avaiable from ``PyPI``.

To install ``cazy_webscraper`` and all required packages using `pip <https://pypi.python.org/pypi/cazy_webscraper>`_, you can use the following command:

.. code-block:: bash

   python -m pip install cazy-webscraper

----------------------
Installing from source
----------------------

``cazy_webscraper`` can be installed from the source code available at the `GitHub repository <https://github.com/cazy-project/cazy_webscraper>`_.

.. WARNING::
   The ``cazy_webscraper`` repository provides the development version of ``cazy_webscraper``. This is the version that is most recently updated, but it may not be the latest stable version. In particular, the development version may contain features that are not yet in a stable version, and it may contain bugs.

.. TIP::
   To obtain the most recent *stable* source code from the ``cazy_webscraper`` repository, download a release from the `releases page <https://github.com/cazy-project/cazy_webscraper/releases>`_ and extract the archive.

.. ATTENTION::
   If you are using `conda`, you can use the ``Makefile`` in the ``cazy_webscraper`` repository to install ``cazy_webscraper`` from source, with the command ``make setup_env``.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Obtaining ``cazy_webscraper`` source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The **development** version of ``cazy_webscraper`` is available at the `GitHub repository <https://github.com/cazy-project/cazy_webscraper>`_ and can be downloaded in the following ways.

To clone the ``cazy_webscraper`` repository using ``git``, you can use the following command:

.. code-block:: bash

   git clone https://github.com/HobnobMancer/cazy_webscraper 

This will download the most recent development version of ``cazy_webscraper`` into a new directory called ``cazy_webscraper``.

Alternatively, the **development** version can be downloaded as an archive file from the link below.. button:: 

* `source code .zip file <https://github.com/HobnobMancer/cazy_webscraper/archive/refs/heads/master.zip>`_

This file can be extracted and will create the directory ``cazy_webscraper``.

^^^^^^^^^^^^^^^^^
Required packages
^^^^^^^^^^^^^^^^^

``cazy_webscraper`` requires several packages to be installed, in order to run. You can install these packages using the ``requirements`` files in the ``cazy_webscraper`` directory. The commands needed depend on your working environment.

******************
Using ``pip``/PyPI
******************

All third-party packages can be installed using ``pip``/PyPI:

.. code-block:: bash

   pip install -r requirements.txt
   pip install -r requirements-dev.txt  # only needed if you are developing the code
   pip install -r requirements-pip.txt  # only needed if you are developing the code

***************
Using ``conda``
***************

The ``conda`` package manager can be used to install all required packages for running ``cazy_webscraper``, but the ``sphinx`` package is not available in ``conda`` and must be installed using ``pip``:

.. code-block:: bash

   conda install --file requirements.txt
   conda install --file requirements-dev.txt  # only needed if you are developing the code
   pip install -r requirements-pip.txt        # only needed if you are developing the code

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Installing ``cazy_webscraper``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``cazy_webscraper`` package can be installed from source using the following command, issued from the root directory of the ``cazy_webscraper`` repository:

.. code-block:: bash

   python setup.py install

If you are intending to edit or develop the code, you can use the ``develop`` option instead of ``install``:

.. code-block:: bash

   pip install -e .
