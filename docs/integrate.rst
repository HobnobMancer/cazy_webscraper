.. _integrate:

===============================================================
Integrate a local CAZyme database into into downstream analyses
===============================================================

To facilitate integrating the dataset compiled using ``cazy_webscraper``, the data 
downloaded from external databases is compiled into a database using the wildly used local SQLite3 database format.

To facilitate integratting the local CAZyme database into third-party ``Python`` applications, use the ``get_db_connection`` function from ``cazy_webscraper``, which will return an open connection to the CAZyme database from ``sqlalchemy``.

``get_db_connection`` takes 2 required args and one optional arg:  
**Required:**  
- Path to the local CAZyme database (provided as a ``pathlib.Path`` object)
- Bool to set the ``sqlalchemy.create_engine`` param ``sql_echo``. When set to ``True``, the SQL log will be printed to the terminal
**Optional:** 
- Bool to reflect if a new database. Default is ``False``, i.e. connecting to an existing local CAZyme database

Import the function into the ``Python`` script using:

.. code-block:: python

    from cazy_webscraper.sql.sql_orm import get_db_connection
