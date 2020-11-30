==============================
Retrieving Structures from PDB
==============================

The CAZy webscraper supports the automated retrieval of 3D structures from PDB for the scraped CAZymes. 
This is performed by using the ``BioPython`` module ``Bio.PDB``. All protein structures that are listed 
for each CAZyme are retrieved from the [RSCB PDB database](https://www.rcsb.org/).

For specific information of the ``Bio.PDB`` module please see the 
[BioPython documentation](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ).


..warning::
    If many requests are going to be made in a series to PDB (for example a series of 100 
    requests), then it is expected practise to perform the scraper at the **weekend** or 
    **outside USA peak times**.


Enabling structure retrieval
-----------------------------

To enable retrieval of protein structures from PDB then use the option ``-p`` or ``--pdb_output`` 
followed by the desired file format. The options of file format are (as specified in the BioPython 
[documentation](https://biopython.org/docs/1.75/api/Bio.PDB.PDBList.html)):

* mmCif (default, PDBx/mmCif file),
* pdb (format PDB),
* xml (PDBML/XML format),
* mmtf (highly compressed),
* bundle (PDB formatted archive for large structure}


Output
------

The retrieved protein structures can be written to a specified directory by using the ``-pdb_output`` 
option followed by the path to desired output directory. This directory does not already need to be 
existing, the webscraper can make the directory.

The ``BioPython`` ``Bio.PDB`` module does not support writing out the protein structures to standard 
out. Therefore, if ``--pdb_output`` is not given the webscraper will write out the protein structures 
to the output directory specified by the ``-o``/``--output`` option. If the ``-o``/``--output`` option 
is set as standard out (STDOUT), for example if the ``-o``/``--output`` option is not called, then 
the protein structure files will be written to the current working directory.
