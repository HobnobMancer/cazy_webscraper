===================================
The Local CAZyme Database Structure
===================================

To facilitate the thorough interrogation of data retrieved from CAZy and minimise storing duplicate and redundant data, data retrieved from CAZy is stored in a local SQL database. 
Every CAZyme scraped from CAZy has the following data:

* Protein name
* CAZy (sub)family
* GenBank accession(s)

Each CAZyme may or may not have the following data, depending on the entry:

* EC number(s)
* UniProt acession(s)
* PDB accession(s)

.. NOTE::
    EC numbers, UniProt accessions and PDB accessions can be retrieved from UniProt for CAZymes 
    in the local CAZyme database using ``cazy_webscraper``.

---------------
Database Schema
---------------

The database ORM (schema) can be viewed `Here <https://hobnobmancer.github.io/cazy_webscraper/database_schema.pdf>`_

--------------
The Logs table
--------------

The database built by ``cazy_webscraper`` contains a table called 'Logs'. This table logs every 
scrape of CAZy, UniProt and GenBank which added data to the database.

The table contains the following columns and data:

* **log_id:** Autoincrement ID number
* **date:** Date scrape was initated (in ISO format)
* **time:** Time scrape was initated (in ISO format)
* **database:** Name of the external database from which data was retrieved (i.e. 'CAZy', 'UniProt' or 'GenBank')
* **retrieved_annotations:** List of annotation types retrieved (e.g. 'EC number, PDB accession, Sequence')
* **classes:** CAZy classes for which data was retrieved
* **families:** CAZy families for which data was retrievedclasses
* **kingdoms:** taxonomy Kingdoms filteres applied
* **genera_filter:** taxonomy Kingdoms filteres applied
* **species_filter:** taxonomy Kingdoms filteres applied
* **strains_filter:** taxonomy Kingdoms filteres applied
* **ec_filter:** EC fliters applied to retrieve data for CAZymes annotated with the specified EC numbers (only applies to retrieval of data from UniProt, GenBank and PDB)
* **cmd_line:** Reproduction of the command line arguments passed to ``cazy_webscraper``.

The 'Logs' table allows any one who uses the database to see how the dataset was compiled.

--------
Genbanks
--------

The Genbanks table contains data retrieved from CAZy and NCBI GenBank. From CAZy the GenBank protein accession. 
The protein sequence for the protein can be retrieved from NCBI GenBank using ``cazy_webscraper``, and stored 
in the Genbanks table, as well as the date the sequence was updated in NCBI. The sequence update date is used 
to check if a newly retrieved protein sequence from GenBank is newer than the sequence stored in the database, and the 
sequence in the database should be updated.

---------------------
Genbanks_CazyFamilies
---------------------

The Genbanks_CazyFamilies table is a relationship. The table defines which protein is 
assigned to which CAZy family.

------------
CazyFamilies
------------

The CazyFamilies lists of CAZy families retrieved from CAZy. If CAZy subfamilies are retrieved 
each CAZy subfamily is associated with its parent CAZy family.

----
Taxs
----

The Taxs table stores taxonomy database, storing the genus and species of the source organisms of CAZymes 
retrieved from CAZy. Each source organism is associated with a taxonomic class.

--------
Kingdoms
--------

The Kingdoms table lists all taxonomic kingdoms of the source organisms downloaded from CAZy.

--------
UniProts
--------

The UniProts table contains protein data retrieved from UniProt using ``cazy_webscraper``. This includes: 

* UniProt ID
* Protein name
* Protein sequence
* Date the protein sequence was last upated in UniProt

------------
Genbanks_Ecs
------------

The Genbanks_Ecs table is a relationship table, and defines which proteins are annotated with which EC numbers. 

---
Ecs
---

The Ecs table lists all EC numbers retrieved from UniProt using ``cazy_webscraper``. The EC numbers stored in the 
local CAZome database are do not have the 'EC' prefix.

-------------
Genbanks_Pdbs
-------------

The Genbanks_Pdbs is a relationship table, and defines which PDB accessions belong to which protein.

----
Pdbs
----

The Pdbs table contains all PDB accessions retrieved from UniProt using ``cazy_webscraper``.

.. NOTE::
    Not all PDB accessions represented in a CAZyme record at CAZy are necessarily present in PDB. For example, some accessions are placeholders while structures are under embargo.

.. NOTE::
    PDB/RCSB protein structures are not recorded in the local SQLite3 database. They are written to disk in a user-specified directory.
