
Retrieving Protein Sequences from GenBank
==========================================

The CAZy webscraper supports the automated retrieval of protein sequences from GenBank for the s
craped CAZymes. This is performed by using the **BioPython** ``Bio.entrez`` module.

.. note::
    For specific information of the ``bio.entrez`` module please see the 
    [entrez documentation](https://biopython.org/docs/1.75/api/Bio.Entrez.html).

.. warning:
    Before using Entrez to access the NCBI’s online resources (via any modules), please read the 
    [NCBI’s Entrez User Requirements](https://www.ncbi.nlm.nih.gov/books/NBK25497/).  
    **If you are found to be abusing NCBI's systems, they can and will ban your access!** 

    In summary:

    For scrapes that would perform more than 100 requests to NCBI, perform the scrape at the 
    **weekend** or **outside USA peak times**.

    Do not perform more than 10 queries to NCBI per second, this is already handled by the 
    webscraper.

    Provide a user email so that NCBI can contact you if there is a problem.


.. note::
    For large requests NCBI suggests using the session history feature of Entrez. At the moment the 
    scraper performs a single request per GenBank accession retrieved from CAZy. Overtime this will 
    may be replaced with using the session history feature to better meet expected practise of the NCBI.


Enabling sequence retrieval
-----------------------------

To enable retrieval of protein sequences from GenBank use the ``-g`` or ``--genbank`` option 
followed by an email address. The provision of the email address is a requirement of Entrez, the 
search system of NCBI.

For example:  
``cazy_webscraper -g dummy_email@domain.com``


Output
------

Phe protein sequence of the scraped CAZymes are retrieved from GenBank are retrieved in 
the FASTA file format and can be written to STDOUT to facilitate piping to a subsequent program or 
written out to disk, within a specified directory.

To write out the FASTA files to a specified directory use the option ``-genbank_output`` followed 
by the path to specified directory.


GenBank synonyms
----------------------------

Often multiple GenBank accession numbers are listed for a given CAZyme within CAZy. However, only the 
first listed accession number is hyperlinked to the GenBank database. Examination of the other listed 
synonyms (referred to as **genbank synonyms** in the webscraper) shows that these GenBank synonyms are 
the result of submission of identical protein sequences, splice site and protein isoforms. It has been 
interpreted that it is the first GenBank accession that is listed and hyperlinked to GenBank is the accession 
number of protein sequence which was used by CAZy to catalogue the CAZyme and the GenBank synonyms were 
identified and listed by having extremely high sequence identity to the catalogued CAZyme. Therefore, the 
webscraper only retrieves the FASTA file of first listed GenBank accession of each CAZyme. 

.. note::
    If you wish to see a feature where the FASTA files of the GenBank synonyms are also retrieved, please 
    raise an issue in the ``cazy_webscraper`` GitHub repository to let us know that this is a feature 
    you would like added.
