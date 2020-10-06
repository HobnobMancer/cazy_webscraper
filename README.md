# cazy_webscrapper
The `cazy_webscraper` is a Python3 package for the automated retrieval of protein data from the [CAZy](http://wwww.cazy.org/) database. This program is free to use under the MIT license when proper recognition is given.

The CAZy database is one of the largest databases for cataloging Carbohydrate Active enZymes, and is one of the most frequently visited protein data bases for researchers in this field. However, CAZy provides not method of automated retrieval of data, especially for large dataset requests. The `cazy_webscraper` provides a way for the automated retrieval of large datasets from the CAZy database. The retrieved protein data is written out into dataframes, and protein sequences written out to FASTA files.

Although it has not been finalised, the way protein data is grouped into the output files will be determined by the user selecting one of three options:
1. One output: All data retrieved from CAZy is written out to a singel dataframe, and all protein sequences written out to one FASTA file
2. One per class: A dataframe and FASTA file is written out for each CAZy class, for example one dataframe and FASTA file is written out for glycoside hydrolases and another written out for polysaccharide lyases.
3. One per family: 

**Note:** This webscraper is still in early development and is not fully opperational yet. Furthermore, additional features may be added later, such as ability to configure the webscraper to only scrape the data for a single CAZy class instead of scraping the entire database.

For more information on planning/development please see the Wiki.
