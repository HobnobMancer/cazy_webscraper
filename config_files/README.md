# Configuration files

To faciliate the reproducbility of compiling a local CAZyme database, `cazy_webscraper` can be configured using a YAML file.

This directory contains:
- [`template_cazy_config.yaml`]() - an empty template cofig file for configuring scraping CAZy
- [`template_get)data_config.yaml`]() - an empty template config file for configuring the retrieval of data from UniProt, GenBank and PDB
- [`example_cazy_config.yaml`]() - An example configuration file for scraping CAZy, to scrape all proteins from Bacteria and Viruses in GH and CE, PL families PL1 and PL2
- [`example_get_data_config.yaml`]() - An example configuration file for retrieving data from UniProt, GenBank or PDB for proteins from genera _Trichoderma_ and _Aspergillus_, and the species _Pythium ultimum_, and belonging for GT families GT1, GT5 and GT6
