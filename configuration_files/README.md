# Configuration files

To faciliate the reproducbility of compiling a local CAZyme database, `cazy_webscraper` can be configured using a YAML file.

This directory contains:
- [`template_cazy_config.yaml`]() - an empty template cofig file for configuring scraping CAZy
- [`template_get)data_config.yaml`]() - an empty template config file for configuring the retrieval of data from UniProt, GenBank and PDB
- [`example_cazy_config.yaml`]() - An example configuration file for scraping CAZy, to scrape all proteins from Bacteria and Viruses in GH and CE, PL families PL1 and PL2
- [`example_get_data_config.yaml`]() - An example configuration file for retrieving data from UniProt, GenBank or PDB for proteins from genera _Trichoderma_ and _Aspergillus_, and the species _Pythium ultimum_, and belonging for GT families GT1, GT5 and GT6

## Using configuration files

Whenever cazy_webscraper is invoked and adds data to a database, the configuration of cazy_webscraper (this is the kingdoms, genera, species, strains, CAZy classes and CAZy family filters which were applied) and the data and time the scrape was initiated is logged in the database. However, for optimal reproduction of how cazy_webscraper was used in your research, you can create shareable documentation that others can use to reproduce your CAZy dateset. This is achieved by creating a configuration file rather than configuring the performance of cazy_webscraper at the command line.

### Creating a configuration file
cazy_webscraper uses the YAML file type for its configuraiton file; if you are new to YAML files please find more detailed information on YAML files [here](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html).

A template and example configuration file for scrapping CAZy using cazy_webscraper can be found in the repo, in the configuration_files directory.

The configuration YAML must contain the following tags/headings (identical to how they are presented below):

- classes
- Glycoside Hydrolases (GHs)
- GlycosylTransferases (GTs)
- Polysaccharide Lyases (PLs)
- Carbohydrate Esterases (CEs)
- Auxiliary Activities (AAs)
- Carbohydrate-Binding Modules (CBMs)
- genera
- species
- strains
- kingdoms

**Note**
The order of the tags/headings does not matter.

### Scraping specific CAZy classes
Under the classes heading list any classes to be scrapped. For each CAZy class listed under ‘classes’, CAZymes will be retrieved for every CAZy family within the CAZy class.

Each class must be listed on a separate line, indented by 4 spaces, and the class name encapsulated with single or double quotation marks. For example:

```yaml
classes:
    - "GH"
    - "PL"
```
The same CAZy class name synonyms used for the command line are accepted for the configuration file.

Scraping specific CAZy families
Under the each of the class names listed in the configuration file, list the names of specific families to be scraped from that class. The respective classes of the specificed families do not need to be added to the ‘classes’ list.

Write the true name of the family not only it’s number, for example GH1 is excepted by 1 is not.

Name families using the standard CAZy nomenclature, such as “GT2” and NOT “GlycosylTransferases_2”. Additionally, use the standard CAZy notation for subfamilies (GH3_1).

**Warning**
If any subfamilies are listed within the configuration file, the retrieval of subfamilies must be enabled at the command line uisng --subfamilies.

Each family must be listed on a separate line and the name surrounded by double or single quotation marks. For example:

```yaml
Glycoside Hydrolases (GHs):
    - "GH1"
    - "GH2"
    - "GH3_1"
```

### Example configuration file
Below is an example of the content you may wish to put in a configuration file. Using this file will retrieve all CAZymes in CAZy class AA, CAZy families GH1, GH3 and PL9 that are either derived from a bacterial or Trichoderma species.

```yaml
classes:
   - "AA"
Glycoside Hydrolases (GHs):
   - "GH1"
   - "GH3"
GlycosylTransferases (GTs):
Polysaccharide Lyases (PLs):
   - "PL9"
Carbohydrate Esterases (CEs):
Auxiliary Activities (AAs):
Carbohydrate-Binding Modules (CBMs):
genera:
   - "Trichoderma"
species:
strains:
kingdoms:
   - "Bacteria"
```

**Note**
Indentations consist of 4 spaces.

You can add ‘comments’ to configuration file. Comments are section of text that are not read by cazy_webscraper and allow you to add notes to your configuration file. For example:

### Using a configuration file
Once you have created a configuration file (we recommend modifying the template one provided with cazy_webscraper you then need to invoke cazy_webscraper and tell it you are using a configuration file. To do this we add the --config/-c flag to the cazy_webscraper command, followed by the path to the configuration file.

The path we pass to cazy_webscraper is a relative path. This means cazy_webscraper will start in the directory the terminal is currently pointed out, and follow the path from there. For example, if we used the command:

cazy_webscraper -c scraper/scraper_config.yaml
Then the computer will look for a directory called scraper in the directory the terminal is looking at, then look within the scraper directory for a yaml file called scraper_config.yaml.

**Note**
To check which directory cazy_webscraper is pointed at type pwd into the terminal and hit enter. This is the ‘Present Working Directory’ command, which will print the path to the directory the terminal is presently looking at.

**Warning**
Your path must point directly to the YAML file. Don’t forget the ‘.yaml’ file extension!
