#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Web scraper to scrape CAZy website and retrieve all protein data.

:cmd_args --config: path to configruration file
:cmd_args --classes: specify CAZy classes to scrape
:cmd_args --database: provide path to a local SQL database to add additional data to
:cmd_args --families: specify CAZy families to retrieve CAZymes from
:cmd_args --force: force overwriting content in exisiting output directory
:cmd_args --genera: specify Genera to retrieve CAZymes from
:cmd_args --kingdoms: specify taxonomy Kingdoms to scrape proteins from
:cmd_args --log: path to log file, enables writing out log messages to a log file
:cmd_args --nodelete: if true does not delete content in pre-existing output directory
:cmd_args --output: path to output directory
:cmd_args --retries: specify the number of times to try scraping a page if connection fails
:cmd_args --subfamilies: enable retrieval of subfamilies from CAZy
:cmd_args --species: specify species to retrieve CAZymes from
:cmd_args --strains: specify specific strains of species to retrieve CAZymes from
:cmd_args --timeout: specify the maximum time (in seconds) before determining connection timed out
:cmd_args --verbose: change logger level from warning to info, verbose logging
"""


import json
import logging
import os
import sys

import pandas as pd

from datetime import datetime
from pathlib import Path
from typing import List, Optional

from tqdm import tqdm

from scraper import crawler
from scraper.sql import sql_orm, sql_interface
from scraper.utilities import (
    build_logger,
    config_logger,
    file_io,
    parsers,
    parse_configuration,
    termcolour,
)


# Define constants

__version__ = "1.0.0-beta"

VERSION_INFO = [
    termcolour(
        f"cazy_webscraper version: {__version__}",
        "cyan",
    ),
]

CITATION_INFO = [
    termcolour(
        "If you use cazy_webscraper in your work, please cite the following publication:",
        "green",
    ),
    termcolour(
        "\tHobbs, E. E. M., Pritchard, L., Chapman, S., Gloster, T. M.,",
        "yellow",
    ),
    termcolour(
        "\t(2021) cazy_webscraper Microbiology Society Annual Conference 2021 poster. ",
        "yellow",
    ),
    termcolour(
        "\tFigShare. Poster.",
        "yellow",
    ),
    termcolour(
        "\thttps://doi.org/10.6084/m9.figshare.14370860.v7",
        "yellow",
    ),
]


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, logger and coordinate overal scrapping of CAZy.

    The collected data can be stored as a singel dataframe containing (not split), split into
    separate dataframes by class or by family. Excluded classes are CAZy classes not specified in
    the configuration file and thus, will not be scraped. User_cazy_families is the list of CAZy
    families specified to be scraped in the configration file.
    """
    cazy_home_url = "http://www.cazy.org"

    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # used in naming files
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # used in terminating message
    start_time = pd.to_datetime(start_time)

    # Program preparation
    if argv is None:
        parser = parsers.cazy_webscraper_parser.build_parser()
        args = parser.parse_args()
    else:
        parser = parsers.cazy_webscraper_parser.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)
    
    # check if printing out version or citation information
    if args.version:
        sys.stderr.write("\n".join(VERSION_INFO) + "\n")
        return
    
    if args.citation:
        sys.stderr.write("\n".join(CITATION_INFO) + "\n")
        return

    # check correct output was provided, exit if not operable
    check_user_input(args)

    logger.info("Parsing configuration")
    (
        excluded_classes,
        config_dict,
        cazy_class_synonym_dict,
        kingdoms_filters,
        taxonomy_filters,
    ) = parse_configuration.parse_configuration(args)

    scrape_config_message = (
        "Configuration:\n"
        f"Classes to scrape: {config_dict['classes']}\n"
        f"GH fams to scrape: {config_dict['Glycoside Hydrolases (GHs)']}\n"
        f"GT fams to scrape: {config_dict['GlycosylTransferases (GTs)']}\n"
        f"PL fams to scrape: {config_dict['Polysaccharide Lyases (PLs)']}\n"
        f"CE fams to scrape: {config_dict['Carbohydrate Esterases (CEs)']}\n"
        f"AA fams to scrape: {config_dict['Auxiliary Activities (AAs)']}\n"
        f"CBM fams to scrape: {config_dict['Carbohydrate-Binding Modules (CBMs)']}\n"
        f"Scraping subfamilies: {args.subfamilies}"
    )

    if len(taxonomy_filters) != 0:
        scrape_config_message += "\nTaxonomy filters applied."
    
    if len(kingdoms_filters) < 5:
        scrape_config_message += f"\nScraping only tax kingdoms: {kingdoms_filters}"

    logger.info(termcolour(scrape_config_message, "cyan"))

    if args.database:  # adding data to an EXISTING database
        connection, logger_name, cache_dir = connect_existing_db(args, time_stamp)
    
    else:  # build a new database
        connection, logger_name, cache_dir = connect_to_new_db(args, time_stamp)

    logger.info("Adding log of scrape to the local CAZyme database")
    sql_interface.log_scrape_in_db(
        time_stamp,
        config_dict,
        kingdoms_filters,
        taxonomy_filters,
        connection,
        args,
    )

    if args.dict_output is not None:
        logger.info("Building dict of CAZy family annotations")
        logger.info(f"Output directory: {(args.output).parent}")
        logger.info(f"Force overwriting exiting output file: {args.force}")

        if str((args.dict_output).parent) != '.':
            # dirs defined in output put
            file_io.make_target_directory(args.dict_output, args.force)

        cazy_dict = {}  # {protein_accession: [CAZy fams]}

    else:
        cazy_dict = None
    
    if args.log is not None:  # write additional log files to user specified dir
        logger_name = args.log.split(".")[0]
    else:
        # write the additional log files to the .cazy_webscraper/log dir
        logger_dir = Path(f"{str(cache_dir.parent)}/logs")
        file_io.make_output_directory(logger_dir, args.force, args.nodelete)
        # add logger dir path to the logger name
        logger_name = f"{logger_dir}/{str(Path(logger_name).name)}"
        
    # create dir to cache downloaded text files and logs of failed scrapes, connections and data errs
    file_io.make_output_directory(cache_dir, args.force, args.nodelete)

    logger.info(f"Created cache dir: {cache_dir}")

    logger.info("Starting retrieval of data from CAZy")
    get_cazy_data(
        cazy_home_url,
        excluded_classes,
        config_dict,
        cazy_class_synonym_dict,
        taxonomy_filters,
        connection,
        cazy_dict,
        cache_dir,
        args,
        logger_name,
    )

    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = pd.to_datetime(end_time)
    total_time = end_time - start_time

    logger.info(
        "Finished scraping CAZy. Terminating program.\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}"
        f"Version: {VERSION_INFO}\n"
        f"Citation: {CITATION_INFO}"
    )

    print(
        "=====================cazy_webscraper=====================\n"
        "Finished scraping CAZy\n"
        f"Scrape initated at {start_time}\n"
        f"Scrape finished at {end_time}\n"
        f"Total run time: {total_time}"
        f"Version: {VERSION_INFO}\n"
        f"Citation: {CITATION_INFO}"
    )


def get_cazy_data(
    cazy_home_url,
    excluded_classes,
    config_dict,
    cazy_class_synonym_dict,
    taxonomy_filters,
    connection,
    cazy_dict,
    cache_dir,
    args,
    logger_name,
):
    """Coordinate retrieval of data from the CAZy website.

    This function coordinates the crawling through the CAZy website by calling the appropriate
    functions, and then retrieving the protein data by calling to the appropriate data again.

    :param cazy_home_url: str, url of CAZy home page
    :param excluded_classes: list, list of classes to not scrape from CAZy
    :param config_dict: dict, user defined configuration of the scraper
    :param cazy_class_synonym_dict: dict, dictionary of excepct CAZy synonyms for CAZy classes
    :param taxonomy_filters: set of genera, species and strains to restrict the scrape to
    :param connection: connection, open database connection (if opted to produce a CAZyme database)
    :param cazy_dict: dict to add store protein accession and fam annotations (if args.dict_output is True)
    :param cache_dir: Path to dir to write out downloaded family txt files
    :param args: cmd args parser
    :param logger_name: str, name used for additional logger files

    Return nothing.
    """
    # define paths for additional logs files
    connection_failures_logger = build_logger(Path(f"{logger_name}_connection_failures.log"))
    sql_failures_logger = build_logger(Path(f"{logger_name}_SQL_errors.log"))
    format_failures_logger = build_logger(Path(f"{logger_name}_format_and_parsing_errors.log"))

    # retrieve links to CAZy class pages, return list of CazyClass instances
    cazy_classes = crawler.get_cazy_classes(
        cazy_home_url,
        excluded_classes,
        cazy_class_synonym_dict,
        args,
    )

    # scrape each retrieved class page
    for cazy_class in tqdm(cazy_classes, desc="Parsing CAZy classes"):

        # first attempt of scraping, retrieve URLs to CAZy families
        if len(list(cazy_class.failed_families.keys())) == 0:

            # retrieve Family class instances (in class_families), containing the Family page URL
            class_families, error_message, incorrect_urls = crawler.get_cazy_family_urls(
                cazy_class.name,
                cazy_class.url,
                cazy_home_url,
                cache_dir,
                args,
            )

            if incorrect_urls is not None:  # log families for which compiled URL is incorrect
                [connection_failures_logger.warning(url) for url in incorrect_urls]

            if class_families is None:  # couldn't retrieve URLs to families for working CAZy class
                cazy_class.tries += 1  # add one to the number of scrapping attempts

                # check if maximum number of attempts to connect have been met

                if cazy_class.tries == (args.retries + 1):  # Maximum number of tries met
                    connection_failures_logger.warning(
                        f"{cazy_class.url}\t"
                        f"{cazy_class.name}\t"
                        "No CAZy familes from this class were scraped\t"
                        f"{error_message}"
                    )

                else:
                    cazy_classes.append(cazy_class)  # retry scraping Class page later Classes

                continue  # go onto next CAZy class

        # Not first try, scrape only the families for which a connections to CAZy previously failed
        else:
            class_families = list(cazy_class.failed_families.keys())
        
        continue

        # Scrape the familes of the current CAZy class, retrieving protein data
        for family in tqdm(class_families, desc=f"Parsing {cazy_class.name} families"):

            if (config_dict is not None) or (config_dict[cazy_class.name] is not None):
            # Specific (sub)families were specified for scraping

                # Allow retrieval of subfamilies when only the parent CAZy family was named in the
                # config file, by searching by the family not subfamily in the config file
                if (args.subfamilies is True) and (family.name.find("_") != -1):
                    name_check = family.name[: (family.name.find("_"))]
                else:
                    name_check = family.name

                if name_check not in config_dict[cazy_class.name]:
                    continue
            
            # else:
            # No (sub)families were specified
            # Don't perform a name check and scrape all families of the CAZy class

            (
                family,
                failed_url_connections,
                family_sql_failures,
                format_errors,
                connection,
                cazy_dict,
            ) = crawler.parse_family(
                family,
                cazy_home_url,
                taxonomy_filters,
                args,
                connection,
                cazy_dict,
                args,
            )

            # check if family file was successfully retrieved
            if family.path.exists() is False:
                # file not successfully 
                try:
                    cazy_class.failed_families[family] += 1
                except KeyError:
                    cazy_class.failed_families[family] = 1  # first attempt
                
                # check if max number of attempts to connect family pages has been met
                if cazy_class.failed_families[family] == (args.retries + 1):
                    del cazy_class.failed_families[family]  # do not try another scrape
                    continue

                if len(list(cazy_class.failed_families.keys())) != 0:
                    # if there are families with previously failed connection attempts
                    # and remaining tries, retry connection after working through other classes
                    cazy_classes.append(cazy_class)

            # write out errors to their respective log files
            for error in failed_url_connections:
                connection_failures_logger.warning(error)

            for error in family_sql_failures:
                sql_failures_logger.warning(error)

            for error in format_errors:
                format_failures_logger.warning(error)

    if cazy_dict is not None:
        with open(args.dict_output, 'w') as f:
            json.dump(connection, f)

    return


def check_user_input(args):
    """Check cmd-line args are suitable.
    
    :param args: cmd-line args parser.
    
    Return nothing."""
    logger = logging.getLogger(__name__)

    if args.database and args.db_output:
        warning_message = (
            "Target path for a NEW database (--db_output, -d) and\n"
            "a path to an EXISTING database (--database, -D) were provided."
            "Please provide one OR the other.\n"
            "Terminating program."
        )
        logger.warning(termcolour(warning_message, "red")
        )
        sys.exit(1)
    
    if (args.db_output is None) and (args.database is None) and (args.dict_output is None) and (not args.no_db):
        warning_message = (
            "No target path for an output database out JSON file provided.\n"
            "Terminating program."
        )
        logger.warning(termcolour(warning_message, "red")
        )
        sys.exit(1)

    if args.no_db and (args.dict_output is None):
        warning_message = (
            "Opted to not write a database and no target path for a JSON file provided.\n"
            "Please provide at least one method for writing out the retrieved data.\n"
            "Terminating program."
        )
        logger.warning(termcolour(warning_message, "red")
        )
        sys.exit(1)

    return


def connect_existing_db(args, time_stamp):
    """Coordinate connecting to an existing local CAZyme database, define logger name and cache dir
    
    :param args: cmd-line args parser
    :param time_stamp: str, time cazy_webscraper was invoked

    Return connection to local CAZyme database, logger file name, and path to cache dir
    """
    logger = logging.getLogger(__name__)
    
    logger.info("Adding data to an existing local CAZyme database")

    if os.path.isfile(args.database) is False:
        logger.error(
            "Could not find local CAZy database.\n"
            "Check path is correct.\n"
            "Terminating programme."
        )
        sys.exit(1)

    try:
        connection = sql_orm.get_db_connection(args)
        logger.info("Opened connection to local CAZyme database")
    except Exception:
        logger.error(
            "Failed to open connection to an exiting local CAZyme database\n."
            "Terminating program\n",
            exc_info=True,
        )
        sys.exit(1)
    
    # used for naming additional log files
    logger_name = str(args.database).split('.')[0]

    # define path to cache family txt files
    cache_dir = Path(f"{str(args.database.parent)}/.cazy_webscraper_{time_stamp}/cache")

    return connection, logger_name, cache_dir
    

def connect_to_new_db(args, time_stamp):
    """Build and connect to a new local CAZyme database.
    
    :param args: cmd-line args parser
    :param time_stamp: str, time cazy_Webscraper was invoked
    
    Return connection to the database, name of the logger, and path to the cache dir
    """
    logger = logging.getLogger(__name__)

    if args.db_output is not None:  # user defined target output for the NEW database
        
        if (os.path.isfile(args.db_output)):  # target file exists
            if args.force:
                logger.warning(
                    "Overwriting existing local CAZyme database at:\n"
                    f"{args.db_output}"
                )

            else:
                logger.warning(
                    "Target path for new database already exists.\n"
                    "Either enable forced overwriting (-f) or add data this data (-D).\n"
                    "Terminating program."
                )
                sys.exit(1)
        
        else:  # may need to build dirs
            logger.info(
                "Building new local CAZyme database\n"
                f"Output directory: {(args.output).parent}\n"
                f"Force overwriting exiting output file: {args.force}"
            )

        if str((args.db_output).parent) != '.':  # dirs defined in output put
            file_io.make_target_directory(args.db_output, args.force)
            cache_dir = Path(f"{str(args.db_output.parent)}/.cazy_webscraper_{time_stamp}/cache")
            
        else:  # writing to cwd
            cache_dir = Path(f".cazy_webscraper_{time_stamp}/cache")

        logger_name = args.db_output.split('.')[0]
        db_path = args.db_output
    
    else:
        logger.info("Using default database name and writing to cwd")
        db_path = Path(f"cazy_webscraper_{time_stamp}")
        cache_dir = Path(f".cazy_webscraper_{time_stamp}/cache")
        logger_name = f'cazy_webscraper_{time_stamp}'
    
    try:
        connection = sql_orm.build_db(db_path, time_stamp)
        logger.info(f"Built new local CAZyme database at\n{db_path}")
    except Exception:
        logger.error(
            "Failed to build new SQL database\n."
            "Terminating program",
            exc_info=True,
        )
        sys.exit(1)

    return connection, logger_name, cache_dir


if __name__ == "__main__":
    main()


