#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Download the CAZy database dump (txt file)."""


import argparse
import logging
import os
import time

from pathlib import Path
from socket import timeout

from tqdm import tqdm
from urllib.error import HTTPError, URLError
from urllib3.exceptions import HTTPError, RequestError
from urllib.request import urlopen
from requests.exceptions import ConnectionError, MissingSchema

from cazy_webscraper import DOWNLOAD_URL


logger = logging.getLogger(__name__)


def get_cazy_db_dump(cache_dir: Path, time_stamp: str, args: argparse.ArgumentParser):
    """Retrieve txt file of CAZy db dump from CAZy or the local disk.

    :param cache_dir: Path(), path to directory where cache is written to
    :param time_stamp: str, date and time cazy_webscraper was intiated
    :param args: cmd-line args parser

    Returns the path to the CAZy DB dump (txt file)
    """
    if args.cazy_data:   # retrieve lines from predownloaded CAZy txt file
        return args.cazy_data

    # download cazy db dump
    cazy_txt_path = cache_dir / f"cazy_db_{time_stamp}.zip"
    tries, retries, success, err_message = 0, (args.retries + 1), False, None

    err_message = None
    while (tries <= retries) and (not success):
        err_message = download_cazy(cazy_txt_path, args, max_tries=(args.retries + 1))

        if err_message is None:
            break

        else:
            tries += 1

    return cazy_txt_path


def download_file_decorator(func):
    """Decorator to re-invoke the wrapped function up to 'args.retries' times."""

    def wrapper(*args, **kwargs):
        logger = logging.getLogger(__name__)
        tries, success, err = 0, False, None

        while not success and (tries < kwargs['max_tries']):
            # reset storing error messsage
            err_message = None

            try:
                func(*args, **kwargs)

            except (
                IOError,
                HTTPError,
                URLError,
                timeout,
                ConnectionError,
                OSError,
                MissingSchema,
                RequestError,
            ) as err_message:
                success = False
                err = err_message

            if err is None:
                success = True

            tries += 1

            if (not success) and (tries < kwargs['max_tries']):
                logger.warning(
                    'Failed to connect to CAZy on try %s/%s\n'
                    'Error raised: %s\n'
                    'Retrying connection to CAZy in 10s',
                    tries, kwargs["max_tries"], err
                )
                time.sleep(10)

        if success is False:
            logger.warning(
                'Failed to connect to CAZy after %s tries\n'
                'Error raised: %s\n',
                kwargs["max_tries"], err
            )
            return err
        else:
            return None

    return wrapper


@download_file_decorator
def download_cazy(out_path: Path, args: argparse.ArgumentParser, **kwargs):
    """Download plain text file database dumb from the CAZy website

    :param out_path: Path, target path to write out downloaded txt file
    :param args: cmd-line args parser
    :param max_tries: int, max number of times connection to CAZy can be attempted
    """
    logger = logging.getLogger(__name__)

    # HTTPError, URLError or timeout error may be raised, handled by wrapper
    with urlopen(DOWNLOAD_URL, timeout=args.timeout) as response:
        file_size = int(response.info().get("Content-length"))
        bsize = 1_048_576

        # IOError may be raised, handled by wrapper
        with open(out_path, 'wb') as fh:
            with tqdm(
                total=file_size,
                desc="Downloading CAZy txt file",
            ) as pbar:
                while True:
                    buffer = response.read(bsize)
                    if not buffer:
                        break
                    pbar.update(len(buffer))
                    fh.write(buffer)

        if os.path.isfile(out_path) is False:
            logger.error('CAZy txt file not created locally.')
            raise IOError
