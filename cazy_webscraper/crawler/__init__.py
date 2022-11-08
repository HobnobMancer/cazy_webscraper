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
"""Retrieve CAZy text file, extract data, apply user scraping filters."""


import logging
import os
import time

from socket import timeout
from tqdm import tqdm
from urllib.error import HTTPError, URLError
from urllib3.exceptions import HTTPError, RequestError
from urllib.request import urlopen
from requests.exceptions import ConnectionError, MissingSchema


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
                    f'Failed to connect to CAZy on try {tries}/{kwargs["max_tries"]}\n'
                    f'Error raised: {err}\n'
                    'Retrying connection to CAZy in 10s'
                )
                time.sleep(10)
            
        if success is False:
            logger.warning(
                f'Failed to connect to CAZy after {kwargs["max_tries"]} tries\n'
                f'Error raised: {err}\n'
            )
            return err
        else:
            return None

    return wrapper


@download_file_decorator
def get_cazy_file(out_path, args, **kwargs):
    """Download plain text file database dumb from the CAZy website

    :param out_path: Path, target path to write out downloaded txt file
    :param args: cmd-line args parser
    :param max_tries: int, max number of times connection to CAZy can be attempted

    Return nothing
    """
    logger = logging.getLogger(__name__)
    download_url = 'http://www.cazy.org/IMG/cazy_data/cazy_data.zip'

    # HTTPError, URLError or timeout error may be raised, handled by wrapper
    response = urlopen(download_url, timeout=args.timeout)

    file_size = int(response.info().get("Content-length"))
    bsize = 1_048_576

    # IOError may be raised, handled by wrapper
    with open(out_path, 'wb') as fh:
        with tqdm(
            total=file_size,
            desc=f"Downloading CAZy txt file",
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

    return
