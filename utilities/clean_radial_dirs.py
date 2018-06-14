#!/usr/bin/env python
"""
@file clean_radial_dirs.py
@author Mike Smith
@email michaesm@marine.rutgers.edu
@brief Clean up radial directories
@purpose In a generic ~/data/radials/ directory, this script recursively searches for site subfolders and moves radial
files into the format ~/data/radials/site/yyyy-mm/. Script also leaves last 30 days in the main radial site folder
"""
import logging
import os

import shutil
import sys
from codar_processing.common import create_dir, timestamp_from_lluv_filename
from datetime import timedelta, datetime

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

source_dir = '/home/codaradm/data/radials/'  # main directory that contains radial files
types = ('.ruv', '.euv')  # File types we want to move
base_date = datetime.utcnow() - timedelta(days=30)

# list the subdirectories from source_dir. Only add directories
site_dirs = [os.path.join(source_dir, o) for o in os.listdir(source_dir) if os.path.isdir(os.path.join(source_dir, o))]

# iterate through the site directories
for site in site_dirs:

    # grab a list of all files in each site directory that has all extension in the tuple, 'types', created above
    files = sorted([f.path for f in os.scandir(site) if f.name.endswith(types)])

    # if files_grabbed is not empty
    if files:
        # iterate through each file in the list files_grabbed
        for filename in files:
            basename = os.path.basename(filename)
            # grab the filename from the path/filename combo, filename
            timestamp = timestamp_from_lluv_filename(basename)

            # create directory name in the format yyyy_mm
            new_dir = os.path.join(site, timestamp.strftime('%Y_%m'))

            # try to create the directory.
            create_dir(new_dir)

            if base_date > timestamp:  # if file is older than 30 days
                if os.path.isfile(os.path.join(new_dir, basename)): # file archived already?
                    # if it is, remove the file from the main folder of the site directory
                    logging.debug('{} is older than 30 days and already archived in {}. Deleting from source directory.'.format(basename, new_dir))
                    os.remove(filename)
                else:
                    # move it to the subdirectory, yyyy_mm
                    logging.info('{} is older than 30 days and is not archived. Moving to {}'.format(basename,new_dir))
                    shutil.move(filename, new_dir)

            else:  # file is newer than 30 days old.
                if os.path.isfile(os.path.join(new_dir, basename)): # is the file archived already?
                    # if it is, leave it in the main folder of the site directory
                    logging.debug('{} already archived.'.format(filename))
                    continue
                else:  # file is not archived already
                    # copy the file to the yyyy_mm directory, but leave a copy in the main site directory
                    logging.info('Archiving {} to {}'.format(filename, new_dir))
                    shutil.copy(filename, new_dir)

