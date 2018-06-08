#!/usr/bin/env python
"""
@file clean_radial_dirs.py
@author Mike Smith
@email michaesm@marine.rutgers.edu
@brief Clean up radial directories
@purpose In a generic ~/data/radials/ directory, this script recursively searches for site subfolders and moves radial
files into the format ~/data/radials/site/yyyy-mm/. Script also leaves last 30 days in the main radial site folder
"""
import glob
import logging
import os
import shutil
import sys
from codar_processing.common import create_dir
from datetime import timedelta, datetime

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

source_dir = '/Volumes/home/codaradm/data/radials/'  # main directory that contains radial files
types = ('*.ruv', '*.euv')  # File types we want to move


# check if the date is greater than 30 days old.
def check_date(file):
    base_date = datetime.utcnow() - timedelta(days=30)
    yy = file.split("_")[2]
    yy = int(yy)
    mm = file.split("_")[3]
    mm = int(mm)
    dd = file.split("_")[4]
    dd = int(dd)
    HH = file.split("_")[5]
    HH = int(HH[:2])
    last_used = datetime(yy, mm, dd, HH, 0, 0, 0)
    if last_used < base_date:
        return True
    return False


# list the subdirectories from source_dir. Only add directories
site_dirs = [os.path.join(source_dir, o) for o in os.listdir(source_dir) if os.path.isdir(os.path.join(source_dir, o))]

# iterate through the site directories
for site in site_dirs:
    # create empty list for files_grabbed
    files_grabbed = []

    # grab a list of all files in each site directory that has all extension in the tuple, 'types', created above
    for files in types:
        files_grabbed.extend(glob.glob(os.path.join(site, files)))

    # if files_grabbed is not empty
    if files_grabbed:
        # iterate through each file in the list files_grabbed
        for filename in files_grabbed:
            # grab the filename from the path/filename combo, filena,e
            basename = os.path.basename(filename)
            base_name = basename.split('_')
            yy = base_name[2]
            mm = base_name[3]

            # create directory name in the format yyyy_mm
            new_dir = os.path.join(site, yy + '_' + mm)

            # try to create the directory. if it exists, pass
            create_dir(new_dir)

            # pass the filename to the check_date function to see if the file is older than 30 days
            found_file = check_date(basename)

            # if file is older than 30 days
            if found_file:
                # is the file archived already?
                if os.path.isfile(os.path.join(new_dir, basename)):
                    # if it is, remove the file from the main folder of the site directory
                    logging.debug('{} is older than 30 days and already archived in {}. Deleting from source directory.'.format(basename, new_dir))
                    os.remove(filename)

                # if the file isn't archived already,
                else:
                    # move it to the subdirectory, yyyy_mm
                    logging.info('{} is older than 30 days and is not archived. Moving to {}'.format(basename,new_dir))
                    shutil.move(filename, new_dir)

            # if the file is newer than 30 days old.
            else:
                # is the file archived already?
                if os.path.isfile(os.path.join(new_dir, basename)):
                    # if it is, leave it in the main folder of the site directory
                    logging.debug('{} already archived.'.format(filename))
                    continue
                # if the file isn't archived already
                else:
                    # copy the file to the yyyy_mm directory, but leave a copy in the main site directory
                    logging.info('Archiving {} to {}'.format(filename, new_dir))
                    shutil.copy(filename, new_dir)