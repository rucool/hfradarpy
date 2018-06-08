#!/usr/bin/env python
"""
@file update_frequencies.py
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Update CODAR site center frequencies in the hfrSites table of the RUCOOL MySQL database
"""
import codar_processing.database_common as dbc
import logging
import os.path
import sys
from codar_processing.radials import Radial
from decimal import *
from glob import glob

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

# Initialize sqlalchemy session with codar MySQL database
session = dbc.db_session()


def update_frequency(session, site, frequency):
    """
    Check frequency of a site. It it changes, update it.
    :param session: sqlalchemy session binding
    :param site: four digit site code
    :return:
    """
    result = dbc.site_check(session, site)
    site_freq_db = result.transmitCenterFrequency
    center_freq = Decimal(frequency)

    if not site_freq_db == center_freq:
        logging.info('{}: Updating frequency from {} MHz to {} MHz'.format(site, site_freq_db,  frequency))
        result.transmitCenterFrequency = center_freq
        session.commit()
    else:
        logging.info('{}: No frequency change required'.format(site))
    return


sourceDir = '/home/codaradm/data/radials/*/'
site_paths = glob(sourceDir)

for site_dir in site_paths:
    file_list = [os.path.join(site_dir, f) for f in os.listdir(site_dir) if f.startswith('RDL')]
    if file_list:
        latest_file = max(file_list, key=os.path.getctime)

        r = Radial(latest_file)
        update_frequency(session, r.header['Site'].strip(' ""'), r.header['TransmitCenterFreqMHz'])