#!/usr/bin/env python
import datetime as dt
import codar_processing.database_common as db
import codar_processing.database_waves as dbw
import glob
import logging
import os
import sys
from configs import database_tables
from codar_processing.waves import Waves

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

wave_dir='/home/codaradm/data/waves/'
sites = ['BRAD', 'BRMR', 'BRNT', 'RATH', 'SEAB', 'SPRK', 'WOOD', 'BISL', 'CMPT', 'CAPE', 'CPHN', 'GCAP', 'HLPN', 'MISQ',
         'MNTK', 'OLDB', 'PALM', 'PORT', 'SILD', 'STLI', 'SUNS', 'VIEW']

global session
session = db.db_session()


def parse_wave_file(wave_file, upload_file_id=None):
    """

    :param wave_file: wave file name and location
    :param upload_file_id: optional file id upload.
    :return:
    """
    mode = 1
    w = Waves(wave_file)

    if w.data['DIST'].isnull().all():
        mode = 2
        w.data = w.data.drop(['DIST'], axis=1)

    w.remove_bad_data()  # Drop columns containing 999 or 1080 values which are replaced by NaN in the Wave class.
    w.flag_wave_heights(wave_min=0.2, wave_max=5)  # Flag wave heights between these two values
    w.clean_wave_header()  # Clean up header information for entry into mysql database
    w.header['filename'] = os.path.splitext(os.path.basename(wave_file))[0]
    w.header['TableWaveMode'] = mode

    # Check to see if the site has been uploaded to the HfrSites table of the MySQL database and get site_id info
    site_info = db.site_check(session, w.header['Site'], w.header['TransmitCenterFreqMHz'], w.header['Origin'])
    w.header['Site'] = site_info.id

    # Upload file data to database
    if not upload_file_id:
        # Upload file information to database
        ref = database_tables.WaveFile(**w.header)
        session.add(ref)
        session.commit()
        session.flush()
        ref_id = ref.id
        dbw.data_route(session, w.data, w.header['Site'], ref_id, fname, initial_upload=True)
    else:
        dbw.data_route(session, w.data, w.header['Site'], uploaded.id, fname, initial_upload=False)

    logging.info('{} - Upload complete.'.format(fname))


for site in sites:
    site_dir = os.path.join(wave_dir, site)
    for fname in sorted(glob.glob(os.path.join(site_dir, '*.wls')), key=os.path.getmtime):
        logging.debug('{} - Checking if file is uploaded to MySQL database.'.format(fname))
        uploaded = db.check_file_upload(session, fname, database_tables.WaveFile)
        if not uploaded:  # Check if the file has been uploaded already. If it hasn't, upload it completely.
            logging.debug('{} - Loading'.format(fname))
            parse_wave_file(fname)
        else:  # If it has, check to see if the file has been updated since the last time it was uploaded.
            logging.debug('{} - Already uploaded.'.format(fname))
            one_month_ago = dt.datetime.now() - dt.timedelta(days=30)
            mod_time = dt.datetime.utcfromtimestamp(os.path.getmtime(fname))

            if mod_time > one_month_ago:
                logging.info('{} - Modified within past month ({}). Updating database'.format(fname, mod_time))
                parse_wave_file(fname, upload_file_id=uploaded.id)