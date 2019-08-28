#!/usr/bin/env python
import datetime as dt
import logging
import os
import sys

import codar_processing.src.database_common as db
import codar_processing.src.database_waves as dbw
from codar_processing.src.waves import Waves
from codar_processing.configs import database_tables

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

# Initialize sqlalchemy session with codar MySQL database
session = db.db_session()

# Load some relational ids upon initial execution of script. Saves time because less database calls.
radar_sites = db.get_sites(session)
pattern_types = db.get_pattern_types(session)


def parse_wave_file(wave_file, upload_file_id=None):
    """

    :param wave_file: wave file name and location
    :param upload_file_id: optional file id upload.
    :return:
    """
    mode = 1
    w = Waves(wave_file, multi_dimensional=False)

    if w.data['DIST'].isnull().all():
        mode = 2
        w.data = w.data.drop(['DIST'], axis=1)

    w.data = w.data.dropna()  # Drop columns containing 999 or 1080 values which are replaced by NaN in the Wave class.
    w.flag_wave_heights(wave_min=0.2, wave_max=5)  # Flag wave heights between these two values
    w.clean_wave_header()  # Clean up header information for entry into mysql database
    w.metadata['filename'] = os.path.splitext(os.path.basename(wave_file))[0]
    w.metadata['TableWaveMode'] = mode

    # Fill certain table columns with relational ids
    # Check to see if the site has been uploaded to the HfrSites table of the MySQL database
    try:
        site_info = radar_sites[radar_sites.site == w.metadata['Site']]
        site_id = int(site_info.id.iloc[0])
    except IndexError:
        logging.info('{} not found. Uploading site to hfrSites table'.format(w.metadata['Site']))
        site_info = db.update_site_table(session, w.metadata['Site'], w.metadata['TransmitCenterFreqMHz'], w.metadata['Origin'])
        site_id = int(site_info)

    w.metadata['Site'] = site_id

    # Upload file data to database
    if not upload_file_id:
        # Upload file information to database
        ref = database_tables.WaveFileHeader(**w.metadata)
        session.add(ref)
        session.commit()
        session.flush()
        ref_id = ref.id
        dbw.data_route(session, w.data, w.metadata['Site'], ref_id, fname, initial_upload=True)
    else:
        dbw.data_route(session, w.data, w.metadata['Site'], uploaded.id, fname, initial_upload=False)

    logging.info('{} - Upload complete.'.format(fname))


if __name__ == '__main__':
    from glob import glob

    wave_dir = '/home/codaradm/data/waves/'
    time_delta = 30
    site_codes = []

    ago = dt.datetime.now() - dt.timedelta(days=time_delta)

    if site_codes:
        paths = [os.path.join(wave_dir, x) for x in site_codes]
    else:
        paths = [os.path.join(wave_dir, o) for o in os.listdir(wave_dir) if os.path.isdir(os.path.join(wave_dir, o))]

    for site_path in paths:
        logging.info('Checking {} wave archive for files modified in the past {} days'.format(site_path.split(wave_dir)[1], time_delta))
        file_list = sorted(glob(os.path.join(site_path, '*.wls')), key=os.path.getmtime)

        for fname in file_list:
            st = os.stat(fname)
            mtime = dt.datetime.fromtimestamp(st.st_mtime)
            if mtime > ago:
                logging.debug('{} - Checking if file is uploaded to MySQL database.'.format(fname))
                uploaded = db.check_file_upload(session, fname, database_tables.WaveFileHeader)
                if not uploaded:  # Check if the file has been uploaded already. If it hasn't, upload it completely.
                    logging.debug('{} - Loading'.format(fname))
                    parse_wave_file(fname)
                else:  # If it has, check to see if the file has been updated since the last time it was uploaded.
                    logging.debug('{} - Already uploaded.'.format(fname))
                    one_month_ago = dt.datetime.now() - dt.timedelta(days=time_delta)
                    mod_time = dt.datetime.utcfromtimestamp(os.path.getmtime(fname))

                    if mod_time > one_month_ago:
                        logging.info('{} - Modified within past month ({}). Updating database'.format(fname, mod_time))
                        parse_wave_file(fname, upload_file_id=uploaded.id)
