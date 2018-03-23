#!/usr/bin/env python
import datetime as dt
import functions.database_common as db
import functions.database_waves as dbw
import functions.common as cf
import functions.waves as waves
import glob
import logging
import os
import pandas as pd
import sys
from configs import database_tables

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


def parse_wave_file(fname):
    """
    Function to parse the wave file
    :param fname: directory/filename of the wave file
    :return: Nothing
    """
    mode = 1
    with open(fname, 'r') as wave_file:
        header_data, distance_data = cf.parse_header(wave_file)  # Parse out header metadata above the data table
        header_data = waves.clean_wave_header(header_data)  # Clean the header data for input into MySQL database.
        header_data['filename'] = os.path.splitext(os.path.basename(fname))[0]
        header = header_data['TableColumnTypes'].split()

        try:
            data = pd.read_csv(fname, comment='%', sep=' ', header=None, names=header, skipinitialspace=True)

            if data['DIST'].isnull().all():
                mode = 2
                data = data.drop(['DIST'], axis=1)
        except BaseException:
            logger.warning('{} - Could not load. Skipping.'.format(fname))

        data['datetime'] = data[['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)

        header_data['TableWaveMode'] = mode
        site_info = db.site_check(session, header_data['Site'], header_data['TransmitCenterFreqMHz'], header_data['Origin'])
        header_data['Site'] = site_info.id
        ref = database_tables.WaveFile(**header_data)
        session.add(ref)
        session.commit()
        session.flush()
        ref_id = ref.id

        dbw.data_route(session, data, site_info, ref_id, fname, initial_upload=True)
        logging.info('{} - Upload complete.'.format(fname))


def main(wave_dir='/home/codaradm/data/waves/', sites=None, database=True):
    if database:
        global session
        session = db.db_session()

    if not sites:
        # Use the following sites as the defaults
        sites = ['BRAD', 'BRMR', 'BRNT', 'RATH', 'SEAB', 'SPRK', 'WOOD', 'BISL', 'CMPT', 'CAPE', 'CPHN', 'GCAP', 'HLPN',
                 'MISQ', 'MNTK', 'OLDB', 'PALM', 'PORT', 'SILD', 'STLI', 'SUNS', 'VIEW']

    for site in sites:
        site_dir = os.path.join(wave_dir, site)
        for fname in sorted(glob.glob(os.path.join(site_dir, '*.wls')), key=os.path.getmtime):
            logging.debug('{} - Checking if file is uploaded to MySQL database.'.format(fname))
            uploaded = db.check_file_upload(session, fname, database_tables.WaveFile)
            if not uploaded:
                logging.debug('{} - Loading'.format(fname))
                parse_wave_file(fname)
            else:
                logging.debug('{} - Already uploaded.'.format(fname))
                one_month_ago = dt.datetime.now() - dt.timedelta(days=30)
                mod_time = dt.datetime.utcfromtimestamp(os.path.getmtime(fname))
                site_info = db.site_check(session, os.path.splitext(os.path.basename(fname))[0].split('_')[1])

                if mod_time > one_month_ago:
                    logging.info('{} - Modified within past month ({}). Updating database'.format(fname, mod_time))

                    with open(fname) as wave_file:
                        wave_data = wave_file.readlines()
                        for line in wave_data:
                            if 'TableColumnTypes' in line:
                                split_line = cf.parse_header_line(line)
                                header = split_line[1].split()
                                break
                    try:
                        data = pd.read_csv(fname, comment='%', sep=' ', header=None, names=header, skipinitialspace=True)

                        if data['DIST'].isnull().all():
                            data = data.drop(['DIST'], axis=1)

                        data['datetime'] = data[['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)

                        dbw.data_route(session, data, site_info, uploaded.id, fname, initial_upload=False)
                    except BaseException:
                        logger.error('{} - Could not load. Skipping.'.format(fname))


if __name__ == '__main__':
    wave_data_directory = cf.path_within_module('data/waves')
    try:
        exit(main(wave_data_directory))
    except Exception:
        logger.exception('Exception in main(): ')
        exit(1)