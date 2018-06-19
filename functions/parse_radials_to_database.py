#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR radial files utilizing the Radial class and upload to MySQL database.
"""
import codar_processing.database_common as db
import codar_processing.database_radials as dbr
import concurrent.futures
import datetime as dt
from glob import glob
import logging
import os
import pandas as pd
import sys
from configs.database_tables import RadialMetadata, RadialDiagnostics, HardwareDiagnostics
from codar_processing.common import timestamp_from_lluv_filename
from codar_processing.radials import Radial

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

# Initialize sqlalchemy session with codar MySQL database
session = db.db_session()


def parse_radial_file(radial_file):
    """
    Parse CODAR radial files utilizing the Radial class and upload to MySQL database.
    :param radial_file: Path to CODAR Radial File
    """
    basename = os.path.basename(radial_file).split('.')[0]
    logging.debug('{} - Checking if file is uploaded to MySQL database.'.format(basename))
    uploaded = db.check_file_upload(session, basename, RadialMetadata)
    if not uploaded:  # Check if the file has been uploaded already. If it hasn't, upload it completely.
        logging.debug('{} - Loading'.format(radial_file))
        try:
            r = Radial(radial_file)

            if not r.is_valid():
                return

            r.validate_header()  # Clean up header information for entry into mysql database
            r.header['filename'] = os.path.splitext(os.path.basename(radial_file))[0]
            r.header['fileModTime'] = dt.datetime.fromtimestamp(os.stat(radial_file).st_mtime)

            # Fill certain table columns with relational ids
            # Check to see if the site has been uploaded to the HfrSites table of the MySQL database and get site_id info
            site_info = db.site_check(session, r.header['Site'], r.header['TransmitCenterFreqMHz'], r.header['Origin'])
            r.header['Site'] = site_info.id
            r.header['PatternType'] = dbr.get_pattern_type_id(session, r.header['PatternType'])

            # Add extra information to header
            r.header['TableType'] = r.tables['1']['TableType']
            r.header['TableColumns'] = r.tables['1']['TableColumns']
            r.header['TableColumnTypes'] = r.tables['1']['TableColumnTypes']
            r.header['TableRows'] = r.tables['1']['TableRows']

            # Upload radial header information and update latest radials table
            radial_id = dbr.upload_radial_header(session, r.header)
            dbr.update_latest_radials(session, r.header['filename'], r.header['TimeStamp'], r.header['Site'], radial_id)

            try:
                # Upload radial diagnostic data
                r.diags_radial = r.diags_radial.drop(['%%', 'TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'],
                                                     axis=1)
                r.diags_radial['id_site'] = r.header['Site']
                r.diags_radial['id_radial'] = radial_id

                dbr.upload_diagnostics(session, RadialDiagnostics, r.diags_radial, r.header['Site'])
                logging.info(
                    '{} - Table `{}` - Diagnostic data uploaded '.format(r.header['filename'], 'hfrRadialDiagnostics'))
            except:
                pass

            try:
                # Upload hardware diagnostic data
                r.diags_hardware = r.diags_hardware.drop(['%%', 'TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'],
                                                         axis=1)
                r.diags_hardware['id_site'] = r.header['Site']
                r.diags_hardware['id_radial'] = radial_id
                dbr.upload_diagnostics(session, HardwareDiagnostics, r.diags_hardware, r.header['Site'])
            except:
                pass
            return logging.info('{} - Uploaded successfully'.format(r.header['filename']))
        except:
            return logging.error('{} - File failed to upload'.format(basename))


if __name__ == '__main__':
    radial_dir = '/home/codaradm/data/radials/'
    initial_loading = True
    time_delta = 535  # days

    paths = [os.path.join(radial_dir, o) for o in os.listdir(radial_dir) if os.path.isdir(os.path.join(radial_dir, o))]
    now = dt.datetime.now()
    ago = now - dt.timedelta(days=time_delta)

    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        for site_path in paths:
            file_list = glob(os.path.join(site_path, '**', '*.ruv'), recursive=True)
            # Convert to dataframe for faster index searching
            df = pd.DataFrame(file_list, columns=['path'])
            df['timestamp'] = df['path'].apply(timestamp_from_lluv_filename)
            tdf = df[df['timestamp'] > ago]

            # Convert back to list and sort by filename
            recents = sorted(tdf['path'].tolist())

            res = executor.map(parse_radial_file, recents)