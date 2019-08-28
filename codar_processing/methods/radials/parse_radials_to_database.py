#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR radial files utilizing the Radial class and upload to MySQL database.
"""
import datetime as dt
import logging
import os
import sys

import pandas as pd

import codar_processing.src.database_common as db
import codar_processing.src.database_radials as dbr
from codar_processing.src.common import timestamp_from_lluv_filename
from codar_processing.src.radials import Radial
from codar_processing.configs.database_tables import RadialMetadata, RadialDiagnostics, HardwareDiagnostics

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

# Initialize sqlalchemy session with codar MySQL database
session = db.db_session()

# Load some relational ids upon initial execution of script. Saves time because less database calls.
sites = db.get_sites(session)
pattern_types = db.get_pattern_types(session)


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

            r.clean_header()  # Clean up header information for entry into mysql database
            r.metadata['filename'] = os.path.splitext(os.path.basename(radial_file))[0]
            r.metadata['fileModTime'] = dt.datetime.fromtimestamp(os.stat(radial_file).st_mtime)

            # Fill certain table columns with relational ids
            # Check to see if the site has been uploaded to the HfrSites table of the MySQL database
            try:
                site_info = sites[sites.site == r.metadata['Site']]
                site_id = int(site_info.id.iloc[0])
            except IndexError:
                logging.info('{} not found. Uploading site to hfrSites table'.format(r.metadata['Site']))
                site_info = db.update_site_table(session, r.metadata['Site'], r.metadata['TransmitCenterFreqMHz'], r.metadata['Origin'])
                site_id = int(site_info)

            r.metadata['Site'] = site_id

            try:
                patt_type = pattern_types[pattern_types.type == r.metadata['PatternType']]
                pattern_id = int(patt_type.id.iloc[0])
            except IndexError:
                logging.error('{} not found. Pattern type invalid'.format(r.metadata['PatternType']))
                return

            r.metadata['PatternType'] = pattern_id

            # Add extra information to header
            r.metadata['TableType'] = r._tables['1']['TableType']
            r.metadata['TableColumns'] = r._tables['1']['TableColumns']
            r.metadata['TableColumnTypes'] = r._tables['1']['TableColumnTypes']
            r.metadata['TableRows'] = r._tables['1']['TableRows']

            # Upload radial header information and update latest radials table
            r.metadata = dbr.upload_radial_header(session, r.metadata)
            dbr.update_latest_radials(session, r.metadata)

            try:
                # Upload radial diagnostic data
                r.diagnostics_radial = r.diagnostics_radial.drop(['TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'], axis=1)
                r.diagnostics_radial['id_site'] = r.metadata['Site']
                r.diagnostics_radial['id_radial'] = r.metadata['radial_id']
                dbr.upload_diagnostics(session, RadialDiagnostics, r.diagnostics_radial, r.metadata['Site'])
                logging.debug('{} - Table `{}` - Diagnostic data uploaded '.format(r.metadata['filename'], 'hfrRadialDiagnostics'))
            except:
                pass

            try:
                # Upload hardware diagnostic data
                r.diagnostics_hardware = r.diagnostics_hardware.drop(['TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'], axis=1)
                r.diagnostics_hardware['id_site'] = r.metadata['Site']
                r.diagnostics_hardware['id_radial'] = r.metadata['radial_id']
                dbr.upload_diagnostics(session, HardwareDiagnostics, r.diagnostics_hardware, r.metadata['Site'])
                logging.debug('{} - Table `{}` - Diagnostic data uploaded '.format(r.metadata['filename'], 'hfrHardwareDiagnostics'))
            except:
                pass
            logging.info('{} - File uploaded successfully'.format(radial_file))
        except:
            logging.error('{} - File failed to upload'.format(radial_file))


if __name__ == '__main__':
    from glob import glob

    radial_dir = '/home/codaradm/data/radials/'
    time_delta = 1  # days. Can be in fractions of a day.
    site_codes = []  # Leave empty if you want the script to determine radial folders in root directory
    recursive = False  # If radials are in a subdirectory of the main radial folder, change to True

    if site_codes:
        # sites must be explicitly defined above
        paths = [os.path.join(radial_dir, x) for x in site_codes]
    else:
        # Do every site folder in radial_dir
        paths = [os.path.join(radial_dir, o) for o in os.listdir(radial_dir) if os.path.isdir(os.path.join(radial_dir, o))]

    ago = dt.datetime.now() - dt.timedelta(days=time_delta)

    for site_path in paths:
        if recursive:
            file_list = sorted(glob(os.path.join(site_path, '**', '*.ruv')))
        else:
            file_list = sorted(glob(os.path.join(site_path, '*.ruv')))

        # Convert to dataframe for faster index searching
        df = pd.DataFrame(file_list, columns=['path'])
        df['timestamp'] = df['path'].apply(timestamp_from_lluv_filename)
        tdf = df[df['timestamp'] > ago]

        # Convert back to list and sort by filename
        recents = tdf['path'].tolist()

        for recent in recents:
            parse_radial_file(recent)
