#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR radial files utilizing the Radial class and upload to MySQL database.
"""
import codar_processing.database_common as db
import codar_processing.database_radials as dbr
import datetime as dt
import logging
import os
import sys
from configs.database_tables import RadialMetadata, RadialDiagnostics, HardwareDiagnostics
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
    r = Radial(radial_file)

    if not r.is_valid():
        return

    r.validate_header()  # Clean up header information for entry into mysql database
    r.header['filename'] = os.path.splitext(os.path.basename(radial_file))[0]
    r.header['fileModTime'] = dt.datetime.fromtimestamp(os.stat(fname).st_mtime)

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
        r.diags_radial = r.diags_radial.drop(['%%', 'TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'], axis=1)
        r.diags_radial['id_site'] = r.header['Site']
        r.diags_radial['id_radial'] = radial_id

        dbr.upload_diagnostics(session, RadialDiagnostics, r.diags_radial, r.header['Site'])
        logging.info('{} - Table `{}` - Diagnostic data uploaded '.format(r.header['filename'], 'hfrRadialDiagnostics'))
    except:
        pass

    try:
        # Upload hardware diagnostic data
        r.diags_hardware = r.diags_hardware.drop(['%%', 'TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'], axis=1)
        r.diags_hardware['id_site'] = r.header['Site']
        r.diags_hardware['id_radial'] = radial_id
        dbr.upload_diagnostics(session, HardwareDiagnostics, r.diags_hardware, r.header['Site'])
    except:
        pass


if __name__ == '__main__':
    radial_dir = '/home/codaradm/data/radials/*/'
    initial_loading = True
    time_delta = 30  # days

    from glob import glob
    paths = glob(radial_dir)
    now = dt.datetime.now()
    ago = now - dt.timedelta(days=time_delta)

    for site_path in paths:
        for fname in sorted(glob(os.path.join(site_path, '**', '*.ruv'), recursive=True)):
            st = os.stat(fname)
            mtime = dt.datetime.fromtimestamp(st.st_mtime)
            if mtime > ago:
                logging.debug('{} modified during the past {} days'.format(fname, mtime))
                logging.debug('{} - Checking if file is uploaded to MySQL database.'.format(fname))
                uploaded = db.check_file_upload(session, fname, RadialMetadata)
                if not uploaded:  # Check if the file has been uploaded already. If it hasn't, upload it completely.
                    logging.debug('{} - Loading'.format(fname))
                    try:
                        parse_radial_file(fname)
                    except:
                        session.rollback()
                        continue