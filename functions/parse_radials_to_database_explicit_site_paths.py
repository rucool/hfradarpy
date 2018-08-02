#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR radial files utilizing the Radial class and upload to MySQL database.
"""
import codar_processing.database_common as db
import codar_processing.database_radials as dbr
# import concurrent.futures
import datetime as dt
import logging
import os
import pandas as pd
import sys
from configs.database_tables import RadialMetadata, RadialDiagnostics, HardwareDiagnostics
from codar_processing.common import timestamp_from_lluv_filename
from codar_processing.radials import Radial
from glob import glob

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

            r.validate_header()  # Clean up header information for entry into mysql database
            r.header['filename'] = os.path.splitext(os.path.basename(radial_file))[0]
            r.header['fileModTime'] = dt.datetime.fromtimestamp(os.stat(radial_file).st_mtime)

            # Fill certain table columns with relational ids
            # Check to see if the site has been uploaded to the HfrSites table of the MySQL database
            try:
                site_info = sites[sites.site == r.header['Site']]
                site_id = int(site_info.id.iloc[0])
            except IndexError:
                logging.info('{} not found. Uploading site to hfrSites table'.format(r.header['Site']))
                site_info = db.update_site_table(session, r.header['Site'], r.header['TransmitCenterFreqMHz'], r.header['Origin'])
                site_id = int(site_info)

            r.header['Site'] = site_id

            try:
                patt_type = pattern_types[pattern_types.type == r.header['PatternType']]
                pattern_id = int(patt_type.id.iloc[0])
            except IndexError:
                logging.info('{} not found. Pattern type invalid'.format(r.header['PatternType']))
                return

            r.header['PatternType'] = pattern_id

            # Add extra information to header
            r.header['TableType'] = r.tables['1']['TableType']
            r.header['TableColumns'] = r.tables['1']['TableColumns']
            r.header['TableColumnTypes'] = r.tables['1']['TableColumnTypes']
            r.header['TableRows'] = r.tables['1']['TableRows']

            # Upload radial header information and update latest radials table
            radial_id = dbr.upload_radial_header(session, r.header)

            try:
                # Upload radial diagnostic data
                r.diags_radial = r.diags_radial.drop(['%%', 'TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'],
                                                     axis=1)
                r.diags_radial['id_site'] = r.header['Site']
                r.diags_radial['id_radial'] = radial_id
                dbr.upload_diagnostics(session, RadialDiagnostics, r.diags_radial, r.header['Site'])
                logging.debug(
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
                logging.debug('{} - Table `{}` - Diagnostic data uploaded '.format(r.header['filename'],
                                                                                   'hfrHardwareDiagnostics'))
            except:
                pass
            logging.info('{} - File uploaded successfully'.format(radial_file))
        except:
            logging.error('{} - File failed to upload'.format(radial_file))


if __name__ == '__main__':
    time_delta = 190  # days
    paths = ['/Volumes/home/codaradm/data/radials/AMAG',
             '/Volumes/home/codaradm/data/radials/AMHE',
             '/Volumes/home/codaradm/data/radials/ASSA',
             '/Volumes/home/codaradm/data/radials/ASVT',
             '/Volumes/home/codaradm/data/radials/BELM',
             '/Volumes/home/codaradm/data/radials/BEOP',
             '/Volumes/home/codaradm/data/radials/BESE',
             '/Volumes/home/codaradm/data/radials/BESP',
             '/Volumes/home/codaradm/data/radials/BISL',
             '/Volumes/home/codaradm/data/radials/BLCK',
             '/Volumes/home/codaradm/data/radials/BRAD',
             '/Volumes/home/codaradm/data/radials/BRBR',
             '/Volumes/home/codaradm/data/radials/BRBY',
             '/Volumes/home/codaradm/data/radials/BRIG',
             '/Volumes/home/codaradm/data/radials/BRLO',
             '/Volumes/home/codaradm/data/radials/BRMR',
             '/Volumes/home/codaradm/data/radials/BRNT',
             '/Volumes/home/codaradm/data/radials/BRRA',
             '/Volumes/home/codaradm/data/radials/BRSE',
             '/Volumes/home/codaradm/data/radials/BRSP',
             '/Volumes/home/codaradm/data/radials/BRWI',
             '/Volumes/home/codaradm/data/radials/BRZY',
             '/Volumes/home/codaradm/data/radials/BSWP',
             '/Volumes/home/codaradm/data/radials/CAPE',
             '/Volumes/home/codaradm/data/radials/CBBT',
             '/Volumes/home/codaradm/data/radials/CDDO',
             '/Volumes/home/codaradm/data/radials/CEDR',
             '/Volumes/home/codaradm/data/radials/CLUB',
             '/Volumes/home/codaradm/data/radials/CMPT',
             '/Volumes/home/codaradm/data/radials/CORE',
             '/Volumes/home/codaradm/data/radials/CPHN',
             '/Volumes/home/codaradm/data/radials/CStM',
             '/Volumes/home/codaradm/data/radials/DUCK',
             '/Volumes/home/codaradm/data/radials/ERRA',
             '/Volumes/home/codaradm/data/radials/FARO',
             '/Volumes/home/codaradm/data/radials/GCAP',
             '/Volumes/home/codaradm/data/radials/GMNB',
             '/Volumes/home/codaradm/data/radials/GRNI',
             '/Volumes/home/codaradm/data/radials/HATY',
             '/Volumes/home/codaradm/data/radials/HEAM',
             '/Volumes/home/codaradm/data/radials/HEMP',
             '/Volumes/home/codaradm/data/radials/HLPN',
             '/Volumes/home/codaradm/data/radials/HOMR',
             '/Volumes/home/codaradm/data/radials/HOOK',
             '/Volumes/home/codaradm/data/radials/HOSR',
             '/Volumes/home/codaradm/data/radials/JOUB',
             '/Volumes/home/codaradm/data/radials/LISL',
             '/Volumes/home/codaradm/data/radials/LOBR',
             '/Volumes/home/codaradm/data/radials/LOHO',
             '/Volumes/home/codaradm/data/radials/LOOK',
             '/Volumes/home/codaradm/data/radials/LOVE',
             '/Volumes/home/codaradm/data/radials/LPWR',
             '/Volumes/home/codaradm/data/radials/MABO',
             '/Volumes/home/codaradm/data/radials/METS',
             '/Volumes/home/codaradm/data/radials/MISQ',
             '/Volumes/home/codaradm/data/radials/MNTK',
             '/Volumes/home/codaradm/data/radials/MRAM',
             '/Volumes/home/codaradm/data/radials/MRCH',
             '/Volumes/home/codaradm/data/radials/MRHE',
             '/Volumes/home/codaradm/data/radials/MVBL',
             '/Volumes/home/codaradm/data/radials/MVCO',
             '/Volumes/home/codaradm/data/radials/MVNA',
             '/Volumes/home/codaradm/data/radials/NANT',
             '/Volumes/home/codaradm/data/radials/NAUS',
             '/Volumes/home/codaradm/data/radials/OLDB',
             '/Volumes/home/codaradm/data/radials/P125',
             '/Volumes/home/codaradm/data/radials/P313',
             '/Volumes/home/codaradm/data/radials/PALM',
             '/Volumes/home/codaradm/data/radials/POOL',
             '/Volumes/home/codaradm/data/radials/PORT',
             '/Volumes/home/codaradm/data/radials/POSI',
             '/Volumes/home/codaradm/data/radials/PYFA',
             '/Volumes/home/codaradm/data/radials/PYFC',
             '/Volumes/home/codaradm/data/radials/PYMA',
             '/Volumes/home/codaradm/data/radials/RATH',
             '/Volumes/home/codaradm/data/radials/RAWO',
             '/Volumes/home/codaradm/data/radials/SEAB',
             '/Volumes/home/codaradm/data/radials/SEOP',
             '/Volumes/home/codaradm/data/radials/SESP',
             '/Volumes/home/codaradm/data/radials/SET1',
             '/Volumes/home/codaradm/data/radials/SILD',
             '/Volumes/home/codaradm/data/radials/SLTR',
             '/Volumes/home/codaradm/data/radials/SPAD',
             '/Volumes/home/codaradm/data/radials/SPNT',
             '/Volumes/home/codaradm/data/radials/SPOP',
             '/Volumes/home/codaradm/data/radials/SPRK',
             '/Volumes/home/codaradm/data/radials/SQUB',
             '/Volumes/home/codaradm/data/radials/STLI',
             '/Volumes/home/codaradm/data/radials/SUNS',
             '/Volumes/home/codaradm/data/radials/TEST',
             '/Volumes/home/codaradm/data/radials/VIEW',
             '/Volumes/home/codaradm/data/radials/WAUW',
             '/Volumes/home/codaradm/data/radials/WILD',
             '/Volumes/home/codaradm/data/radials/WOOD']

    now = dt.datetime.now()
    ago = now - dt.timedelta(days=time_delta)

    for site_path in paths:
        file_list = glob(os.path.join(site_path, '**', '*.ruv'), recursive=True)
        # Convert to dataframe for faster index searching
        df = pd.DataFrame(file_list, columns=['path'])
        df['timestamp'] = df['path'].apply(timestamp_from_lluv_filename)
        tdf = df[df['timestamp'] > ago]

        # Convert back to list and sort by filename
        recents = sorted(tdf['path'].tolist())

        for f in recents:
            parse_radial_file(f)

        # TODO
        # Add multiprocessing functionality that actually works.....
        # max_workers = 8
        # with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        #     for recent, result in zip(recents, executor.map(parse_radial_file, recents)):
        #         logging.info('success')