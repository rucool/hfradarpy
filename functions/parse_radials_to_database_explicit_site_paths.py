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
            return 'File uploaded successfully'
        except:
            return 'File failed to upload'


if __name__ == '__main__':
    time_delta = 30  # days
    paths = ['/home/codaradm/data/radials/AMAG',
             '/home/codaradm/data/radials/AMHE',
             '/home/codaradm/data/radials/ASSA',
             '/home/codaradm/data/radials/ASVT',
             '/home/codaradm/data/radials/BELM',
             '/home/codaradm/data/radials/BEOP',
             '/home/codaradm/data/radials/BESE',
             '/home/codaradm/data/radials/BESP',
             '/home/codaradm/data/radials/BISL',
             '/home/codaradm/data/radials/BLCK',
             '/home/codaradm/data/radials/BRAD',
             '/home/codaradm/data/radials/BRBR',
             '/home/codaradm/data/radials/BRBY',
             '/home/codaradm/data/radials/BRIG',
             '/home/codaradm/data/radials/BRLO',
             '/home/codaradm/data/radials/BRMR',
             '/home/codaradm/data/radials/BRNT',
             '/home/codaradm/data/radials/BRRA',
             '/home/codaradm/data/radials/BRSE',
             '/home/codaradm/data/radials/BRSP',
             '/home/codaradm/data/radials/BRWI',
             '/home/codaradm/data/radials/BRZY',
             '/home/codaradm/data/radials/BSWP',
             '/home/codaradm/data/radials/CAPE',
             '/home/codaradm/data/radials/CBBT',
             '/home/codaradm/data/radials/CDDO',
             '/home/codaradm/data/radials/CEDR',
             '/home/codaradm/data/radials/CLUB',
             '/home/codaradm/data/radials/CMPT',
             '/home/codaradm/data/radials/CORE',
             '/home/codaradm/data/radials/CPHN',
             '/home/codaradm/data/radials/CStM',
             '/home/codaradm/data/radials/DUCK',
             '/home/codaradm/data/radials/ERRA',
             '/home/codaradm/data/radials/FARO',
             '/home/codaradm/data/radials/GCAP',
             '/home/codaradm/data/radials/GMNB',
             '/home/codaradm/data/radials/GRNI',
             '/home/codaradm/data/radials/HATY',
             '/home/codaradm/data/radials/HEAM',
             '/home/codaradm/data/radials/HEMP',
             '/home/codaradm/data/radials/HLPN',
             '/home/codaradm/data/radials/HOMR',
             '/home/codaradm/data/radials/HOOK',
             '/home/codaradm/data/radials/HOSR',
             '/home/codaradm/data/radials/JOUB',
             '/home/codaradm/data/radials/LISL',
             '/home/codaradm/data/radials/LOBR',
             '/home/codaradm/data/radials/LOHO',
             '/home/codaradm/data/radials/LOOK',
             '/home/codaradm/data/radials/LOVE',
             '/home/codaradm/data/radials/LPWR',
             '/home/codaradm/data/radials/MABO',
             '/home/codaradm/data/radials/METS',
             '/home/codaradm/data/radials/MISQ',
             '/home/codaradm/data/radials/MNTK',
             '/home/codaradm/data/radials/MRAM',
             '/home/codaradm/data/radials/MRCH',
             '/home/codaradm/data/radials/MRHE',
             '/home/codaradm/data/radials/MVBL',
             '/home/codaradm/data/radials/MVCO',
             '/home/codaradm/data/radials/MVNA',
             '/home/codaradm/data/radials/NANT',
             '/home/codaradm/data/radials/NAUS',
             '/home/codaradm/data/radials/OLDB',
             '/home/codaradm/data/radials/P125',
             '/home/codaradm/data/radials/P313',
             '/home/codaradm/data/radials/PALM',
             '/home/codaradm/data/radials/POOL',
             '/home/codaradm/data/radials/PORT',
             '/home/codaradm/data/radials/POSI',
             '/home/codaradm/data/radials/PYFA',
             '/home/codaradm/data/radials/PYFC',
             '/home/codaradm/data/radials/PYMA',
             '/home/codaradm/data/radials/RATH',
             '/home/codaradm/data/radials/RAWO',
             '/home/codaradm/data/radials/SEAB',
             '/home/codaradm/data/radials/SEOP',
             '/home/codaradm/data/radials/SESP',
             '/home/codaradm/data/radials/SET1',
             '/home/codaradm/data/radials/SILD',
             '/home/codaradm/data/radials/SLTR',
             '/home/codaradm/data/radials/SPAD',
             '/home/codaradm/data/radials/SPNT',
             '/home/codaradm/data/radials/SPOP',
             '/home/codaradm/data/radials/SPRK',
             '/home/codaradm/data/radials/SQUB',
             '/home/codaradm/data/radials/STLI',
             '/home/codaradm/data/radials/SUNS',
             '/home/codaradm/data/radials/TEST',
             '/home/codaradm/data/radials/VIEW',
             '/home/codaradm/data/radials/WAUW',
             '/home/codaradm/data/radials/WILD',
             '/home/codaradm/data/radials/WOOD']

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
        # max_workers = 16
        # with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        #     for recent, result in zip(recents, executor.map(parse_radial_file, recents)):
        #         logging.info('{} - {}'.format(recent, result))