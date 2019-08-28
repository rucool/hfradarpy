#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR radial files utilizing the Radial subclass and download qc values from MySQL database
"""
import codar_processing.src.database_common as db
import concurrent.futures
import datetime as dt
import glob
import logging
import os
import sys
import time
from codar_processing.configs.database_tables import Sites, QCValues
from codar_processing.methods.radials.qc_radial_file import main as qc_radial

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'ERROR'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


def qc_data(radial):
    qc_arguments['radial_file'] = radial
    st = os.stat(radial)
    mtime = dt.datetime.fromtimestamp(st.st_mtime)
    if mtime > ago:
        logging.info('{} modified during the past {} days: {}'.format(radial, days_to_check, mtime))
        qc_radial(**qc_arguments)


# List of sites to check. If left empty, we will find all available sites on the fileserver and run on everything
# sites = ['AMAG', 'BLCK', 'BRIG','HEMP', 'HOOK', 'LISL', 'LOVE', 'MRCH', 'MVCO', 'NANT', 'NAUS',' WILD', 'CEDR', 'CORE', 'DUCK', 'HATY', ]
sites = ['AMAG', 'CEDR', 'CORE', 'DUCK', 'HATY', 'HEMP', 'HOOK', 'LISL', 'LOVE', 'MRCH', 'MVCO', 'NANT', 'NAUS', 'WILD']
radial_dir = '/Volumes/home/codaradm/data_reprocessed/radials/'
save_dir = '/Users/mikesmith/Documents/radials_qc/'
# parallel = True
days_to_check = 5000

# Open up database connection. Database configuration is in ~/configs/configs.py
global session
session = db.db_session()

# If the list of sites is empty, we will run qc on every site code
if not sites:
    paths = glob.glob('{}/*/'.format(radial_dir))
    sites = [x.strip('/').split('/')[-1] for x in paths]

now = dt.datetime.now()
ago = now-dt.timedelta(days=days_to_check)

# Query the MySQL database for QC Values
results = session.query(QCValues, Sites).join(Sites).filter(Sites.site.in_(sites)).all()

# Create dictionary of qc values for easy
qc_values = {}
for _q,_s in results:
    qc_values[_s.site] = dict(
        qc_values=dict(
            qc_qartod_radial_count=dict(radial_min_count=_q.radial_min_count, radial_low_count=_q.radial_low_count),
            qc_qartod_maximum_velocity=dict(radial_max_speed=_q.radial_max_speed),
            qc_qartod_spatial_median=dict(radial_smed_range_cell_limit=_q.radial_smed_range_cell_limit,
                                          radial_smed_angular_limit=_q.radial_smed_angular_limit,
                                          radial_smed_current_difference=_q.radial_smed_current_difference)))

start_time = time.time()
for site in qc_values.keys():
    qc_arguments = qc_values[site]
    qc_arguments['save_path'] = os.path.join(save_dir, site)
    site_dir = os.path.join(radial_dir, site)
    files = sorted(glob.glob(os.path.join(site_dir, '*/*.ruv'), recursive=True))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        zip(files, executor.map(qc_data, files))

elapsed_time = time.time() - start_time
logging.info('Radial QC Complete. {} - seconds elapsed from start to finish.'.format(elapsed_time))