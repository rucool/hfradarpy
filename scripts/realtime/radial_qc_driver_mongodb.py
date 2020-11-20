#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR radial files utilizing the Radial subclass and download qc values from MySQL database
"""
import concurrent.futures
import glob
import logging
import os
import sys
import time
import datetime as dt
from pymongo import MongoClient
from hfradar.methods.radials.qc_radial_file import main as qc_radial
from hfradar.src.common import list_to_dataframe

try:
    from hfradar.configs.configs import mongodb_configs
except ModuleNotFoundError:
    from hfradar.configs.configs_default import mongodb_configs

start_time = time.time()

# Set up the logger
logger = logging.getLogger(__name__)
log_level = 'ERROR'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


def qc_data(radial, qc_arguments):
    qc_arguments['radial_file'] = radial
    qc_radial(**qc_arguments)


if __name__ == '__main__':
    radial_dir = '/home/codaradm/data/radials/'
    save_dir = '/home/codaradm/data/radials_qc/'

    # List of sites to check. If empty, script finds all available sites on the fileserver and runs on all
    sites = ['AMAG', 'ASSA', 'BLCK', 'BRIG', 'CEDR', 'CORE', 'DUCK', 'HATY', 'HEMP', 'HOOK', 'LISL', 'LOVE', 'MRCH',
             'NANT', 'NAUS', 'WILD', 'MVCO']

    # If start_time, end_time are set to None, script will instead process radial files from the most recent radial
    # n number of days in the past defined in the variable days_to_check
    start_datetime = None  #pd.Timestamp(2020, 3, 11, 18, 0, 0)
    end_datetime = None  #pd.Timestamp(2020, 3, 11, 20, 0, 0)
    days_to_check = 30
    parallel = False
    month_subfolders = True

    if (start_datetime is None) or (end_datetime is None):
        end_datetime = dt.datetime.utcnow()
        start_datetime = end_datetime - dt.timedelta(days=days_to_check)

    # If the list of sites is empty, we will run qc on every site code
    if not sites:
        paths = glob.glob('{}/*/'.format(radial_dir))
        sites = [x.strip('/').split('/')[-1] for x in paths]

    client = MongoClient(mongodb_configs()['uri'])
    db = client.codar

    if month_subfolders:
        dir_str = '*/{}*.ruv'
    else:
        dir_str = '*.ruv'

    for site in sites:
        try:
            qc_arguments = dict(qc_values=list(db.sites.find({'code': site}))[0]['qc_settings'],
                                save_path=os.path.join(save_dir, site))
        except IndexError:
            print('{} does not exist in `site` collection'.format(site))
            pass

        site_dir = os.path.join(radial_dir, site)
        ptype = list(db.radial_status.find({'_id': site}))[0]['preferred_type']

        if ptype == 'N/A':
            files = glob.glob(os.path.join(site_dir, dir_str), recursive=True)
        else:
            files = glob.glob(os.path.join(site_dir, dir_str.format(ptype)), recursive=True)
        df = list_to_dataframe(files).sort_index()
        df = df[start_datetime: end_datetime]
        file_list = df.file.to_list()
        if parallel:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                zip(file_list, executor.map(qc_data, file_list))
        else:
            for f in file_list:
                qc_data(f, qc_arguments)

    elapsed_time = time.time() - start_time
    logging.info('Radial QC Complete. {} - seconds elapsed from start to finish.'.format(elapsed_time))
