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
import pandas as pd
from pymongo import MongoClient
from codar_processing.methods.radials.qc_radial_file import main as qc_radial
from codar_processing.src.common import list_to_dataframe

start_time = time.time()

# Set up the logger
logger = logging.getLogger(__name__)
log_level = 'ERROR'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


def qc_data(radial):
    qc_arguments['radial_file'] = radial
    qc_radial(**qc_arguments)


if __name__ == '__main__':
    radial_dir = '/home/codaradm/data/radials/'
    save_dir = '/home/codaradm/data/radials_qc/'

    # List of sites to check. If empty, script finds all available sites on the fileserver and runs on all
    sites = ['AMAG', 'ASSA', 'BLCK', 'BRIG', 'CEDR', 'CORE', 'DUCK', 'HATY', 'HEMP', 'HOOK', 'LISL', 'LOVE', 'MRCH',
             'NANT', 'NAUS', 'WILD']

    # If start_time, end_time are set to None, script will instead process radial files from the most recent radial
    # n number of days in the past defined in the variable days_to_check
    start_datetime = None
    end_datetime = None
    days_to_check = 7
    parallel = True

    if (start_datetime is None) or (end_datetime is None):
        end_datetime = pd.Timestamp.utcnow()
        start_datetime = end_datetime - pd.Timedelta(days=days_to_check)

    # If the list of sites is empty, we will run qc on every site code
    if not sites:
        paths = glob.glob('{}/*/'.format(radial_dir))
        sites = [x.strip('/').split('/')[-1] for x in paths]

    client = MongoClient('mongodb://michaesm:sslctob22vDH@mongodb.marine.rutgers.edu:27017')
    db = client.codar

    for site in sites:
        try:
            qc_arguments = dict(qc_values=list(db.sites.find({'code': site}))[0]['qc_settings'],
                                save_path=os.path.join(radial_dir, site))
        except IndexError:
            print('{} does not exist in `site` collection'.format(site))
            pass

        site_dir = os.path.join(radial_dir, site)
        files = glob.glob(os.path.join(site_dir, '*/*.ruv'), recursive=True)
        df = list_to_dataframe(files).sort_index()
        df = df[start_datetime: end_datetime]
        if parallel:
            file_list = df.file.to_list()
            with concurrent.futures.ProcessPoolExecutor() as executor:
                zip(file_list, executor.map(qc_data, file_list))
        else:
            for row in df.itertuples():
                qc_data(row.file, qc_arguments)

    elapsed_time = time.time() - start_time
    logging.info('Radial QC Complete. {} - seconds elapsed from start to finish.'.format(elapsed_time))