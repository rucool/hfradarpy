import logging
import os
import sys
from codar_processing.src.radials import Radial
from pymongo import MongoClient, InsertOne, IndexModel, ASCENDING, DESCENDING
from pymongo.errors import BulkWriteError
from codar_processing.src.common import divide_chunks


def frequency_check(freq):
    if 0 < freq <= 10:
        type = 'long_range'
    elif 10 < freq <= 20:
        type = 'mid_range'
    elif 20 < freq <= 50:
        type = 'standard_range'
    else:
        logger.warning('Center Frequency not known. Filling with unknown type value')
        type = 'unknown_range'
    return type


# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


# Define indexes for collections. This optimizes searching through the
index1 = IndexModel([('filename', ASCENDING)], unique=True)  # filename name needs to be unique
index2 = IndexModel([('Site', ASCENDING),
                     ('PatternType', ASCENDING),
                     ('TimeStamp', DESCENDING),
                     ])


def main(file_list):
    client = MongoClient()
    db = client.codar
    db.radials.create_indexes([index1, index2])

    bulk_info = []
    for radial in file_list:  # TODO Add multiprocessing here.
        # print(radial)
        r = Radial(radial)
        # print(r)
        if r.is_valid():
            r.metadata['Site'] = r.metadata['Site']
            try:
                r.metadata['PatternType'] = r.metadata['PatternType'].lower()
            except KeyError:
                pass
            # print(r.file_name)
            r.clean_header(split_origin=True)
            r.metadata['filename'] = r.file_name

            # assign a system type so we can sort on this
            r.metadata['SystemType'] = frequency_check(r.metadata['TransmitCenterFreqMHz'])
            r.metadata['RadialSolutions'] = r.data.__len__()

            # Try statements in case a radial file doesn't contain a diagnostic table.
            try:
                r.metadata['diagnostics_hardware'] = r.diagnostics_hardware.to_dict(orient='r')
            except AttributeError as ae:
                logging.error(ae)
            try:
                r.metadata['diagnostics_radial'] = r.diagnostics_radial.to_dict(orient='r')
            except AttributeError as ae:
                logging.error(ae)

            bulk_info.append(InsertOne(r.metadata))

    try:
        db.radials.bulk_write(bulk_info, ordered=False)
        logging.info('Bulk radial insert successful. {} radials inserted.'.format(len(bulk_info)))
    except BulkWriteError as bwe:
        logging.error(bwe.details)


if __name__ == '__main__':
    from glob import glob
    import datetime as dt
    import pandas as pd
    from datetime import datetime

    radial_dir = '/Volumes/boardwalk/codaradm/data/radials/'
    time_delta = 365  # days. Can be in fractions of a day.
    site_codes = ['BRAD', 'BRMR', 'BRNT', 'RATH', 'SEAB', 'SPRK', 'WOOD', 'ASSA', 'BLCK', 'BRIG', 'CEDR', 'SILD', 'PORT']  # Leave empty if you want the script to determine radial folders in root directory
    # site_codes = ['SILD']  # Leave empty if you want the script to determine radial folders in root directory
    sub_dirs = True  # If radials are in a subdirectory of the main radial folder, change to True
    sub_dir_format = '%Y_%m'
    load_file_limit = 1000
    multiprocessing = True
    workers = 16  # only applies to multiprocessing

    startTime = datetime.now()  # Grabbing this so we can determine script completion time at the end of the script

    if site_codes:
        # sites must be explicitly defined above
        paths = [os.path.join(radial_dir, x) for x in site_codes]
    else:
        # Do every site folder in radial_dir
        paths = [os.path.join(radial_dir, o) for o in os.listdir(radial_dir) if os.path.isdir(os.path.join(radial_dir, o))]

    ago = (dt.datetime.now() - dt.timedelta(days=time_delta))
    ago_folder = ago.replace(day=1, hour=0, minute=0, second=0, microsecond=0)
    ago_file = ago.replace(hour=0, minute=0, second=0, microsecond=0)

    if multiprocessing:
        import concurrent.futures

    for site_path in paths:
        if sub_dirs:
            dirs = [os.path.join(site_path, x) for x in pd.date_range(ago_folder, dt.datetime.utcnow(), freq='MS').strftime(sub_dir_format).tolist()]

            file_paths = []
            for item in dirs:
                file_paths = file_paths + glob(os.path.join(item, '*.ruv'))
        else:
            file_paths = glob(os.path.join(site_path, '*.ruv'))

        # Sort file paths
        file_paths = sorted(file_paths)

        # Chunk out lists into separate lists
        chunked_file_paths = divide_chunks(file_paths, load_file_limit)
        # chunked_file_paths = divide_chunks(['/Users/mikesmith/Downloads/RDLi_LOVE_2019_04_15_1200.ruv'], load_file_limit)

        if multiprocessing:
            # ~ 30 minutes for 12 months og long range data with single process
            # ~ 3.1 minutes for 12 months of long range data with ProcessPoolExecutor of 16 workers
            # ~ 13.0 minutes for 12 months of long range data with ThreadPoolExecutor of 16 workers
            with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                res = executor.map(main, chunked_file_paths)
        else:
            for chunk in chunked_file_paths:
                main(chunk)
            # main(chunked_file_paths)

    print(datetime.now() - startTime)
