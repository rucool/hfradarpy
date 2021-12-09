import datetime as dt
import glob
import os
import re
import pandas as pd

import logging
logger = logging.getLogger(__name__)


# Yield successive n-sized chunks from l.
# Taken from https://www.geeksforgeeks.org/break-list-chunks-size-n-python/
def divide_chunks(l, n):
    """
    Yield successive n-sized chunks from a list

    :param l: list to be broken down
    :param n: integer of the size chunks you want from the list
    """
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]


def create_dir(new_dir):
    """
    Helper function to create a new directory. It will try to make the directory even if it already exists, without overwriting the existing directory.
    """
    os.makedirs(new_dir, exist_ok=True)


def list_files(types, main_dir, sub_directories=()):
    """

    :param types: file extension that you want to find
    :param main_dir: main directory that you want to recursively search for files
    :param sub_directories: Tuple containing strings of subdirectories you want to avoid
    :return:  file list
    """
    file_list = []  # create empty list for finding files

    sub_dirs = [os.path.join(main_dir, o) for o in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, o)) and o in sub_directories]

    for sub in sub_dirs:
        for ext in types:
            file_list.extend(glob.glob(os.path.join(sub, ext)))
    file_list = sorted(file_list)
    return file_list


def list_to_dataframe(file_list):
    df = pd.DataFrame(sorted(file_list), columns=['file'])
    try:
        df['time'] = df['file'].str.extract(r'(\d{4}_\d{2}_\d{2}_\d{4})')
        df['time'] = df['time'].apply(lambda x: dt.datetime.strptime(x, '%Y_%m_%d_%H%M'))
        df = df.set_index(['time']).sort_index()
    except ValueError:
        logging.error('Cannot pass empty file_list to function. Returning empty dataframe.')
    return df


def timestamp_from_lluv_filename(filename):
    timestamp_regex = re.compile(r'\d{4}_\d{2}_\d{2}_\d{4}')
    mat_time = timestamp_regex.search(filename).group()
    timestamp = dt.datetime.strptime(mat_time, '%Y_%m_%d_%H%M')
    return timestamp
