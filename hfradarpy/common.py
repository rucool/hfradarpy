import datetime as dt
import glob
import os
import re
import pandas as pd

import logging

logger = logging.getLogger(__name__)


def divide_chunks(iterable, chunksize):
    """
    Yield successive n-sized chunks from a list

    Link: https://www.geeksforgeeks.org/break-list-chunks-size-n-python/

    Args:
        l (list): list to be broken down
        n (int): integer of the size chunks you want from the list

    Yields:
        list: a list containing lists of size, n
    """
    # looping till length l
    for i in range(0, len(iterable), chunksize):
        yield iterable[i: i + chunksize]


def list_files(types, main_dir, sub_directories=()):
    """
    Return a list of files given the directory of the files and extension type.
    You may also provide a list of sub_directories for the function to avoid.

    Args:
        types (str): file extension that you want to find
        main_dir (_type_): main directory that you want to recursively search for files
        sub_directories (tuple, optional):  Tuple containing strings of subdirectories you want to avoid. Defaults to ().

    Returns:
        list: list of files
    """
    file_list = []  # create empty list for finding files

    sub_dirs = [
        os.path.join(main_dir, o)
        for o in os.listdir(main_dir)
        if os.path.isdir(os.path.join(main_dir, o)) and o in sub_directories
    ]

    for sub in sub_dirs:
        for ext in types:
            file_list.extend(glob.glob(os.path.join(sub, ext)))
    file_list = sorted(file_list)
    return file_list


def list_to_dataframe(file_list):
    """
    Convert a list of ctf files, that are named in the the standard format 'year_month_day_hours' to a pandas dataframe

    Args:
        file_list (list): a list of files

    Returns:
        pd.DataFrame: Pandas DataFrame containing a list of files
    """
    df = pd.DataFrame(sorted(file_list), columns=["file"])
    try:
        df["time"] = df["file"].str.extract(r"(\d{4}_\d{2}_\d{2}_\d{4})")
        df["time"] = df["time"].apply(lambda x: dt.datetime.strptime(x, "%Y_%m_%d_%H%M"))
        df = df.set_index(["time"]).sort_index()
    except ValueError:
        logging.error("Cannot pass empty file_list to function. Returning empty dataframe.")
    return df


def timestamp_from_lluv_filename(filename):
    """
    Convert the string timestamp represented in CTF file names into a dt.datetime.

    Args:
        filename (str): filename

    Returns:
        dt.datetime: a datetime representation of the time included in the ctf filename
    """
    timestamp_regex = re.compile(r"\d{4}_\d{2}_\d{2}_\d{4}")
    mat_time = timestamp_regex.search(filename).group()
    timestamp = dt.datetime.strptime(mat_time, "%Y_%m_%d_%H%M")
    return timestamp
