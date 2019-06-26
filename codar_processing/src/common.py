import datetime as dt
import glob
import io
import logging
import os
import pprint
import re
from abc import ABCMeta, abstractmethod
from collections import OrderedDict

import numpy as np
import pandas as pd
import xarray as xr

logger = logging.getLogger(__name__)

desired_width = 320
pd.set_option('display.width', desired_width)
datetime_format = '%Y%m%dT%H%M%SZ'


def aggregate_netcdfs(files, save_dir, save_filename=None):
    """
    This function allows you to aggregate multiple netcdf files into a single file. It will concatenate on the coordinates of the netcdf files
    :param files: list of files or regular expression to location of files you want to open
    :param save_dir: directory to save aggregated netcdf file
    :param save_filename: filename of aggregated netcdf file
    :return:
    """
    # Create save directory if it doesn't exist
    create_dir(save_dir)

    print('Aggregating the following datasets:')
    pprint.pprint(files)

    # Opening files lazily (not into memory) using xarray
    ds = xr.open_mfdataset(files)
    ds.attrs['time_coverage_start'] = pd.Timestamp(ds['time'].min().data).strftime(datetime_format)
    ds.attrs['time_coverage_end'] = pd.Timestamp(ds['time'].max().data).strftime(datetime_format)

    # Encode variables for efficiency reasons
    encoding = make_encoding(ds)
    for v in list(ds.coords):
        encoding[v] = dict(zlib=False, _FillValue=False)

    print('Saving aggregated datasets as netCDF4 file.')
    ds.load()  # Load lazy arrays of open files into memory. Performance is better once loaded

    if save_filename:
        save_file = '{}/{}-{}_{}'.format(save_dir, ds.time_coverage_start, ds.time_coverage_end, save_filename)
        ds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims=['time'])
    else:
        save_file = '{}/{}-{}_totals_aggregated.nc'.format(save_dir, ds.time_coverage_start, ds.time_coverage_end)
        ds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims=['time'])
    ds.close()
    return save_file


# Yield successive n-sized chunks from l.
# Taken from https://www.geeksforgeeks.org/break-list-chunks-size-n-python/
def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]


def create_dir(new_dir):
    # Check if dir exists.. if it doesn't... create it.
    if not os.path.isdir(new_dir):
        try:
            os.makedirs(new_dir)
        except OSError:
            if os.path.exists(new_dir):
                pass
            else:
                raise


def list_files(types, main_dir, avoid_sub_directories):
    """

    :param types: file extension that you want to find
    :param main_dir: main directory that you want to recursively search for files
    :param avoid_sub_directories: Tuple containing strings of subdirectories you want to avoid
    :return:  file list
    """
    file_list = []  # create empty list for finding files

    sub_dirs = [os.path.join(main_dir, o) for o in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, o)) and o not in avoid_sub_directories]

    for sub in sub_dirs:
        for ext in types:
            file_list.extend(glob.glob(os.path.join(sub, ext)))
    file_list = sorted(file_list)
    return file_list


def list_to_dataframe(list):
    df = pd.DataFrame(list, columns=['file'])
    df['time'] = df['file'].str.extract(r'(\d{4}_\d{2}_\d{2}_\d{4})')
    df['time'] = df['time'].apply(lambda x: dt.datetime.strptime(x, '%Y_%m_%d_%H%M'))
    df = df.set_index(['time'])
    return df


def timestamp_from_lluv_filename(filename):
    timestamp_regex = re.compile('\d{4}_\d{2}_\d{2}_\d{4}')
    mat_time = timestamp_regex.search(filename).group()
    timestamp = dt.datetime.strptime(mat_time, '%Y_%m_%d_%H%M')
    return timestamp


def make_encoding(ds, time_start='days since 2006-01-01 00:00:00', comp_level=4, chunksize=10000, fillvalue=-999.00):
    encoding = {}

    for k in ds.data_vars:
        values = ds[k].values
        shape = values.shape

        encoding[k] = {'zlib': True, 'complevel': comp_level, '_FillValue': np.float32(fillvalue)}

        if 0 not in shape:
            if values.dtype.kind == 'O':
                values = values.astype('str')

            if values.dtype.kind == 'S':
                size = values.dtype.itemsize
                if size > 1:
                    shape = shape + (size,)

            dim0 = min(shape[0], chunksize)
            shape = (dim0,) + shape[1:]
            encoding[k]['chunksizes'] = shape

    # add the encoding for time so xarray exports the proper time
    encoding['time'] = dict(units=time_start, calendar='gregorian', zlib=False, _FillValue=False, dtype=np.double)
    # encoding['site_code_flags'] = dict(zlib=True, _FillValue=int(0))

    return encoding


class CTFParser(object):
    """
    A generic parser for the CODAR CTF file format.
    """

    __metaclass__ = ABCMeta

    def __init__(self, fname):
        """
        Return an LLUVParser object whose
        """
        split_path = os.path.split(fname)
        self.file_path = split_path[0]
        self.file_name = split_path[1]
        self.full_file = os.path.realpath(fname)
        self.metadata = OrderedDict()
        self._tables = OrderedDict()

        # Load the LLUV Data with this generic LLUV parsing routine below
        table_count = 0
        table = False  # Set table to False. Once a table is found, switch to True.
        self.is_wera = False  # Default false. If 'WERA' is detected in the Manufacturer flag. It is set to True
        processing_info = []

        with open(self.full_file, 'r', encoding='ISO-8859-1') as open_file:
            open_lluv = open_file.readlines()
            if any('%End:' in s for s in open_lluv):  # if there is no %End: the file is corrupt!
                # Parse header and footer metadata
                for line in open_lluv:
                    if not table:  # If we are not looking at a table or a tables header information
                        if line.startswith('%%'):
                            continue
                        elif line.startswith('%'):  # Parse the single commented header lines
                            key, value = self._parse_header_line(line)
                            if 'TableType' in line:  # Save this data as global header information
                                table = True  # we found a table
                                table_count = table_count + 1  # this is the nth table
                                table_data = u''
                                # self._data_header[table_count] = []
                                self._tables[str(table_count)] = OrderedDict()
                                self._tables[str(table_count)][key] = value
                                self._tables[str(table_count)]['_TableHeader'] = []
                            elif 'Manufacturer' in line:
                                if 'WERA' in value:
                                    self.is_wera = True
                            elif table_count > 0:
                                if key == 'ProcessingTool':
                                    processing_info.append(value)
                                else:
                                    self.metadata[key] = value
                            else:
                                self.metadata[key] = value
                    elif table:
                        if line.startswith('%'):
                            if line.startswith('%%'):  # table header information
                                self._tables[str(table_count)]['_TableHeader'].append(line)
                            else:  # Table metadata and diagnostic data are prepended by at least 1 % sign
                                if len(line.split(':')) == 1:  # Diagnostic Data
                                    line = line.replace('%', '').strip()
                                    table_data += '{}\n'.format(line)
                                else:  # Table data
                                    key, value = self._parse_header_line(line)
                                    if 'TableEnd' in line:
                                        # use pandas read_csv because it interprets the datatype for each column of the csv
                                        tdf = pd.read_csv(io.StringIO(table_data),
                                                          sep=' ',
                                                          header=None,
                                                          names=self._tables[str(table_count)]['TableColumnTypes'].split(),
                                                          skipinitialspace=True, )

                                        self._tables[str(table_count)]['data'] = tdf
                                        table = False
                                    else:
                                        self._tables[str(table_count)][key] = value
                        else:  # Uncommented lines are the main data table.
                            table_data += '{}'.format(line)
                self.metadata['ProcessingTool'] = processing_info
                if self.is_wera:
                    self._tables['1']['data'] = pd.read_csv(io.StringIO(table_data),
                                                            sep=' ',
                                                            header=None,
                                                            names=['LOND', 'LATD', 'VELU', 'VELV', 'VFLG', 'EACC', 'RNGE', 'BEAR', 'VELO', 'HEAD'],  # WERA has incorrect TableColumnTypes in their files.....
                                                            skipinitialspace=True, )
                self._iscorrupt = False
            else:
                logging.error('{}: File corrupt. Skipping to next file.'.format(self.full_file))
                self._iscorrupt = True

    def is_valid(self, table='1'):
        """
        Check if the data table for the file contains data
        :param table: string containing the table number to validate. Defaults to the primary data table '1'
        :return: True or False
        """
        try:
            return not self._tables[table]['data'].empty
        except:
            return False

    @staticmethod
    def _parse_header_line(line):
        """
        Parse a line into a key, value
        :param line: a line from a text file
        :type line: string
        :return: a tuple containing the key, value for the line
        :rtype: tuple
        """

        line = line.replace('%', '')  # Strip the % sign from the line
        line = line.replace('\n', '')  # Strip the new line character from the end of the line
        line_split = line.split(':')
        key = line_split[0]  # save key variable
        value = line_split[1].strip()  # save value variable and strip whitespace from value
        return key, value

    @abstractmethod
    def file_type(self):
        """Return a string representing the type of file this is."""
        pass

    def replace_invalid_values(self, values=[999.00, 1080.0]):
        """
        Convert invalid CODAR values to NaN
        :param df: dataframe
        :param values: list of CODAR fill values that reflect non calculable values
        :return: dataframe with invalid values set to NaN
        """
        logging.info('Replacing invalid values {} with NaN'.format(values))
        self.data.replace(values, np.nan, inplace=True)
