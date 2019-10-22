#!/usr/bin/env python
"""
@file aggregate_netcdfs.py
@author Mike Smith
@email michaesm@marine.rutgers.edu
@brief This utility aggregates netCDF files in a given folder.
@purpose To convert many netCDF files into a single netCDF file that is aggregated on the same dimension
"""
import logging
import sys
from codar_processing.src.common import aggregate_netcdfs
from glob import glob

logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


def main(data_dir, save_dir, save_name):
    """

    :param data_dir: Directory that contains netCDF files.
    :param save_dir: Directory that you want to save aggregated netCDF file.
    :param save_name: Filename you want to name aggregated netCDF file. This will be prepended by start and end time of data. Leave blank for default filename.
    :return:
    """
    data = sorted(glob(data_dir))
    save_file = aggregate_netcdfs(data, save_dir, save_name)

    logging.info('{} created'.format(save_file))
    logging.info('netCDF4 aggregation successful')


if __name__ == '__main__':
    data_dir = '../data/totals/oi/nc/hourly/*.nc'
    save_dir = '../data/totals/oi/nc/aggregated/'
    save_name = 'aggregated.nc'
    main(data_dir, save_dir, save_name)
