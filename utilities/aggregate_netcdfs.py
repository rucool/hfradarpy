#!/usr/bin/env python
"""
@file aggregate_netcdfs.py
@author Mike Smith
@email michaesm@marine.rutgers.edu
@brief This utility aggregates netCDF files in a given folder.
@purpose To convert many netCDF files into a single netCDF file that is aggregated on the same dimension
"""
import click
import logging
import sys
from codar_processing.common import aggregate_netcdfs
from glob import glob

logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


@click.command()
@click.option('--data_dir', default='../data/totals/nc/*.nc', help='Path to netCDF files')
@click.option('--save_dir', default='../data/totals/nc/monthly/', help='Path to save files')
@click.option('--save_name', default='aggregated.nc', help='Save filename')
def main(data_dir, save_dir, save_name):
    data = sorted(glob(data_dir))
    save_file = aggregate_netcdfs(data, save_dir, save_name)

    logging.info('{} created'.format(save_file))
    logging.info('netCDF4 aggregation successful')


if __name__ == '__main__':
    main()
