#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR radial files utilizing the Radial subclass and run class defined quality control (QC) methods
"""

import logging
import os
import sys
import glob

from codar_processing.src.radials import Radial

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


def main(radial_file, save_path, qc_values):
    """
    Main function to parse and qc radial files
    :param radial_file: Path to radial file
    :param save_path: Path to save quality controlled radial file
    :param qc_values: Dictionary containing thresholds for each QC test
    """
    try:
        r = Radial(radial_file, to_xarray=False)
    except Exception as err:
        logging.error('{} - {}'.format(radial_file, err))
        return

    if r.is_valid():
        # run high frequency radar qartod tests on open radial file
        r.initialize_qc()
        r.qc_qartod_syntax()
        r.qc_qartod_maximum_velocity(**qc_values['qc_qartod_maximum_velocity'])
        r.qc_qartod_valid_location()
        r.qc_qartod_radial_count(**qc_values['qc_qartod_radial_count'])
        r.qc_qartod_spatial_median(**qc_values['qc_qartod_spatial_median'])
        # r.qc_qartod_avg_radial_bearing(qc_values['average_bearing_threshold'])

        # Export radial file to either a radial or netcdf
        try:
            r.export(os.path.join(save_path, r.file_name), 'radial')
        except ValueError as err:
            logging.error('{} - {}'.format(radial_file, err))
            pass


if __name__ == '__main__':
    radial_path = '../../data/radials/SEAB'
    radials = glob.glob(os.path.join(radial_path, '*.ruv'))
    save_path = '../../data/radials_qc/SEAB/'

    qc_values = dict(
        qc_qartod_radial_count=dict(radial_min_count=50, radial_low_count=140),
        qc_qartod_maximum_velocity=dict(radial_max_speed=300),
        qc_qartod_spatial_median=dict(radial_smed_range_cell_limit=2.1,
                                      radial_smed_angular_limit=10,
                                      radial_smed_current_difference=30))

    for radial in radials:
        main(radial, save_path, qc_values)
