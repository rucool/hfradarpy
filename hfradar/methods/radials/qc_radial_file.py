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
import datetime as dt
from hfradar.src.radials import Radial

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


def main(radial_file, save_path, qc_values, export_type='radial'):
    """
    Main function to parse and qc radial files
    :param radial_file: Path to radial file
    :param save_path: Path to save quality controlled radial file
    :param qc_values: Dictionary containing thresholds for each QC test
    """
    try:
        r = Radial(radial_file, mask_over_land=False)
    except Exception as err:
        logging.error('{} - {}'.format(radial_file, err))
        return

    if r.is_valid():
        t0 = r.time - dt.timedelta(hours=1)
        previous_radial = '{}_{}'.format('_'.join(r.file_name.split('_')[:2]), t0.strftime('%Y_%m_%d_%H00.ruv'))
        previous_full_file = os.path.join(os.path.dirname(r.full_file), previous_radial)
        qc_keys = qc_values.keys()

        # run high frequency radar qartod tests on open radial file
        r.initialize_qc()
        r.qc_qartod_syntax()

        if 'qc_qartod_maximum_velocity' in qc_keys:
            r.qc_qartod_maximum_velocity(**qc_values['qc_qartod_maximum_velocity'])

        r.qc_qartod_valid_location()

        if 'qc_qartod_radial_count' in qc_keys:
            r.qc_qartod_radial_count(**qc_values['qc_qartod_radial_count'])

        if 'qc_qartod_spatial_median' in qc_keys:
            r.qc_qartod_spatial_median(**qc_values['qc_qartod_spatial_median'])

        if 'qc_qartod_temporal_gradient' in qc_keys:
            r.qc_qartod_temporal_gradient(previous_full_file)

        if 'qc_qartod_avg_radial_bearing' in qc_keys:
            r.qc_qartod_avg_radial_bearing(**qc_values['qc_qartod_avg_radial_bearing'])

        if 'qc_qartod_primary_flag' in qc_keys:
            r.qc_qartod_primary_flag(**qc_values['qc_qartod_primary_flag'])

        # Export radial file to either a radial or netcdf
        try:
            r.export(os.path.join(save_path, r.file_name), export_type)
        except ValueError as err:
            logging.error('{} - {}'.format(radial_file, err))
            pass


if __name__ == '__main__':
    radial_path = '../../data/radials/ruv/SEAB/'
    radials = glob.glob(os.path.join(radial_path, '*.ruv'))
    save_path = '../../data/radials_qc/ruv/SEAB/'
    export_type = 'radial'

    qc_values = dict(
        qc_qartod_avg_radial_bearing=dict(reference_bearing=151, warning_threshold=15, failure_threshold=30),
        qc_qartod_radial_count=dict(radial_min_count=75.0, radial_low_count=225.0),
        qc_qartod_maximum_velocity=dict(radial_max_speed=300.0, radial_high_speed=100.0),
        qc_qartod_spatial_median=dict(radial_smed_range_cell_limit=2.1, radial_smed_angular_limit=10,
                                      radial_smed_current_difference=30),
        qc_qartod_temporal_gradient=dict(gradient_temp_fail=32, gradient_temp_warn=25),
        qc_qartod_primary_flag=dict(include=['qc_qartod_syntax', 'qc_qartod_valid_location', 'qc_qartod_radial_count',
                                             'qc_qartod_maximum_velocity', 'qc_qartod_spatial_median'])
    )

    radials = glob.glob(os.path.join(radial_path, '*.ruv'))

    for radial in sorted(radials):
        main(radial, save_path, qc_values, export_type)
