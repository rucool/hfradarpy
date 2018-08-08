#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR radial files utilizing the Radial subclass and run class defined quality control (QC) methods
"""

import logging
import os
import sys
from codar_processing.radials import Radial

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
        r = Radial(radial_file, n_dimensional=True)
    except Exception as err:
        logging.error('{} - {}'.format(radial_file, err))
        return

    if r.is_valid():
        # run high frequency radar qartod tests on open radial file
        r.qc_qartod_location()
        r.qc_qartod_speed(qc_values['radial_max_speed'])
        r.qc_qartod_radial_count(qc_values['radial_min_count'], qc_values['radial_low_count'])

        # Export radial file to either a radial or netcdf
        try:
            r.export(os.path.join(save_path, r.file_name), 'radial')
        except ValueError as err:
            logging.error('{} - {}'.format(radial_file, err))
            pass


if __name__ == '__main__':
    radial = '../data/radials/SEAB/2018_03/RDLi_SEAB_2018_03_01_0200.ruv'
    save_path = '../data/radials_qc/SEAB/2018_03/'
    qc_values = dict(radial_max_speed=30, radial_min_count=50, radial_low_count=140)
    main(radial, save_path, qc_values)