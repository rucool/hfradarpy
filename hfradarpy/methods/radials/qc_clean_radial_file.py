#!/usr/bin/env python
"""
@author Teresa Updyke
@email garner@ccpo.odu.edu
@purpose Remove radial data from a CODAR file if quality control tests resulted in a fail code in the primary flag
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


def main(radial_file, radial_qc_file, save_path, export_type='radial'):
    """
    Main function to clean radial files
    :param radial_file: Path to radial file
    :param save_path: Path to save cleaned radial file
    """
    try:
        r = Radial(radial_file, mask_over_land=False)
    except Exception as err:
        logging.error('{} - {}'.format(radial_file, err))
        return
    try:
        rqc = Radial(radial_qc_file, mask_over_land=False)
    except Exception as err:
        logging.error('{} - {}'.format(radial_file, err))
        return

    if r.is_valid() & rqc.is_valid():   # this passes if there is data in the table, other checks?

        d = r.data
        dqc = rqc.data
        if 'PRIM' in rqc.data:
            rt = d[dqc['PRIM'] != 4]
            r.data = rt

            for key in r._tables.keys():
                table = r._tables[key]
                if 'LLUV' in table['TableType']:
                    r._tables['1']['TableRows'] = rt.shape[0]
            #else:
            #   warning that it didn't update number of table rows
        #else:
        # warning of failure to update file, the original will be exported

        # Export radial file to either a radial or netcdf
        try:
            r.export(os.path.join(save_path, r.file_name), export_type)
        except ValueError as err:
            logging.error('{} - {}'.format(radial_file, err))
            pass


if __name__ == '__main__':
    radial_path = '../../data/radials/ruv/SEAB/'
    qc_radial_path = '../../data/radials_qc/ruv/SEAB/'
    save_path = '../../data/radials_clean/ruv/SEAB/'
    export_type = 'radial'

    n = len(qc_radial_path)
    qc_radials = glob.glob(os.path.join(qc_radial_path, '*.ruv'))

    for qc_radial in sorted(qc_radials):
        radial = radial_path + qc_radial[n:]
        main(radial, qc_radial, save_path, export_type)
