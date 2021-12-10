#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR radial files utilizing the Radial subclass and convert to CF 1.6 NetCDF4 files
"""
import os
import glob
from hfradarpy.radials import Radial
from pathlib import Path
import logging
import sys

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


def main(radial_file, save_dir, types=['tabular']):
    """
    Main function to parse and qc radial files
    :param radial_file: Path to radial file
    :param save_path: Path to save quality controlled radial file
    """
    save_dir = Path(save_dir)

    try:
        r = Radial(radial_file)
    except Exception:
        return

    if r.is_valid():
        sname = save_dir / r.file_name
        try:
            for t in types:
                r.export(sname, 'netcdf-{}'.format(t), prepend_ext=True)
        except ValueError:
            pass


if __name__ == '__main__':
    data_root = (Path(__file__).parent.with_name('examples') / 'data').resolve()
    output_path = (Path(__file__).parent.with_name('examples') / 'output').resolve()

    radial_path = data_root / 'radials' / 'ruv' / 'SEAB'
    radial_output = output_path / 'radials' / 'nc' / 'SEAB'

    types = ['tabular', 'multidimensional']

    for radial in sorted(radial_path.glob('*.ruv')):
        try:
            print(str(radial))
            main(radial, radial_output, types)
        except Exception:
            logger.exception('Exception in main(): ')
