#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR radial files utilizing the Radial subclass and convert to CF 1.6 NetCDF4 files
"""
import os
import glob

from codar_processing.src.radials import Radial


def main(radial_file, save_path):
    """
    Main function to parse and qc radial files
    :param radial_file: Path to radial file
    :param save_path: Path to save quality controlled radial file
    """
    try:
        r = Radial(radial_file)
    except Exception:
        return

    if r.is_valid():
        try:
            r.export(
                os.path.join(save_path, r.file_name.replace('.ruv', '.nc')),
                'netcdf'
            )
        except ValueError:
            pass


if __name__ == '__main__':
    radial_path = '../../data/radials/SEAB'
    radials = glob.glob(os.path.join(radial_path, '*.ruv'))
    save_path = '../../data/radials_nc/SEAB'

    for radial in radials:
        main(radial, save_path)
