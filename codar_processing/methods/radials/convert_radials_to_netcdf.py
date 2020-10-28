#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR radial files utilizing the Radial subclass and convert to CF 1.6 NetCDF4 files
"""
import os
import glob

from codar_processing.src.radials import Radial


def main(radial_file, save_path, types=['tabular']):
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
            for t in types:
                r.export(
                    os.path.join(save_path, t, r.metadata['Site'].split()[0], r.file_name.strip('.ruv')),
                    'netcdf-{}'.format(t)
                )
        except ValueError:
            pass


if __name__ == '__main__':
    radial_path = '../../data/radials/ruv/SEAB/'
    radials = glob.glob(os.path.join(radial_path, '*.ruv'))
    save_path = '../../data/radials/nc/'
    types = ['tabular', 'multidimensional']

    for radial in radials:
        main(radial, save_path, types)
