#!/usr/bin/env python
import os
from codar_processing.radials import Radial


def main(radial_file, save_path):
    r = Radial(radial_file)

    # run high frequencvy radar qartod tests on open radial file
    r.qc_qartod_location()
    r.qc_qartod_speed(250)
    r.qc_qartod_radial_count(150, 300)

    # Export radial file to either a radial or netcdf
    r.export(save_path, 'radial')


if __name__ == '__main__':
    radial = '../data/radials/SEAB/2018_03/RDLi_SEAB_2018_03_01_0200.ruv'
    save_path = '../data/radials_qc/SEAB/2018_03/{}'.format(os.path.basename(radial))
    main(radial, save_path)