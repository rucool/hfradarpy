#!/usr/bin/env python
import glob
import logging
import os
import sys
from codar_processing.common import aggregate_netcdfs

logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

save_dir = '/Users/mikesmith/Documents/2018_codar_reprocess/nc_agg/'
data_dir = '/Users/mikesmith/Documents/2018_codar_reprocess/lsq/'
save_filename = 'RU_MARA_{}_lsq_aggregated.nc'

sub_dirs = ['2017_01', '2017_02', '2017_03', '2017_04', '2017_05', '2017_06']

for sub in sub_dirs:
    path = os.path.join(data_dir, sub, '*.nc')
    save_name = save_filename.format(sub)
    save_file = aggregate_netcdfs(path, save_dir, save_name)

    print('{} created'.format(save_file))
    print('netCDF4 aggregation successful')