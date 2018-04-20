#!/usr/bin/env python
import glob
import logging
import sys
from codar_processing.common import aggregate_netcdfs

logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

save_dir = '../data/totals/nc/monthly/'
data_dir = '../data/totals/nc/'
save_filename = 'RU_MARA_2017_01_aggregated_filtered.nc'
regex = '*.nc'

data = sorted(glob.glob('{}/{}'.format(data_dir, regex)))

save_file = aggregate_netcdfs(data, save_dir, save_filename)

print('{} created'.format(save_file))
print('netCDF4 aggregation successful')