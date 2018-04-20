#!/usr/bin/env python
import glob
import logging
import sys
from codar_processing.common import aggregate_netcdfs

logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

save_dir = '/Users/mikesmith/Documents/2018_codar_reprocess/nc/2017_02/'
data_dir = '/Users/mikesmith/Documents/2018_codar_reprocess/nc/2017_02/'
save_filename = 'RU_MARA_2017_02_aggregated.nc'
regex = '*.nc'

data = sorted(glob.glob('{}/{}'.format(data_dir, regex)))

save_file = aggregate_netcdfs(data, save_dir, save_filename)

print('{} created'.format(save_file))
print('netCDF4 aggregation successful')