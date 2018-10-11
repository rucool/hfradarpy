#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Realtime parsing of CODAR wave files utilizing the Wave subclass and convert to CF Compliant NetCDF4 files
"""
import datetime as dt
import glob
import logging
import os
import sys
from functions.waves.wave_to_netcdf import main as wave_to_netcdf

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

# Script settings
wave_dir = '/home/codaradm/data/waves/'
save_dir = '/home/codaradm/data/waves/'
days_to_check = 30
sites = ['BRAD', 'BRMR', 'BRNT', 'RATH', 'SEAB', 'SPRK', 'WOOD', 'BISL', 'CMPT', 'CAPE', 'CPHN', 'GCAP', 'HLPN', 'MISQ',
         'MNTK', 'OLDB', 'PALM', 'PORT', 'SILD', 'STLI', 'SUNS', 'VIEW']

# Calculate datetime of time delta set in days_to_check from now
now = dt.datetime.now()
ago = now-dt.timedelta(days=days_to_check)

for site in sites:
    logging.debug('Checking {} wave archive for files modified in the past {} days'.format(site, days_to_check))
    site_dir = os.path.join(wave_dir, site)
    save_dir_site = os.path.join(save_dir, site, 'nc')

    for fname in sorted(glob.glob(os.path.join(site_dir, 'WVLM*.wls'))):
        st = os.stat(fname)
        mtime = dt.datetime.fromtimestamp(st.st_mtime)
        if mtime > ago:
            logging.info('{} modified {} days ago'.format(fname, mtime))
            wave_to_netcdf(fname, save_dir_site)