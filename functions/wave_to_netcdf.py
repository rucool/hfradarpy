#!/usr/bin/env python
import logging
import numpy as np
import os
import sys
from codar_processing.calc import reckon
from configs.configs import netcdf_global_attributes
from codar_processing.common import make_encoding, create_dir
from codar_processing.waves import Waves

# time_string = '2018-03-01T00:00:00Z'

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

required_attributes = dict(title='MARACOOS Wave Heights',
                           naming_authority='edu.rutgers.marine.rucool',
                           comment='Network maintained by MARACOOS. See references attribute',
                           acknowledgment='This data is provided by the Mid-Atlantic Regional Association Coastal Ocean Observing System (MARACOOS). Funding is provided by the U.S. Integration Ocean Observing System (IOOS).',
                           standard_name_vocabulary='CF Standard Name Table v41',
                           creator_name='Michael Smith',
                           creator_email='michaesm@marine.rutgers.edu',
                           creator_url='rucool.marine.rutgers.edu',
                           institution='Center for Ocean Observing and Leadership, Department of Marine & Coastal Sciences, Rutgers University',
                           project='Mid-Atlantic Regional Association Coastal Ocean Observing System - High Frequency Radar Sea Surface Current Mapping',
                           sea_name='Mid-Atlantic Bight',
                           creator_type='person',
                           creator_institution='Rutgers University',
                           contributor_name='Scott Glenn, Josh Kohut, Hugh Roarty, Ethan Handel, Michael Smith, Laura Nazzaro, Teresa Updyke, Larry Atkinson, Rich Arena, Wendell Brown, Mike Muglia, Harvey Seim',
                           contributor_role='Principal Investigator, Principal Investigator, Principal Investigator, Hardware Maintenance, Data Manager, Data Manager, Hardware Maintenance, Principal Investigator, Hardware Maintenance, Principal Investigator, Hardware Maintenance, Principal Investigator',
                           platform='MARACOOS HF Radar 5MHz Network',
                           instrument='CODAR SeaSonde High Frequency Radar',
                           references='http://maracoos.org/node/146 https://rucool.marine.rutgers.edu/facilities https://rucool.marine.rutgers.edu/data',
                           summary='Ocean Wave Heights',
                           history='Half-hourly codar wave data from one site into one daily file containing wave heights.',
                           cdm_data_type='Timeseries',
                           source='CODAR SeaSonde Surface Current Mapping Device',
                           processing_level='Level 2',
                           keywords='Environmental Advisories > Marine Advisories > Marine Weather/Forecast, Oceans > Coastal Processes, Oceans > Ocean Circulation, Oceans > Ocean Waves, Oceans > Ocean Winds, Oceans > Ocean Tides, Spectral/Engineering > Radar',
                           publisher_name='Center for Ocean Observing and Leadership, Department of Marine & Coastal Sciences, Rutgers University',
                           publisher_email='michaesm@marine.rutgers.edu',
                           publisher_url='rucool.marine.rutgers.edu')


def main(wave_file, save_dir):
    """
    Function to parse the wave file
    :param fname: directory/filename of the wave file
    :return: Nothing
    """
    w = Waves(wave_file)
    old_style = False

    # Clean up wave header
    w.clean_wave_header()

    if w.data['DIST'].isnull().all():
        w.data = w.data.drop(['TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC', 'ACNT', 'RCLL', 'WDPT', 'MTHD', 'FLAG', 'WHNM', 'WHSD', 'PMWH', 'DIST'], axis=1)
    else:
        old_style = True
        w.data = w.data.drop(['TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC', 'ACNT', 'RCLL', 'WDPT', 'MTHD', 'FLAG', 'WHNM', 'WHSD', 'PMWH'], axis=1)

    w.data = w.data.set_index('datetime')

    # Convert pandas dataframe to xarray dataset
    ds = w.data.to_xarray()

    # if old_style:
    #     # rename variables to something meaningful
    #     rename = dict(datetime='time',
    #                   DIST='distance_from_origin',
    #                   MWHT='wave_height',
    #                   MWPD='wave_period',
    #                   WAVB='wave_bearing',
    #                   WNDB='wind_bearing')
    # else:
    #     # rename variables to something meaningful
    rename = dict(datetime='time',
                  MWHT='wave_height',
                  MWPD='wave_period',
                  WAVB='wave_bearing',
                  WNDB='wind_bearing')
    ds.rename(rename, inplace=True)

    mean_bearing = ds['wave_bearing'].where((ds['wave_bearing'] != 1080) & (ds['wave_bearing'] != 999)).mean().data.tolist()

    # Grab min and max time in dataset for entry into global attributes for cf compliance
    time_start = ds['time'].min().data
    time_end = ds['time'].max().data

    # Assign global attributes for CF compliant time series files
    global_attr = netcdf_global_attributes(ds, time_start, time_end)
    global_attr['Site'] = w.header['Site']
    ds = ds.assign_attrs(global_attr)

    # if old_style:
    #     ds['distance_from_origin'].attrs['long_name'] = 'distance of result from origin (transmitter) in kilometers'
    #     ds['distance_from_origin'].attrs['standard_name'] = 'sea_surface_wave_significant_height'
    #     ds['distance_from_origin'].attrs['units'] = 'km'
    #     ds['distance_from_origin'].attrs['comment'] = 'distance of result from origin in kilometers'
    #     ds['distance_from_origin'].attrs['valid_min'] = np.float32(0)
    #     ds['distance_from_origin'].attrs['valid_max'] = np.float32(100)

    # Set wave_height attributes
    ds['wave_height'].attrs['long_name'] = 'wave model height in meters'
    ds['wave_height'].attrs['standard_name'] = 'sea_surface_wave_significant_height'
    ds['wave_height'].attrs['units'] = 'm'
    ds['wave_height'].attrs['comment'] = 'wave model height in meters for every one of three waves'
    ds['wave_height'].attrs['valid_min'] = np.float32(0)
    ds['wave_height'].attrs['valid_max'] = np.float32(100)

    # Set wave_period attributes
    ds['wave_period'].attrs['long_name'] = 'wave spectra period in seconds'
    ds['wave_period'].attrs['standard_name'] = 'sea_surface_wave_mean_period'
    ds['wave_period'].attrs['units'] = 's'
    ds['wave_period'].attrs['comment'] = 'wave spectra period in seconds'
    ds['wave_period'].attrs['valid_min'] = np.float32(0)
    ds['wave_period'].attrs['valid_max'] = np.float32(100)

    # Set wave_bearing attributes
    ds['wave_bearing'].attrs['long_name'] = 'wave from direction in degrees'
    ds['wave_bearing'].attrs['standard_name'] = 'sea_surface_wave_from_direction'
    ds['wave_bearing'].attrs['units'] = 'degrees'
    ds['wave_bearing'].attrs['comment'] = 'wave from direction in degrees'
    ds['wave_bearing'].attrs['valid_min'] = np.float32(0)
    ds['wave_bearing'].attrs['valid_max'] = np.float32(360)

    # Set wind_bearing attributes
    ds['wind_bearing'].attrs['long_name'] = 'wind from direction in degrees'
    ds['wind_bearing'].attrs['standard_name'] = 'sea_surface_wind_wave_from_direction'
    ds['wind_bearing'].attrs['units'] = 'degrees'
    ds['wind_bearing'].attrs['comment'] = 'wind from direction in degrees'
    ds['wind_bearing'].attrs['valid_min'] = np.float32(0)
    ds['wind_bearing'].attrs['valid_max'] = np.float32(360)

    encoding = make_encoding(ds)

    create_dir(save_dir)
    nc_file = '{}.nc'.format(os.path.join(save_dir, os.path.basename(wave_file).split('.')[0]))

    # Convert files to netcdf
    ds.to_netcdf(nc_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims=['time'])


if __name__ == '__main__':
    wave_file = '../data/waves/wls/SEAB/WVLM_SEAB_2018_01_01_0000.wls'
    save_dir = '../data/waves/nc/SEAB'

    # main(wave_file, save_dir)
    try:
        exit(main(wave_file, save_dir))
    except Exception:
        logger.exception('Exception in main(): ')
        exit(1)
#
