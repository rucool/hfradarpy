#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Parse CODAR wave files utilizing the Wave subclass and convert to CF Compliant NetCDF4 files
"""
# import click
import logging
import numpy as np
import os
import sys
from codar_processing.calc import reckon
from configs.configs import netcdf_global_attributes
from codar_processing.common import make_encoding, create_dir
from codar_processing.waves import Waves

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

required_attributes = dict(ncei_template_version='NCEI_NetCDF_Point_Template_v2.0',
                           title='MARACOOS Wave Heights',
                           naming_authority='edu.rutgers.marine.rucool',
                           comment='Network maintained by MARACOOS. See references attribute',
                           acknowledgement='This data is provided by the Mid-Atlantic Regional Association Coastal Ocean Observing System (MARACOOS). Funding is provided by the U.S. Integration Ocean Observing System (IOOS).',
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
                           cdm_data_type='Point',
                           source='CODAR SeaSonde Surface Current Mapping Device',
                           processing_level='Level 2',
                           keywords='Environmental Advisories > Marine Advisories > Marine Weather/Forecast, Oceans > Coastal Processes, Oceans > Ocean Circulation, Oceans > Ocean Waves, Oceans > Ocean Winds, Oceans > Ocean Tides, Spectral/Engineering > Radar',
                           publisher_name='Center for Ocean Observing and Leadership, Department of Marine & Coastal Sciences, Rutgers University',
                           publisher_email='michaesm@marine.rutgers.edu',
                           publisher_url='rucool.marine.rutgers.edu')


# @click.command()
# @click.option('--wave_file', default='../data/waves/wls/SEAB/WVLM_SEAB_2018_01_01_0000.wls', help='Path to Wave file')
# @click.option('--save_dir', default='../data/waves/nc/SEAB', help='Path to save files')
def main(wave_file, save_dir):
    """
    Function to parse the wave file and output as NetCDF4 file
    :param wave_file: Path to wave file
    :param save_dir: Path to save directory for generated NetCDF4 files
    """
    w = Waves(wave_file)

    # Clean up wave header
    w.clean_wave_header()

    # Remove wave heights less than 0.2 and greater than 5 m
    w.remove_bad_wave_heights(0.2, 5)

    # Calculate point location
    radar_median = np.median([float(x) for x in w.header['CoastlineSector'].split(', ')])
    origin = [float(x) for x in w.header['Origin'].split('  ')]
    point = reckon(origin[0], origin[1], radar_median, 2)

    # create lat and lon columns
    w.data['lon'] = point[1]
    w.data['lat'] = point[0]

    if w.data['DIST'].isnull().all():
        w.data = w.data.drop(['TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC', 'ACNT', 'RCLL', 'WDPT', 'MTHD', 'FLAG', 'WHNM', 'WHSD', 'PMWH', 'DIST'], axis=1)
    else:
        # w.data = w.data.drop(['TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC', 'ACNT', 'RCLL', 'WDPT', 'MTHD', 'FLAG', 'WHNM', 'WHSD', 'PMWH'], axis=1)
        return

    w.data = w.data.set_index('datetime')

    # Convert pandas dataframe to xarray dataset
    ds = w.data.to_xarray()

    # rename variables to something meaningful
    rename = dict(datetime='time',
                  MWHT='wave_height',
                  MWPD='wave_period',
                  WAVB='wave_bearing',
                  WNDB='wind_bearing')

    ds.rename(rename, inplace=True)

    # Grab min and max time in dataset for entry into global attributes for cf compliance
    time_start = ds['time'].min().data
    time_end = ds['time'].max().data

    # Assign global attributes for CF compliant time series files
    global_attr = netcdf_global_attributes(required_attributes, time_start, time_end)
    global_attr['Site'] = w.header['Site']
    ds = ds.assign_attrs(global_attr)

    # set time attribute
    ds['time'].attrs['standard_name'] = 'time'
    ds['time'].attrs['long_name'] = 'Universal Time Coordinated (UTC) Time'

    # Set wave_height attributes
    ds['wave_height'].attrs['long_name'] = 'wave model height in meters'
    ds['wave_height'].attrs['standard_name'] = 'sea_surface_wave_significant_height'
    ds['wave_height'].attrs['units'] = 'm'
    ds['wave_height'].attrs['comment'] = 'wave model height in meters for every one of three waves'
    ds['wave_height'].attrs['valid_min'] = np.double(0)
    ds['wave_height'].attrs['valid_max'] = np.double(100)
    ds['wave_height'].attrs['coordinates'] = 'time'
    ds['wave_height'].attrs['grid_mapping'] = 'crs'
    ds['wave_height'].attrs['coverage_content_type'] = 'physicalMeasurement'

    # Set wave_period attributes
    ds['wave_period'].attrs['long_name'] = 'wave spectra period in seconds'
    ds['wave_period'].attrs['standard_name'] = 'sea_surface_wave_mean_period'
    ds['wave_period'].attrs['units'] = 's'
    ds['wave_period'].attrs['comment'] = 'wave spectra period in seconds'
    ds['wave_period'].attrs['valid_min'] = np.double(0)
    ds['wave_period'].attrs['valid_max'] = np.double(100)
    ds['wave_period'].attrs['coordinates'] = 'time'
    ds['wave_period'].attrs['grid_mapping'] = 'crs'
    ds['wave_period'].attrs['coverage_content_type'] = 'physicalMeasurement'

    # Set wave_bearing attributes
    ds['wave_bearing'].attrs['long_name'] = 'wave from direction in degrees'
    ds['wave_bearing'].attrs['standard_name'] = 'sea_surface_wave_from_direction'
    ds['wave_bearing'].attrs['units'] = 'degrees'
    ds['wave_bearing'].attrs['comment'] = 'wave from direction in degrees'
    ds['wave_bearing'].attrs['valid_min'] = np.double(0)
    ds['wave_bearing'].attrs['valid_max'] = np.double(360)
    ds['wave_bearing'].attrs['coordinates'] = 'time'
    ds['wave_bearing'].attrs['grid_mapping'] = 'crs'
    ds['wave_bearing'].attrs['coverage_content_type'] = 'physicalMeasurement'

    # Set wind_bearing attributes
    ds['wind_bearing'].attrs['long_name'] = 'wind from direction in degrees'
    ds['wind_bearing'].attrs['standard_name'] = 'sea_surface_wind_wave_from_direction'
    ds['wind_bearing'].attrs['units'] = 'degrees'
    ds['wind_bearing'].attrs['comment'] = 'wind from direction in degrees'
    ds['wind_bearing'].attrs['valid_min'] = np.double(0)
    ds['wind_bearing'].attrs['valid_max'] = np.double(360)
    ds['wind_bearing'].attrs['coordinates'] = 'time'
    ds['wind_bearing'].attrs['grid_mapping'] = 'crs'
    ds['wind_bearing'].attrs['coverage_content_type'] = 'physicalMeasurement'

    # Set lon attributes
    ds['lon'].attrs['long_name'] = 'Longitude'
    ds['lon'].attrs['standard_name'] = 'longitude'
    ds['lon'].attrs['short_name'] = 'lon'
    ds['lon'].attrs['units'] = 'degrees_east'
    ds['lon'].attrs['axis'] = 'X'
    ds['lon'].attrs['valid_min'] = np.double(-180.0)
    ds['lon'].attrs['valid_max'] = np.double(180.0)
    ds['lon'].attrs['grid_mapping'] = 'crs'

    # Set lat attributes
    ds['lat'].attrs['long_name'] = 'Latitude'
    ds['lat'].attrs['standard_name'] = 'latitude'
    ds['lat'].attrs['short_name'] = 'lat'
    ds['lat'].attrs['units'] = 'degrees_north'
    ds['lat'].attrs['axis'] = 'Y'
    ds['lat'].attrs['valid_min'] = np.double(-90.0)
    ds['lat'].attrs['valid_max'] = np.double(90.0)
    ds['lat'].attrs['grid_mapping'] = 'crs'

    encoding = make_encoding(ds)

    # add container variables that contain no data
    kwargs = dict(crs=None, instrument=None)
    ds = ds.assign(**kwargs)

    # Set crs attributes
    ds['crs'].attrs['grid_mapping_name'] = 'latitude_longitude'
    ds['crs'].attrs['inverse_flattening'] = 298.257223563
    ds['crs'].attrs['long_name'] = 'Coordinate Reference System'
    ds['crs'].attrs['semi_major_axis'] = '6378137.0'
    ds['crs'].attrs['epsg_code'] = 'EPSG:4326'
    ds['crs'].attrs['comment'] = 'http://www.opengis.net/def/crs/EPSG/0/4326'

    ds['instrument'].attrs['long_name'] = 'CODAR SeaSonde High Frequency Radar'
    ds['instrument'].attrs['sensor_type'] = 'Direction-finding high frequency radar antenna'
    ds['instrument'].attrs['make_model'] = 'CODAR SeaSonde'
    ds['instrument'].attrs['serial_number'] = 1

    create_dir(save_dir)
    nc_file = '{}.nc'.format(os.path.join(save_dir, os.path.basename(wave_file).split('.')[0]))

    # Convert files to netcdf
    ds.to_netcdf(nc_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims=['time'])


if __name__ == '__main__':
    wave_file = '../data/waves/wls/SEAB/WVLM_SEAB_2018_01_01_0000.wls'
    save_dir = '../data/waves/nc/SEAB'
    try:
        main(wave_file, save_dir)
    except Exception:
        logger.exception('Exception in main(): ')
        exit(1)