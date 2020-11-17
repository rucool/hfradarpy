#!/usr/bin/env python
"""
Convert CODAR Totals MATLAB .mat files generated using the HFRProgs toolbox into Climate Forecasting compliant netcdf files
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Convert MAT files created using the hfrProgs MATLAB toolbox into CF-1.6/NCEI Grid 2.0 compliant netCDF4 files
"""
import datetime as dt
import logging
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr
from numpy import matlib as mb
from scipy.io import loadmat

from hfradar.src.calc import gridded_index
from hfradar.src.common import create_dir, timestamp_from_lluv_filename
from hfradar.src.common import make_encoding
from hfradar.configs import configs_default as configs

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


def matlab2datetime(matlab_time):
    """
    Convert Matlab time to Python datetime
    :param matlab_time: MATLAB datenum integer
    :return: Python Datetime
    """
    day = dt.datetime.fromordinal(int(matlab_time))
    day_frac = dt.timedelta(days=matlab_time % 1) - dt.timedelta(days=366)
    return day + day_frac


def main(grid, mat_file, save_dir, user_attributes, flags=None, domain=[], method='oi'):
    """
    Convert MAT files created using the hfrProgs MATLAB toolbox into CF-1.6/NCEI Grid 2.0 compliant netCDF4 files
    :param grid: CSV file containing lon,lat grid information
    :param mat_file: Filepath to MAT file containing HFRProgs
    :param save_dir: Directory to save netCDF files to
    :param user_attributes: User defined dataset attributes for netCDF global attribute. Required for CF/NCEI compliance
    :param flags: Dictionary of thresholds at which we should filter data above
    :param method: 'oi' or 'lsq'. OI is optimal interpolation. LSQ is unweighted least squares
    """
    fname = os.path.basename(mat_file)
    try:
        # load .mat file
        data = loadmat(mat_file, squeeze_me=True, struct_as_record=False)
        logging.debug('{} - MAT file successfully loaded '.format(fname))
    except Exception as err:
        logging.error('{} - {}. MAT file could not be loaded.'.format(fname, err))
        #return

    if not domain:
        domain = data['TUV'].DomainName
        if not domain:
            domain = 'MARA'
    else:
        domain = 'MARA'

    time = timestamp_from_lluv_filename(mat_file)

    # convert matlab time to python datetime
    # time = dt.datetime.strptime(mat_time, '%Y_%m_%d_%H00')
    time_index = pd.date_range(time.strftime('%Y-%m-%d %H:%M:%S'), periods=1)  # create pandas datetimeindex from time
    time_string = time.strftime('%Y%m%dT%H%M%SZ')  # create timestring from time

    file_name = 'RU_{}_{}.nc'.format(domain, time_string)
    file_and_path = os.path.join(save_dir, file_name)

    try:
        logging.debug('{} - Saving file data to variables'.format(fname))
        # load longitude and latitude data associated with variables
        lonlat = data['TUV'].LonLat.astype(np.float32)

        # create variables for eastward and northward velocities
        u = data['TUV'].U.astype(np.float32)
        v = data['TUV'].V.astype(np.float32)
        u_units = data['TUV'].UUnits
        v_units = data['TUV'].VUnits

        #maxspd = data['TUV'].OtherMetadata.cleanTotals.maxspd
        maxspd = data['conf'].Totals.MaxTotSpeed

        if method == 'oi':
            # create variables for associated error values
            u_err = data['TUV'].ErrorEstimates.Uerr.astype(np.float32)
            v_err = data['TUV'].ErrorEstimates.Verr.astype(np.float32)
            uv_covariance = data['TUV'].ErrorEstimates.UVCovariance.astype(np.float32)
            total_errors = data['TUV'].ErrorEstimates[1].TotalErrors.astype(np.float32)

            # Data Processing Information
            num_rads = data['TUV'].OtherMatrixVars.makeTotalsOI_TotalsNumRads.astype(int)
            min_rads = data['TUV'].OtherMetadata.makeTotalsOI.parameters.MinNumRads
            min_sites = data['TUV'].OtherMetadata.makeTotalsOI.parameters.MinNumSites
            mdlvar = data['TUV'].OtherMetadata.makeTotalsOI.parameters.mdlvar
            errvar = data['TUV'].OtherMetadata.makeTotalsOI.parameters.errvar
            sx = data['TUV'].OtherMetadata.makeTotalsOI.parameters.sx
            sy = data['TUV'].OtherMetadata.makeTotalsOI.parameters.sy
            temporal_threshold = data['TUV'].OtherMetadata.makeTotalsOI.parameters.tempthresh
            #processing_parameters = [maxspd, min_sites, min_rads, temporal_threshold, sx, sy, mdlvar, errvar]
            processing_parameters = [min_sites, min_rads, temporal_threshold, sx, sy, mdlvar, errvar]
            #processing_parameters_info = '1) Maximum Total Speed Threshold (cm s-1)\n'
            processing_parameters_info = '1) Minimum number of radial sites\n'
            processing_parameters_info += '2) Minimum number of radial vectors\n'
            processing_parameters_info += '3) Temporal search window for radial solutions (Fraction of a day)\n'
            processing_parameters_info += '4) Decorrelation scales in the north direction\n'
            processing_parameters_info += '5) Decorrelation scales in the east direction\n'
            processing_parameters_info += '6) Signal variance of the surface current fields (cm2 s-2)\n'
            processing_parameters_info += '7) Data error variance of the input radial velocities (cm2 s-2)\n'

            #QC Information
            uerr_testname = data['conf'].Totals.OI.cleanTotalsVarargin[0][0] + ' ' + data['conf'].Totals.OI.cleanTotalsVarargin[0][1]
            uerr_threshold = data['conf'].Totals.OI.cleanTotalsVarargin[0][2]
            verr_testname = data['conf'].Totals.OI.cleanTotalsVarargin[1][0] + ' ' + data['conf'].Totals.OI.cleanTotalsVarargin[1][1]
            verr_threshold = data['conf'].Totals.OI.cleanTotalsVarargin[1][2]
            qc_info = 'Quality control reference: IOOS QARTOD HF Radar ver 1.0 May 2016\n'
            qc_info += 'QCFlagDefinitions: 1 = pass, 2 = not_evaluated, 3 = suspect, 4 = fail, 9 = missing_data\n'
            qc_primary_flag_info = 'QCPrimaryFlagDefinition: Highest flag value of QC16, QC18, QC19, QC20\n'
            qc_primary_flag_info += 'This flag will be set to not_evaluated only if ALL individual tests were not_evaluated.'
            qc_operator_mask_info = qc_info + 'The qc_operator_mask follows QCFlagDefinitions and is set at discretion of the operator or data manager.'
            qc16_info = qc_info + 'QC16 Max Speed Threshold [max_vel = ' + str(maxspd) + ' (cm/s)]'
            qc18_info = qc_info + 'QC18 Valid Location [landmask file = ' + data['conf'].Totals.MaskFile + ']'
            qc19_info = qc_info + 'QC19 OI Uncertainty Threshold [' + uerr_testname + ' ' + str(uerr_threshold) + ']'
            qc20_info = qc_info + 'QC20 OI Uncertainty Threshold [' + verr_testname + ' ' + str(verr_threshold) + ']'
            qc16 = data['TUVqc'].QC16.astype(np.int32)
            qc18 = data['TUVqc'].QC18.astype(np.int32)
            qc19 = data['TUVqc'].QC19.astype(np.int32)
            qc20 = data['TUVqc'].QC20.astype(np.int32)
            qc_primary_flag = data['TUVqc'].PRIM.astype(np.int32)
            qc_operator_mask = data['TUVqc'].qc_operator_mask.astype(np.int32)


        elif method == 'lsq':
            # create variables for associated error values
            u_err = data['TUV'].ErrorEstimates[1].Uerr.astype(np.float32)
            v_err = data['TUV'].ErrorEstimates[1].Verr.astype(np.float32)
            uv_covariance = data['TUV'].ErrorEstimates[1].UVCovariance.astype(np.float32)
            total_errors = data['TUV'].ErrorEstimates[1].TotalErrors.astype(np.float32)

            # Data Processing Information
            num_rads = data['TUV'].OtherMatrixVars.makeTotals_TotalsNumRads.astype(int)
            min_rads = data['TUV'].OtherMetadata.makeTotals.parameters.MinNumRads
            min_sites = data['TUV'].OtherMetadata.makeTotals.parameters.MinNumSites
            spatial_threshold = data['TUV'].OtherMetadata.makeTotals.parameters.spatthresh
            temporal_threshold = data['TUV'].OtherMetadata.makeTotals.parameters.tempthresh
            #processing_parameters = [maxspd, min_sites, min_rads, temporal_threshold, spatial_threshold]
            processing_parameters = [min_sites, min_rads, temporal_threshold, spatial_threshold]
            #processing_parameters_info = '1) Maximum Total Speed Threshold (cm s-1)\n'
            processing_parameters_info = '1) Minimum number of radial sites.\n'
            processing_parameters_info += '2) Minimum number of radial vectors.\n'
            processing_parameters_info += '3) Temporal search window for radial solutions (Fractions of a day)\n'
            processing_parameters_info += '4) Spatial search radius for radial solutions (km)\n'


            #QC Information
            gdoptestname = data['conf'].Totals.cleanTotalsVarargin[0] + ' ' + data['conf'].Totals.cleanTotalsVarargin[1]
            gdopthreshold = data['conf'].Totals.cleanTotalsVarargin[2]
            qc_info = 'Quality control reference: IOOS QARTOD HF Radar ver 1.0 May 2016\n'
            qc_info += 'QCFlagDefinitions: 1 = pass, 2 = not_evaluated, 3 = suspect, 4 = fail, 9 = missing_data\n'
            qc_primary_flag_info = 'QCPrimaryFlagDefinition: Highest flag value of QC16, QC18, QC19, QC20\n'
            qc_primary_flag_info += 'This flag will be set to not_evaluated only if ALL tests were not_evaluated.'
            qc_operator_mask_info = qc_info + 'The qc_operator_mask follows QCFlagDefinitions and is set at discretion of the operator or data manager.'
            qc16_info = qc_info + 'QC16 Max Speed Threshold [max_vel = ' + str(maxspd) + ' (cm/s)]'
            qc18_info = qc_info + 'QC18 Valid Location [landmask file = ' + data['conf'].Totals.MaskFile + ']'
            qc15_info = qc_info + 'QC15 GDOP Threshold [' + gdoptestname + ' ' + str(gdopthreshold) + ']'
            qc15 = data['TUVqc'].QC15.astype(np.int)
            qc16 = data['TUVqc'].QC16.astype(np.int)
            qc18 = data['TUVqc'].QC18.astype(np.int)
            qc_primary_flag = data['TUVqc'].PRIM.astype(np.int)
            qc_operator_mask = data['TUVqc'].qc_operator_mask.astype(np.int)

    except AttributeError as err:
        logging.error('{} - {}. MAT file missing variable needed to create netCDF4 file'.format(fname, err))
        #return

    # Create a grid to shape 1d data
    lon = np.unique(grid['lon'].values.astype(np.float32))
    lat = np.unique(grid['lat'].values.astype(np.float32))
    [x, y] = np.meshgrid(lon, lat)

    # Create a dictionary of variables that we want to grid
    if method == 'oi':
        data_dict = dict(u=u,
                         v=v,
                         u_err=u_err,
                         v_err=v_err,
                         uv_covariance=uv_covariance,
                         total_errors=total_errors,
                         num_radials=num_rads,
                         qc16_maxspeed=qc16,
                         qc18_validlocation=qc18,
                         qc19_uerr=qc19,
                         qc20_verr=qc20,
                         qc_primary_flag=qc_primary_flag,
                         qc_operator_mask=qc_operator_mask,
                         )
    elif method == 'lsq':
        data_dict = dict(u=u,
                         v=v,
                         u_err=u_err,
                         v_err=v_err,
                         uv_covariance=uv_covariance,
                         total_errors = total_errors,
                         num_radials=num_rads,
                         qc15_total_errors=qc15,
                         qc16_maxspeed=qc16,
                         qc18_validlocation=qc18,
                         qc_primary_flag=qc_primary_flag,
                         qc_operator_mask=qc_operator_mask,
                         )

    logging.debug('{} - Gridding data to 2d grid'.format(fname))

    # convert 1d data into 2d gridded form. data_dict must be a dictionary.
    x_ind, y_ind = gridded_index(x, y, lonlat[:, 0], lonlat[:, 1])

    for key in data_dict.keys():
        temp_data = mb.tile(np.nan, x.shape)
        temp_data[(y_ind, x_ind)] = data_dict[key]

        # expand dimensions for time and depth
        count = 0
        while count < 2:  # add two dimensions to from of array for time and z (depth)
            temp_data = np.expand_dims(temp_data, axis=0)
            count = count + 1
            data_dict[key] = temp_data

    logging.debug('{} - Loading data into xarray dataset'.format(fname))

    # initialize xarray dataset. Add variables. Add coordinates
    ds = xr.Dataset()
    coords = ('time', 'z', 'lat', 'lon')
    ds['u'] = (coords, np.float32(data_dict['u']))
    ds['v'] = (coords, np.float32(data_dict['v']))
    ds['u_err'] = (coords, np.float32(data_dict['u_err']))
    ds['v_err'] = (coords, np.float32(data_dict['v_err']))
    ds['uv_covariance'] = (coords, np.float32(data_dict['uv_covariance']))
    ds['num_radials'] = (coords, data_dict['num_radials'])
    ds['qc16_maxspeed'] = (coords, np.int32(data_dict['qc16_maxspeed']))
    ds['qc18_validlocation'] = (coords, np.int32(data_dict['qc18_validlocation']))
    if method == 'oi':
        ds['qc19_uerr'] = (coords, np.int32(data_dict['qc19_uerr']))
        ds['qc20_verr'] = (coords, np.int32(data_dict['qc20_verr']))
    elif method == 'lsq':
        ds['qc15_gdop'] = (coords, np.int32(data_dict['qc15_gdop']))
    ds['qc_primary_flag'] = (coords, np.int32(data_dict['qc_primary_flag']))
    ds['qc_operator_mask'] = (coords, np.int32(data_dict['qc_operator_mask']))

    ds.coords['lon'] = lon
    ds.coords['lat'] = lat
    ds.coords['z'] = np.array([np.float32(0)])
    ds.coords['time'] = time_index

    if flags:
        for k, v in flags.items():
            ds = ds.where(ds[k] <= v)

    ds['processing_parameters'] = (('parameters'), processing_parameters)

    # Grab min and max time in dataset for entry into global attributes for cf compliance
    time_start = ds['time'].min().data
    time_end = ds['time'].max().data

    global_attributes = configs.netcdf_global_attributes(user_attributes, time_start, time_end)

    global_attributes['geospatial_lat_min'] = lat.min()
    global_attributes['geospatial_lat_max'] = lat.max()
    global_attributes['geospatial_lon_min'] = lon.min()
    global_attributes['geospatial_lon_max'] = lon.max()
    if method == 'oi':
        global_attributes['method'] = 'Optimal Interpolation'
    elif method == 'lsq':
        global_attributes['method'] = 'Unweighted Least Squares'

    logging.debug('{} - Assigning global attributes to dataset'.format(fname))
    ds = ds.assign_attrs(global_attributes)

    logging.debug('{} - Assigning local attributes to each variable in dataset'.format(fname))
    # set time attribute
    ds['time'].attrs['standard_name'] = 'time'

    # Set lon attributes
    ds['lon'].attrs['long_name'] = 'Longitude'
    ds['lon'].attrs['standard_name'] = 'longitude'
    ds['lon'].attrs['short_name'] = 'lon'
    ds['lon'].attrs['units'] = 'degrees_east'
    ds['lon'].attrs['axis'] = 'X'
    ds['lon'].attrs['valid_min'] = np.float32(-180.0)
    ds['lon'].attrs['valid_max'] = np.float32(180.0)

    # Set lat attributes
    ds['lat'].attrs['long_name'] = 'Latitude'
    ds['lat'].attrs['standard_name'] = 'latitude'
    ds['lat'].attrs['short_name'] = 'lat'
    ds['lat'].attrs['units'] = 'degrees_north'
    ds['lat'].attrs['axis'] = 'Y'
    ds['lat'].attrs['valid_min'] = np.float32(-90.0)
    ds['lat'].attrs['valid_max'] = np.float32(90.0)

    # Set depth attributes
    ds['z'].attrs['long_name'] = 'Average Depth of Sensor'
    ds['z'].attrs['standard_name'] = 'depth'
    ds['z'].attrs['comment'] = 'Derived from mean value of depth variable'
    ds['z'].attrs['units'] = 'm'
    ds['z'].attrs['axis'] = 'Z'
    ds['z'].attrs['positive'] = 'down'

    # Set u attributes
    ds['u'].attrs['long_name'] = 'Eastward Surface Current (cm/s)'
    ds['u'].attrs['standard_name'] = 'surface_eastward_sea_water_velocity'
    ds['u'].attrs['short_name'] = 'u'
    ds['u'].attrs['units'] = u_units
    ds['u'].attrs['valid_min'] = np.float32(-300)
    ds['u'].attrs['valid_max'] = np.float32(300)
    ds['u'].attrs['coordinates'] = 'lon lat'
    ds['u'].attrs['grid_mapping'] = 'crs'
    if method == 'oi':
        ds['u'].attrs['ancillary_variables'] = 'qc16_maxspeed qc18_validlocation qc19_uerr qc20_verr qc_primary_flag qc_operator_mask'
    elif method == 'lsq':
        ds['u'].attrs['ancillary_variables'] = 'qc15_gdop qc16_maxspeed qc18_validlocation qc_primary_flag qc_operator_mask'

    # Set v attributes
    ds['v'].attrs['long_name'] = 'Northward Surface Current (cm/s)'
    ds['v'].attrs['standard_name'] = 'surface_northward_sea_water_velocity'
    ds['v'].attrs['short_name'] = 'v'
    ds['v'].attrs['units'] = v_units
    ds['v'].attrs['valid_min'] = np.float32(-300)
    ds['v'].attrs['valid_max'] = np.float32(300)
    ds['v'].attrs['coordinates'] = 'lon lat'
    ds['v'].attrs['grid_mapping'] = 'crs'
    if method == 'oi':
        ds['v'].attrs['ancillary_variables'] = 'qc16_maxspeed qc18_validlocation qc19_uerr qc20_verr qc_primary_flag qc_operator_mask'
    elif method == 'lsq':
        ds['v'].attrs['ancillary_variables'] = 'qc15_gdop qc16_maxspeed qc18_validlocation qc_primary_flag qc_operator_mask'


    # Set u_err attributes
    ds['u_err'].attrs['units'] = '1'
    ds['u_err'].attrs['valid_min'] = np.float32(0)
    ds['u_err'].attrs['valid_max'] = np.float32(1)
    ds['u_err'].attrs['coordinates'] = 'lon lat'
    ds['u_err'].attrs['grid_mapping'] = 'crs'

    # Set v_err attributes
    ds['v_err'].attrs['units'] = '1'
    ds['v_err'].attrs['valid_min'] = np.float32(0)
    ds['v_err'].attrs['valid_max'] = np.float32(1)
    ds['v_err'].attrs['coordinates'] = 'lon lat'
    ds['v_err'].attrs['grid_mapping'] = 'crs'

    if method == 'lsq':
        ds['u_err'].attrs['long_name'] = 'Associated GDOP mapping error value associated with eastward velocity component'
        ds['v_err'].attrs['long_name'] = 'Associated GDOP mapping error value associated with northward velocity component'
        ds['total_errors'].attrs['long_name'] = 'Associated GDOP mapping error value, TotalErrors are the square-root of the norm covariance matrix (or the square root of the maximum eigenvalue of the 2x2 UV covariance matrix)'
        ds['u_err'].attrs['comment'] = 'velocity measurements with error values over 1.25 are of questionable quality'
        ds['v_err'].attrs['comment'] = 'velocity measurements with error values over 1.25 are of questionable quality'
        ds['total_errors'].attrs['comment'] = 'velocity measurements with error values over 1.25 are of questionable quality'

    elif method == 'oi':
        ds['u_err'].attrs['long_name'] = 'Normalized uncertainty error associated with eastward velocity component'
        ds['v_err'].attrs['long_name'] = 'Normalized uncertainty error associated with northward velocity component'
        ds['total_errors'].attrs['long_name'] = 'Associated mapping error value, TotalErrors are the square-root of the sum of the squares of U_err and Verr'
        ds['u_err'].attrs['comment'] = 'velocity measurements with error values over 0.6 are of questionable quality'
        ds['v_err'].attrs['comment'] = 'velocity measurements with error values over 0.6 are of questionable quality'
        ds['total_errors'].attrs['comment'] = 'velocity measurements with error values over _ are of questionable quality'


    # Set uv_covariance attributes
    ds['uv_covariance'].attrs['long_name'] = 'Eastward and Northward covariance directional information of u and v'
    ds['uv_covariance'].attrs['units'] = '1'
    ds['uv_covariance'].attrs['comment'] = 'directional information of u and v'
    ds['uv_covariance'].attrs['coordinates'] = 'lon lat'
    ds['uv_covariance'].attrs['grid_mapping'] = 'crs'

    # Set num_radials attributes
    ds['num_radials'].attrs['long_name'] = 'Number of radial measurements used to calculate each totals velocity'
    ds['num_radials'].attrs['comment'] = 'totals are not calculated with fewer than 3 contributing radial measurements from 2 sites'
    ds['num_radials'].attrs['coordinates'] = 'lon lat'
    ds['num_radials'].attrs['grid_mapping'] = 'crs'

    # Set information attributes
    ds['processing_parameters'].attrs['long_name'] = 'General and method specific processing parameter information'
    ds['processing_parameters'].attrs['comment'] = processing_parameters_info
    # ds['processing_parameters'].attrs['coordinates'] = 'parameters'

    # Set qc flag attributes
    ds['qc16_maxspeed'].attrs['standard_name'] = 'sea_water_velocity quality_flag'
    ds['qc16_maxspeed'].attrs['units'] = '1'
    ds['qc16_maxspeed'].attrs['valid_min'] = 1
    ds['qc16_maxspeed'].attrs['valid_max'] = 9
    ds['qc16_maxspeed'].attrs['coordinates'] = 'lon lat'
    ds['qc16_maxspeed'].attrs['grid_mapping'] = 'crs'
    ds['qc16_maxspeed'].attrs['flag_values'] = '1 2 3 4 9'
    ds['qc16_maxspeed'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
    ds['qc16_maxspeed'].attrs['comment'] = qc16_info

    ds['qc18_validlocation'].attrs['standard_name'] = 'sea_water_velocity location_test_quality_flag'
    ds['qc18_validlocation'].attrs['units'] = '1'
    ds['qc18_validlocation'].attrs['valid_min'] = 1
    ds['qc18_validlocation'].attrs['valid_max'] = 9
    ds['qc18_validlocation'].attrs['coordinates'] = 'lon lat'
    ds['qc18_validlocation'].attrs['grid_mapping'] = 'crs'
    ds['qc18_validlocation'].attrs['flag_values'] = '1 2 3 4 9'
    ds['qc18_validlocation'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
    ds['qc18_validlocation'].attrs['comment'] = qc18_info

    if method == 'oi':
        ds['qc19_uerr'].attrs['standard_name'] = 'sea_water_velocity quality_flag'
        ds['qc19_uerr'].attrs['units'] = '1'
        ds['qc19_uerr'].attrs['valid_min'] = 1
        ds['qc19_uerr'].attrs['valid_max'] = 9
        ds['qc19_uerr'].attrs['coordinates'] = 'lon lat'
        ds['qc19_uerr'].attrs['grid_mapping'] = 'crs'
        ds['qc19_uerr'].attrs['flag_values'] = '1 2 3 4 9'
        ds['qc19_uerr'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
        ds['qc19_uerr'].attrs['comment'] = qc19_info

        ds['qc20_verr'].attrs['standard_name'] = 'sea_water_velocity quality_flag'
        ds['qc20_verr'].attrs['units'] = '1'
        ds['qc20_verr'].attrs['valid_min'] = 1
        ds['qc20_verr'].attrs['valid_max'] = 9
        ds['qc20_verr'].attrs['coordinates'] = 'lon lat'
        ds['qc20_verr'].attrs['grid_mapping'] = 'crs'
        ds['qc20_verr'].attrs['flag_values'] = '1 2 3 4 9'
        ds['qc20_verr'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
        ds['qc20_verr'].attrs['comment'] = qc20_info

    elif method == 'lsq':
        ds['qc15_gdop'].attrs['standard_name'] = 'sea_water_velocity quality_flag'
        ds['qc15_gdop'].attrs['units'] = '1'
        ds['qc15_gdop'].attrs['valid_min'] = 1
        ds['qc15_gdop'].attrs['valid_max'] = 9
        ds['qc15_gdop'].attrs['coordinates'] = 'lon lat'
        ds['qc15_gdop'].attrs['grid_mapping'] = 'crs'
        ds['qc15_gdop'].attrs['flag_values'] = '1 2 3 4 9'
        ds['qc15_gdop'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
        ds['qc15_gdop'].attrs['comment'] = qc15_info

    ds['qc_primary_flag'].attrs['standard_name'] = 'sea_water_velocity aggregate_quality_flag'
    ds['qc_primary_flag'].attrs['units'] = '1'
    ds['qc_primary_flag'].attrs['valid_min'] = 1
    ds['qc_primary_flag'].attrs['valid_max'] = 4
    ds['qc_primary_flag'].attrs['coordinates'] = 'lon lat'
    ds['qc_primary_flag'].attrs['grid_mapping'] = 'crs'
    ds['qc_primary_flag'].attrs['flag_values'] = '1 2 3 4'
    ds['qc_primary_flag'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail'
    ds['qc_primary_flag'].attrs['comment'] = qc_primary_flag_info
    if method == 'oi':
        ds['qc_primary_flag'].attrs['ancillary_variables'] = 'qc16_maxspeed qc18_validlocation qc19_uerr qc20_verr'
    elif method == 'lsq':
        ds['qc_primary_flag'].attrs['ancillary_variables'] = 'qc15_gdop qc16_maxspeed qc18_validlocation'

    ds['qc_operator_mask'].attrs['units'] = '1'
    ds['qc_operator_mask'].attrs['valid_min'] = 1
    ds['qc_operator_mask'].attrs['valid_max'] = 9
    ds['qc_operator_mask'].attrs['coordinates'] = 'lon lat'
    ds['qc_operator_mask'].attrs['grid_mapping'] = 'crs'
    ds['qc_operator_mask'].attrs['flag_values'] = '1 2 3 4 9'
    ds['qc_operator_mask'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
    ds['qc_operator_mask'].attrs['comment'] = qc_operator_mask_info

    logging.debug('{} - Setting variable encoding and fill values for netCDF4 output'.format(fname))

    # encode variables for export to netcdf
    encoding = make_encoding(ds)
    encoding['lon'] = dict(zlib=False, _FillValue=False)
    encoding['lat'] = dict(zlib=False, _FillValue=False)
    encoding['z'] = dict(zlib=False, _FillValue=False)

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

    # Create save directory if it doesn't exist.
    create_dir(save_dir)

    logging.debug('{} - Saving dataset to netCDF4 file: {}'.format(fname, file_and_path))
    ds.to_netcdf(file_and_path, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims=['time'])
    logging.info('{} - netCDF4 file successfully created: {}'.format(fname, file_and_path))


if __name__ == '__main__':
    import glob

    # Define test inputs
    files = sorted(glob.glob('../../data/totals/oi/mat/*.mat'))
    grid_file = '../../data/grid_files/maracoos_grid_6km_extended.txt'
    save_dir = '../../data/totals/oi/nc/hourly'
    threshold = dict(u_err=0.6, v_err=0.6, uv_covariance=0.6)

    # load csv file containing the grid
    grid = pd.read_csv(grid_file, sep=',', header=None, names=['lon', 'lat'], delim_whitespace=True)

    user_attributes = dict(title='MARACOOS 6km Sea Surface Currents',
                           naming_authority='edu.rutgers.marine.rucool',
                           comment='Network maintained by MARACOOS. For oi_* global attribute explanations, see references attribute',
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
                           contributor_name='Scott Glenn, Josh Kohut, Hugh Roarty, Ethan Handel, Michael Smith, Laura Nazzaro, Teresa Updyke, Larry Atkinson, Rich Arena, Wendell Brown, Patterson Taylor, Mike Muglia, Harvey Seim',
                           contributor_role='Principal Investigator, Principal Investigator, Principal Investigator, Hardware Maintenance, Data Manager, Data Manager, Principal Investigator, Principal Investigator, Hardware Maintenance, Principal Investigator, Hardware Maintenance, Principal Investigator, Principal Investigator',
                           platform='MARACOOS HF Radar 5MHz Network',
                           instrument='Network includes CODAR sites AMAG, ASSA, BLCK, BRIG, CEDR, CORE, DUCK, FARO, HATY, HEMP, HOOK, LISL, LOVE, MABO, MRCH, MVCO, NANT, NAUS, PYFC, and WILD',
                           references='http://maracoos.org/node/146 https://rucool.marine.rutgers.edu/facilities https://rucool.marine.rutgers.edu/data',
                           summary='Optimally Interpolated Total Vectors calculated by HFRProgs toolbox using MATLAB. Mercator lat/lon projection',
                           ncei_template_version='NCEI_NetCDF_Grid_Template_v2.0',
                           history='Hourly codar radial data combined into one hourly file containing vectors.',
                           cdm_data_type='Grid',
                           source='CODAR SeaSonde Surface Current Mapping Device',
                           processing_level='Level 3',
                           keywords='Environmental Advisories > Marine Advisories > Marine Weather/Forecast, Oceans > Coastal Processes, Oceans > Ocean Circulation, Oceans > Ocean Waves, Oceans > Ocean Winds, Oceans > Ocean Tides, Spectral/Engineering > Radar',
                           publisher_name='NOAA National Centers for Environmental Information',
                           publisher_email='ncei.info@noaa.gov',
                           publisher_url='www.ncei.noaa.gov')

    for file in files:
        main(grid, file, save_dir, user_attributes, threshold)
