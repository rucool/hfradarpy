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
import re

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
    ## time = dt.datetime.strptime(mat_time, '%Y_%m_%d_%H00')

    #timestamp_regex = re.compile('\d{12}')
    #mat_time = timestamp_regex.search(fname).group()
    #time = dt.datetime.strptime(mat_time, '%Y%m%d%H%M')


    time_index = pd.date_range(time.strftime('%Y-%m-%d %H:%M:%S'), periods=1)  # create pandas datetimeindex from time
    ##time_string = time.strftime('%Y%m%dT%H%M%SZ')  # create timestring from time
    #time_string = time.strftime('%Y%m%d%H%M')  # create timestring from time

    #file_name = 'RU_{}_{}.nc'.format(domain, time_string)
    #file_name = '{}_hfr_midatl_6km_rtv_oi_maracoos.nc'.format(time_string,domain)
    file_name = fname[:-3] + 'nc'
    file_and_path = os.path.join(save_dir, file_name)

    try:
        logging.debug('{} - Saving file data to variables'.format(fname))
        # load longitude and latitude data associated with variables
        lonlat = data['TUV'].LonLat.astype(np.float32)

        # create variables for eastward and northward velocities
        u_vel = data['TUV'].U.astype(np.float32)
        v_vel = data['TUV'].V.astype(np.float32)
        # rounds to the whole number for velocity
        # similar to HFRNet we will not report fractions of a cm/s and the scale factors of 0.01 assigned for velocity variables
        # in the code below will indicate to data users that the units are m/s when the values in the array
        # are multiplied by the scale factor
        # HFRNet did this in part so they would be able to assign a variable type of short
        u = np.round(u_vel).astype(np.float32)
        v = np.round(v_vel).astype(np.float32)
        u_units = 'm/s'
        v_units = 'm/s'





        radial_metadata = ','.join(data['TUVmetadata'].radial_metadata)
        #radmeta = '\n,'.join(data['TUVmetadata'].radial_metadata)
        #quote = '"'
        #radial_metadata = quote + radmeta + quote

        radial_num_sites = data['TUVmetadata'].radial_num_sites

        #maxspd = data['TUV'].OtherMetadata.cleanTotals.maxspd
        maxspd = data['TUVmetadata'].conf.Totals.MaxTotSpeed

        flag_values = np.array([1,2,3,4,9],dtype=np.byte)

        if method == 'oi':
            # create variables for associated error values
            u_err = data['TUV'].ErrorEstimates.Uerr.astype(np.float32)
            v_err = data['TUV'].ErrorEstimates.Verr.astype(np.float32)
            uv_covariance = data['TUV'].ErrorEstimates.UVCovariance.astype(np.float32)
            total_errors = data['TUV'].ErrorEstimates.TotalErrors.astype(np.float32)

            # Data Processing Information
            num_rads = data['TUV'].OtherMatrixVars.makeTotalsOI_TotalsNumRads.astype(np.short)
            min_rads = data['TUV'].OtherMetadata.makeTotalsOI.parameters.MinNumRads
            min_rads = np.byte(min_rads)
            min_sites = data['TUV'].OtherMetadata.makeTotalsOI.parameters.MinNumSites
            min_sites = np.byte(min_sites)
            mdlvar = data['TUV'].OtherMetadata.makeTotalsOI.parameters.mdlvar
            mdlvar = np.short(mdlvar)
            errvar = data['TUV'].OtherMetadata.makeTotalsOI.parameters.errvar
            errvar = np.short(errvar)
            sx = data['TUV'].OtherMetadata.makeTotalsOI.parameters.sx
            sx = np.float32(sx)
            sy = data['TUV'].OtherMetadata.makeTotalsOI.parameters.sy
            sy = np.float32(sy)
            temporal_threshold = data['TUV'].OtherMetadata.makeTotalsOI.parameters.tempthresh
            temporal_threshold = np.float32(temporal_threshold)

            #QC Information
            uerr_testname = data['TUVmetadata'].conf.Totals.OI.cleanTotalsVarargin[0][0] + ' ' + data['TUVmetadata'].conf.Totals.OI.cleanTotalsVarargin[0][1]
            uerr_threshold = data['TUVmetadata'].conf.Totals.OI.cleanTotalsVarargin[0][2]
            uerr_threshold = np.float32(uerr_threshold)
            verr_testname = data['TUVmetadata'].conf.Totals.OI.cleanTotalsVarargin[1][0] + ' ' + data['TUVmetadata'].conf.Totals.OI.cleanTotalsVarargin[1][1]
            verr_threshold = data['TUVmetadata'].conf.Totals.OI.cleanTotalsVarargin[1][2]
            verr_threshold = np.float32(verr_threshold)
            qc_info = 'Quality control reference: IOOS QARTOD HF Radar ver 1.0 May 2016\n'
            qc_info += 'QCFlagDefinitions: 1 = pass, 2 = not_evaluated, 3 = suspect, 4 = fail, 9 = missing_data\n'
            qc_primary_flag_info = 'QCPrimaryFlagDefinition: Highest flag value of qc303, qc305, qc306, qc307\n'
            qc_primary_flag_info += 'This flag will be set to not_evaluated only if ALL individual tests were not_evaluated.'
            qc_operator_flag_info = qc_info + 'The qc_operator_flag follows QCFlagDefinitions and is set at discretion of the operator or data manager.'
            qc303_info = qc_info + 'qc303 Max Speed Threshold [max_vel = ' + str(maxspd) + ' (cm/s)]'
            qc305_info = qc_info + 'qc305 Valid Location [landmask file = ' + data['TUVmetadata'].conf.Totals.MaskFile + ']'
            qc306_info = qc_info + 'qc306 OI Uncertainty Threshold [' + uerr_testname + ' ' + str(uerr_threshold) + ']'
            qc307_info = qc_info + 'qc307 OI Uncertainty Threshold [' + verr_testname + ' ' + str(verr_threshold) + ']'
            qc303 = data['TUV'].QC303
            qc305 = data['TUV'].QC305
            qc306 = data['TUV'].QC306
            qc307 = data['TUV'].QC307
            qc_primary_flag = data['TUV'].PRIM
            qc_operator_flag = data['TUV'].qc_operator_flag
            qc303 = np.byte(qc303)
            qc305 = np.byte(qc305)
            qc306 = np.byte(qc306)
            qc307 = np.byte(qc307)
            qc_primary_flag = np.byte(qc_primary_flag)
            qc_operator_flag = np.byte(qc_operator_flag)


        elif method == 'lsq':

            # create variables for associated error values
            uerr_data = data['TUV'].ErrorEstimates[1].Uerr
            verr_data = data['TUV'].ErrorEstimates[1].Verr
            uvcov_data = data['TUV'].ErrorEstimates[1].UVCovariance
            totalerrors_data = data['TUV'].ErrorEstimates[1].TotalErrors
            # rounds errors to the whole number (cm^2/s^2 for u_err, v_err, uv_covariance and cm/s for total_errors)
            # when scale factor is applied the units will be m^2/s^2 and m/s
            u_err = np.round(uerr_data).astype(np.float32)
            v_err = np.round(verr_data).astype(np.float32)
            uv_covariance = np.round(uvcov_data).astype(np.float32)
            total_errors = np.round(totalerrors_data).astype(np.float32)

            # create variables for associated error values
            #u_err = data['TUV'].ErrorEstimates[1].Uerr.astype(np.float32)
            #v_err = data['TUV'].ErrorEstimates[1].Verr.astype(np.float32)
            #uv_covariance = data['TUV'].ErrorEstimates[1].UVCovariance.astype(np.float32)
            #total_errors = data['TUV'].ErrorEstimates[1].TotalErrors.astype(np.float32)

            # create variables for eastward and northward velocities
            u_vel = data['TUV'].U.astype(np.float32)
            v_vel = data['TUV'].V.astype(np.float32)
            # rounds to the whole number for velocity
            u = np.round(u_vel).astype(np.float32)
            v = np.round(v_vel).astype(np.float32)
            u_units = 'm/s'
            v_units = 'm/s'

            # Data Processing Information
            num_rads = data['TUV'].OtherMatrixVars.makeTotals_TotalsNumRads.astype(np.short)
            min_rads = data['TUV'].OtherMetadata.makeTotals.parameters.MinNumRads
            min_rads = np.byte(min_rads)
            min_sites = data['TUV'].OtherMetadata.makeTotals.parameters.MinNumSites
            min_sites = np.byte(min_sites)
            spatial_threshold = data['TUV'].OtherMetadata.makeTotals.parameters.spatthresh
            spatial_threshold = np.float32(spatial_threshold)
            temporal_threshold = data['TUV'].OtherMetadata.makeTotals.parameters.tempthresh
            temporal_threshold = np.float32(temporal_threshold)
            #processing_parameters = [maxspd, min_sites, min_rads, temporal_threshold, spatial_threshold]
            #processing_parameters = [min_sites, min_rads, temporal_threshold, spatial_threshold]
            #processing_parameters_info = '1) Maximum Total Speed Threshold (cm s-1)\n'
            #processing_parameters_info = '1) Minimum number of radial sites.\n'
            #processing_parameters_info += '2) Minimum number of radial vectors.\n'
            #processing_parameters_info += '3) Temporal search window for radial solutions (Fractions of a day)\n'
            #processing_parameters_info += '4) Spatial search radius for radial solutions (km)\n'


            #QC Information
            gdoptestname = data['TUVmetadata'].conf.Totals.cleanTotalsVarargin[0] + ' ' + data['TUVmetadata'].conf.Totals.cleanTotalsVarargin[1]
            gdopthreshold = data['TUVmetadata'].conf.Totals.cleanTotalsVarargin[2]
            gdopthreshold = np.float32(gdopthreshold)
            qc_info = 'Quality control reference: IOOS QARTOD HF Radar ver 1.0 May 2016\n'
            qc_info += 'QCFlagDefinitions: 1 = pass, 2 = not_evaluated, 3 = suspect, 4 = fail, 9 = missing_data\n'
            qc_primary_flag_info = 'QCPrimaryFlagDefinition: Highest flag value of qc303, qc305, qc306, qc307\n'
            qc_primary_flag_info += 'This flag will be set to not_evaluated only if ALL tests were not_evaluated.'
            qc_operator_flag_info = qc_info + 'The qc_operator_flag follows QCFlagDefinitions and is set at discretion of the operator or data manager.'
            qc303_info = qc_info + 'qc303 Max Speed Threshold [max_vel = ' + str(maxspd) + ' (cm/s)]'
            qc305_info = qc_info + 'qc305 Valid Location [landmask file = ' + data['TUVmetadata'].conf.Totals.MaskFile + ']'
            qc302_info = qc_info + 'qc302 GDOP Threshold [' + gdoptestname + ' ' + str(gdopthreshold) + ']'
            qc302 = data['TUV'].QC302
            qc303 = data['TUV'].QC303
            qc305 = data['TUV'].QC305
            qc_primary_flag = data['TUV'].PRIM
            qc_operator_flag = data['TUV'].qc_operator_flag
            qc302 = np.byte(qc302)
            qc303 = np.byte(qc303)
            qc305 = np.byte(qc305)
            qc_primary_flag = np.byte(qc_primary_flag)
            qc_operator_flag = np.byte(qc_operator_flag)

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
                         number_of_radials=num_rads,
                         qc303_maxspeed=qc303,
                         qc305_validlocation=qc305,
                         qc306_uerr=qc306,
                         qc307_verr=qc307,
                         qc_primary_flag=qc_primary_flag,
                         qc_operator_flag=qc_operator_flag,
                         )
    elif method == 'lsq':
        data_dict = dict(u=u,
                         v=v,
                         u_err=u_err,
                         v_err=v_err,
                         uv_covariance=uv_covariance,
                         total_errors = total_errors,
                         number_of_radials=num_rads,
                         qc302_total_errors=qc302,
                         qc303_maxspeed=qc303,
                         qc305_validlocation=qc305,
                         qc_primary_flag=qc_primary_flag,
                         qc_operator_flag=qc_operator_flag,
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
    ds['total_errors'] = (coords, np.float32(data_dict['total_errors']))
    ds['uv_covariance'] = (coords, np.float32(data_dict['uv_covariance']))
    ds['number_of_radials'] = (coords, np.short(data_dict['number_of_radials']))
    ds['qc303_maxspeed'] = (coords, np.byte(data_dict['qc303_maxspeed']))
    ds['qc305_validlocation'] = (coords, np.byte(data_dict['qc305_validlocation']))
    if method == 'oi':
        ds['qc306_uerr'] = (coords, np.byte(data_dict['qc306_uerr']))
        ds['qc307_verr'] = (coords, np.byte(data_dict['qc307_verr']))
    elif method == 'lsq':
        ds['qc302_gdop'] = (coords, np.byte(data_dict['qc302_gdop']))
    ds['qc_primary_flag'] = (coords, np.byte(data_dict['qc_primary_flag']))
    ds['qc_operator_flag'] = (coords, np.byte(data_dict['qc_operator_flag']))

    ds.coords['lon'] = lon
    ds.coords['lat'] = lat
    ds.coords['z'] = np.array([np.float32(0)])
    ds.coords['time'] = time_index

    if flags:
        for k, v in flags.items():
            ds = ds.where(ds[k] <= v)

    #ds['processing_parameters'] = (('parameters'), processing_parameters)

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
    global_attributes['MARACOOS_processing_version'] = str(data['TUV'].MARACOOS_processing_version)
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
    ds['u'].attrs['long_name'] = 'Eastward Surface Current (m/s)'
    ds['u'].attrs['standard_name'] = 'surface_eastward_sea_water_velocity'
    ds['u'].attrs['short_name'] = 'u'
    ds['u'].attrs['units'] = u_units
    ds['u'].attrs['scale_factor'] = np.float32(0.01)
    ds['u'].attrs['valid_min'] = np.float32(-300)
    ds['u'].attrs['valid_max'] = np.float32(300)
    ds['u'].attrs['coordinates'] = 'lon lat'
    ds['u'].attrs['grid_mapping'] = 'crs'
    if method == 'oi':
        ds['u'].attrs['ancillary_variables'] = 'qc303_maxspeed qc305_validlocation qc306_uerr qc307_verr qc_primary_flag qc_operator_flag'
    elif method == 'lsq':
        ds['u'].attrs['ancillary_variables'] = 'qc302_gdop qc303_maxspeed qc305_validlocation qc_primary_flag qc_operator_flag'

    # Set v attributes
    ds['v'].attrs['long_name'] = 'Northward Surface Current (m/s)'
    ds['v'].attrs['standard_name'] = 'surface_northward_sea_water_velocity'
    ds['v'].attrs['short_name'] = 'v'
    ds['v'].attrs['units'] = v_units
    ds['v'].attrs['scale_factor'] = np.float32(0.01)
    ds['v'].attrs['valid_min'] = np.float32(-300)
    ds['v'].attrs['valid_max'] = np.float32(300)
    ds['v'].attrs['coordinates'] = 'lon lat'
    ds['v'].attrs['grid_mapping'] = 'crs'
    if method == 'oi':
        ds['v'].attrs['ancillary_variables'] = 'qc303_maxspeed qc305_validlocation qc306_uerr qc307_verr qc_primary_flag qc_operator_flag'
    elif method == 'lsq':
        ds['v'].attrs['ancillary_variables'] = 'qc302_gdop qc303_maxspeed qc305_validlocation qc_primary_flag qc_operator_flag'

    # Set u_err attributes
    ds['u_err'].attrs['short_name'] = 'uerr'
    ds['u_err'].attrs['units'] = '1'
    if method == 'lsq':
        ds['u_err'].attrs['units'] = 'm^2/s^2'
        ds['u_err'].attrs['scale_factor'] = np.float32(0.0001)
    elif method == 'oi':
        ds['u_err'].attrs['units'] = '1'
    ds['u_err'].attrs['valid_min'] = np.float32(0)
    ds['u_err'].attrs['valid_max'] = np.float32(1)
    ds['u_err'].attrs['coordinates'] = 'lon lat'
    ds['u_err'].attrs['grid_mapping'] = 'crs'

    # Set v_err attributes
    ds['v_err'].attrs['short_name'] = 'verr'
    if method == 'lsq':
        ds['v_err'].attrs['units'] = 'm^2/s^2'
        ds['v_err'].attrs['scale_factor'] = np.float32(0.0001)
    elif method == 'oi':
        ds['v_err'].attrs['units'] = '1'
    ds['v_err'].attrs['valid_min'] = np.float32(0)
    ds['v_err'].attrs['valid_max'] = np.float32(1)
    ds['v_err'].attrs['coordinates'] = 'lon lat'
    ds['v_err'].attrs['grid_mapping'] = 'crs'

    # Set total_errors attributes
    ds['total_errors'].attrs['short_name'] = 'totalerr'
    if method == 'lsq':
        ds['total_errors'].attrs['units'] = 'm/s'
        ds['total_errors'].attrs['scale_factor'] = np.float32(0.01)
    elif method == 'oi':
        ds['total_errors'].attrs['units'] = '1'
    ds['total_errors'].attrs['valid_min'] = np.float32(-300)
    ds['total_errors'].attrs['valid_max'] = np.float32(300)
    ds['total_errors'].attrs['comment'] = 'directional information of u and v'
    ds['total_errors'].attrs['coordinates'] = 'lon lat'
    ds['total_errors'].attrs['grid_mapping'] = 'crs'

    # Set uv_covariance attributes    ds['uv_covariance'].attrs['long_name'] = 'Eastward and Northward covariance directional information of u and v'
    ds['uv_covariance'].attrs['short_name'] = 'uvcov'
    if method == 'lsq':
        ds['uv_covariance'].attrs['units'] = 'm^2/s^2'
        ds['uv_covariance'].attrs['scale_factor'] = np.float32(0.0001)
    elif method == 'oi':
        ds['uv_covariance'].attrs['units'] = '1'
    ds['uv_covariance'].attrs['valid_min'] = np.float32(-300)
    ds['uv_covariance'].attrs['valid_max'] = np.float32(300)
    ds['uv_covariance'].attrs['comment'] = 'directional information of u and v'
    ds['uv_covariance'].attrs['coordinates'] = 'lon lat'
    ds['uv_covariance'].attrs['grid_mapping'] = 'crs'


    if method == 'lsq':
        ds['u_err'].attrs['long_name'] = 'Associated GDOP mapping error value associated with eastward velocity component'
        ds['v_err'].attrs['long_name'] = 'Associated GDOP mapping error value associated with northward velocity component'
        ds['total_errors'].attrs['long_name'] = 'Associated GDOP mapping error value, TotalErrors are the square-root of the norm covariance matrix (or the square root of the maximum eigenvalue of the 2x2 UV covariance matrix)'
        ds['u_err'].attrs['comment'] = 'velocity measurements with high u_err values are of questionable quality'
        ds['v_err'].attrs['comment'] = 'velocity measurements with high v_err values are of questionable quality'
        ds['total_errors'].attrs['comment'] = 'velocity measurements with total_error values over 1.25 are of questionable quality'

    elif method == 'oi':
        ds['u_err'].attrs['long_name'] = 'Normalized uncertainty associated with eastward velocity component, <(u_hat - u)^2>/<u^2>'
        ds['v_err'].attrs['long_name'] = 'Normalized uncertainty associated with northward velocity component, <(v_hat - v)^2>/<v^2>'
        ds['uv_covariance'].attrs['long_name'] = 'Directional information of u and v, <(u_hat -u)(v_hat- v)>/sqrt(<u^2><v^2>)'
        ds['total_errors'].attrs['long_name'] = 'Associated mapping error value, sqrt(U_err^2 + Verr^2)'
        ds['u_err'].attrs['comment'] = 'normalized by constant signal covariance for all gridpoints, velocity measurements with u_err values over 0.6 are of questionable quality'
        ds['v_err'].attrs['comment'] = 'normalized by constant signal covariance for all gridpoints, velocity measurements with v_err values over 0.6 are of questionable quality'
        ds['total_errors'].attrs['comment'] = 'velocity measurements with high total_error values are of questionable quality'

    # Set number_of_radials attributes
    ds['number_of_radials'].attrs['long_name'] = 'Number of radial measurements used to calculate each totals velocity'
    ds['number_of_radials'].attrs['short_name'] = 'num_rads'
    ds['number_of_radials'].attrs['units'] = '1'
    ds['number_of_radials'].attrs['comment'] = 'totals are not calculated with fewer than 3 contributing radial measurements from 2 sites'
    ds['number_of_radials'].attrs['coordinates'] = 'lon lat'
    ds['number_of_radials'].attrs['grid_mapping'] = 'crs'

    # Set information attributes

    # Set qc flag attributes
    ds['qc303_maxspeed'].attrs['standard_name'] = 'surface_eastward_sea_water_velocity status_flag'
    ds['qc303_maxspeed'].attrs['long_name'] = 'maximum speed test flag'
    ds['qc303_maxspeed'].attrs['short_name'] = 'qc303'
    ds['qc303_maxspeed'].attrs['units'] = '1'
    ds['qc303_maxspeed'].attrs['valid_min'] = np.byte(1)
    ds['qc303_maxspeed'].attrs['valid_max'] = np.byte(9)
    ds['qc303_maxspeed'].attrs['coordinates'] = 'lon lat'
    ds['qc303_maxspeed'].attrs['grid_mapping'] = 'crs'
    ds['qc303_maxspeed'].attrs['flag_values'] = flag_values
    ds['qc303_maxspeed'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
    ds['qc303_maxspeed'].attrs['comment'] = qc303_info

    ds['qc305_validlocation'].attrs['standard_name'] = 'surface_eastward_sea_water_velocity status_flag'
    ds['qc305_validlocation'].attrs['long_name'] = 'valid location test flag'
    ds['qc305_validlocation'].attrs['short_name'] = 'qc305'
    ds['qc305_validlocation'].attrs['units'] = '1'
    ds['qc305_validlocation'].attrs['valid_min'] = np.byte(1)
    ds['qc305_validlocation'].attrs['valid_max'] = np.byte(9)
    ds['qc305_validlocation'].attrs['coordinates'] = 'lon lat'
    ds['qc305_validlocation'].attrs['grid_mapping'] = 'crs'
    ds['qc305_validlocation'].attrs['flag_values'] = flag_values
    ds['qc305_validlocation'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
    ds['qc305_validlocation'].attrs['comment'] = qc305_info

    if method == 'oi':
        ds['qc306_uerr'].attrs['standard_name'] = 'surface_eastward_sea_water_velocity status_flag'
        ds['qc306_uerr'].attrs['long_name'] = 'Uerr flag'
        ds['qc306_uerr'].attrs['short_name'] = 'qc306'
        ds['qc306_uerr'].attrs['units'] = '1'
        ds['qc306_uerr'].attrs['valid_min'] = np.byte(1)
        ds['qc306_uerr'].attrs['valid_max'] = np.byte(9)
        ds['qc306_uerr'].attrs['coordinates'] = 'lon lat'
        ds['qc306_uerr'].attrs['grid_mapping'] = 'crs'
        ds['qc306_uerr'].attrs['flag_values'] = flag_values
        ds['qc306_uerr'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
        ds['qc306_uerr'].attrs['comment'] = qc306_info

        ds['qc307_verr'].attrs['standard_name'] = 'surface_northward_sea_water_velocity status_flag'
        ds['qc307_verr'].attrs['long_name'] = 'Verr flag'
        ds['qc307_verr'].attrs['short_name'] = 'qc307'
        ds['qc307_verr'].attrs['units'] = '1'
        ds['qc307_verr'].attrs['valid_min'] = np.byte(1)
        ds['qc307_verr'].attrs['valid_max'] = np.byte(9)
        ds['qc307_verr'].attrs['coordinates'] = 'lon lat'
        ds['qc307_verr'].attrs['grid_mapping'] = 'crs'
        ds['qc307_verr'].attrs['flag_values'] = flag_values
        ds['qc307_verr'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
        ds['qc307_verr'].attrs['comment'] = qc307_info

    elif method == 'lsq':
        ds['qc302_gdop'].attrs['standard_name'] = 'surface_eastward_sea_water_velocity status_flag'
        ds['qc302_gdop'].attrs['long_name'] = 'GDOP mapping error flag'
        ds['qc302_gdop'].attrs['short_name'] = 'qc302'
        ds['qc302_gdop'].attrs['units'] = '1'
        ds['qc302_gdop'].attrs['valid_min'] = np.byte(1)
        ds['qc302_gdop'].attrs['valid_max'] = np.byte(9)
        ds['qc302_gdop'].attrs['coordinates'] = 'lon lat'
        ds['qc302_gdop'].attrs['grid_mapping'] = 'crs'
        ds['qc302_gdop'].attrs['flag_values'] = flag_values
        ds['qc302_gdop'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
        ds['qc302_gdop'].attrs['comment'] = qc302_info

    ds['qc_primary_flag'].attrs['standard_name'] = 'surface_eastward_sea_water_velocity status_flag'
    ds['qc_primary_flag'].attrs['long_name'] = 'Primary Flag'
    ds['qc_primary_flag'].attrs['short_name'] = 'primaryflag'
    ds['qc_primary_flag'].attrs['units'] = '1'
    ds['qc_primary_flag'].attrs['valid_min'] = np.byte(1)
    ds['qc_primary_flag'].attrs['valid_max'] = np.byte(9)
    ds['qc_primary_flag'].attrs['coordinates'] = 'lon lat'
    ds['qc_primary_flag'].attrs['grid_mapping'] = 'crs'
    ds['qc_primary_flag'].attrs['flag_values'] = flag_values
    ds['qc_primary_flag'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
    ds['qc_primary_flag'].attrs['comment'] = qc_primary_flag_info
    if method == 'oi':
        ds['qc_primary_flag'].attrs['ancillary_variables'] = 'qc303_maxspeed qc305_validlocation qc306_uerr qc307_verr'
    elif method == 'lsq':
        ds['qc_primary_flag'].attrs['ancillary_variables'] = 'qc302_gdop qc303_maxspeed qc305_validlocation'

    ds['qc_operator_flag'].attrs['standard_name'] = 'surface_eastward_sea_water_velocity status_flag'
    ds['qc_operator_flag'].attrs['long_name'] = 'Operator Flag'
    ds['qc_operator_flag'].attrs['short_name'] = 'opflag'
    ds['qc_operator_flag'].attrs['units'] = '1'
    ds['qc_operator_flag'].attrs['valid_min'] = np.byte(1)
    ds['qc_operator_flag'].attrs['valid_max'] = np.byte(9)
    ds['qc_operator_flag'].attrs['coordinates'] = 'lon lat'
    ds['qc_operator_flag'].attrs['grid_mapping'] = 'crs'
    ds['qc_operator_flag'].attrs['flag_values'] = flag_values
    ds['qc_operator_flag'].attrs['flag_meanings'] = 'pass not_evaluated suspect fail missing_data'
    ds['qc_operator_flag'].attrs['comment'] = qc_operator_flag_info

    logging.debug('{} - Setting variable encoding and fill values for netCDF4 output'.format(fname))

    # encode variables for export to netcdf
    encoding = make_encoding(ds)
    #encoding['lon'] = dict(zlib=False, _FillValue=False)
    #encoding['lat'] = dict(zlib=False, _FillValue=False)
    #encoding['z'] = dict(zlib=False, _FillValue=False)
    encoding['lon'] = dict(zlib=False, _FillValue=999)
    encoding['lat'] = dict(zlib=False, _FillValue=999)
    encoding['z'] = dict(zlib=False, _FillValue=999)

    # add container variables that contain no data
    kwargs = dict(crs=np.byte(0), instrument=np.byte(0), radial_metadata=np.byte(0), processing_parameters=np.byte(0))
    ds = ds.assign(**kwargs)

    # Set crs attributes
    ds['crs'].attrs['grid_mapping_name'] = 'latitude_longitude'
    ds['crs'].attrs['inverse_flattening'] = 298.257223563
    ds['crs'].attrs['long_name'] = 'Coordinate Reference System'
    ds['crs'].attrs['semi_major_axis'] = '6378137.0'
    ds['crs'].attrs['epsg_code'] = 'EPSG:4326'
    ds['crs'].attrs['comment'] = 'http://www.opengis.net/def/crs/EPSG/0/4326'

    ds['instrument'].attrs['long_name'] = 'CODAR SeaSonde High Frequency Radar'
    ds['processing_parameters'].attrs['units'] = '1'
    ds['instrument'].attrs['sensor_type'] = 'Direction-finding high frequency radar antenna'
    ds['instrument'].attrs['make_model'] = 'CODAR SeaSonde'
    ds['instrument'].attrs['serial_number'] = 1

    ds['radial_metadata'].attrs['long_name'] = 'Metadata on radial velocities used to compute total solutions'
    ds['radial_metadata'].attrs['short_name'] = 'radmeta'
    ds['processing_parameters'].attrs['units'] = '1'
    ds['radial_metadata'].attrs['number_files_loaded'] = str(radial_num_sites)
    ds['radial_metadata'].attrs['number_files_loaded_description'] = 'Number of radial files loaded'
    ds['radial_metadata'].attrs['files_loaded'] = radial_metadata
    ds['radial_metadata'].attrs['files_loaded_description'] = 'Radial file names loaded'

    ds['processing_parameters'].attrs['long_name'] = 'General and method specific processing parameter information'
    ds['processing_parameters'].attrs['short_name'] = 'params'
    ds['processing_parameters'].attrs['units'] = '1'
    ds['processing_parameters'].attrs['parameter1'] = min_sites
    ds['processing_parameters'].attrs['parameter1_description'] = 'Minimum number of radial sites'
    ds['processing_parameters'].attrs['parameter2'] = min_rads
    ds['processing_parameters'].attrs['parameter2_description'] = 'Minimum number of radial vectors'
    ds['processing_parameters'].attrs['parameter3'] = temporal_threshold
    ds['processing_parameters'].attrs['parameter3_description'] = 'Temporal search window for radial solutions (Fraction of a day)'

    if method == 'oi':
        ds['processing_parameters'].attrs['parameter4'] = sx
        ds['processing_parameters'].attrs['parameter4_description'] = 'Decorrelation scales in the east direction'
        ds['processing_parameters'].attrs['parameter5'] = sy
        ds['processing_parameters'].attrs['parameter5_description'] = 'Decorrelation scales in the north direction'
        ds['processing_parameters'].attrs['parameter6'] = mdlvar
        ds['processing_parameters'].attrs['parameter6_description'] = 'Signal variance of the surface current fields (cm2 s-2)'
        ds['processing_parameters'].attrs['parameter7'] = errvar
        ds['processing_parameters'].attrs['parameter7_description'] = 'Data error variance of the input radial velocities (cm2 s-2)'
    elif method == 'lsq':
        ds['processing_parameters'].attrs['parameter4'] = spatial_threshold
        ds['processing_parameters'].attrs['parameter4_description'] = 'Spatial search radius for radial solutions (km)'

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
                           instrument='Network includes CODAR sites AMAG, ASSA, BLCK, BRIG, CEDR, CORE, DUCK, HATY, HEMP, HOOK, LISL, LOVE, MRCH, MVCO, NANT, NAUS, and WILD',
                           references='DOI:10.1029/2007JC004244 DOI:10.1007/s10236-012-0533-9 http://maracoos.org https://rucool.marine.rutgers.edu/facilities https://rucool.marine.rutgers.edu/data',
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
