#!/usr/bin/env python
"""
@file codar_means.py
@author Mike Smith
@email michaesm@marine.rutgers.edu
@brief This function averages CODAR surface current data over a given time period.
@purpose This function is used to create a mean over a given time period of CODAR surface currents. Mean will only be
calculated for a grid point if that grid point has data present more than 50% of the time. U and V errors (GDOP or GDOSA)
can also be used to filter the dataset.
"""
import os
import xarray as xr
import xarray.ufuncs as xu
from codar_processing.src.common import create_dir

coords = ('z', 'lat', 'lon')


def custom_mean(group, min_coverage=0.5):
    """
    Custom mean function.
    :param group: dataset group created by using the groupby builtin function of an xarray dataset.
    :param min_coverage: Percent of time data must be present for a grid point to be averaged (fraction)
    :return: mean of dataset group.
    """
    # percent coverage. count where grid points are not null and get mean
    percent_coverage = group['u'].notnull().mean('time')

    # filter by the min_coverage input variable
    filtered = group.where(percent_coverage > min_coverage)

    # get standard deviations of u, v and uv variables
    u_stdev = filtered['u'].std('time')
    v_stdev = filtered['v'].std('time')
    uv_stdev = xu.sqrt(u_stdev ** 2 + v_stdev ** 2)

    # get the mean of u and v
    result = filtered[['u', 'v']].mean('time')

    # calculate magnitude
    mag = xu.sqrt(result['u'] ** 2 + result['v'] ** 2)

    # add variables  to this group
    result['percent_coverage'] = (coords, percent_coverage)
    result['magnitude'] = (coords, mag)
    result['u_stdev'] = (coords, u_stdev)
    result['v_stdev'] = (coords, v_stdev)
    result['uv_stdev'] = (coords, uv_stdev)

    return result


def filter_by_error(ds, g_u=1, g_v=1):
    """
    Filter out data where both U and V contain errors above a threshold
    :param ds: Xarray dataset containing codar surface current data
    :param g_u: U error threshold
    :param g_v: V error threshold
    :return:
    """
    ds = ds.where((ds.u_err < g_u) & (ds.v_err < g_v))
    return ds


def main(files, save_dir, average, error):
    """
    This function averages CODAR surface current data over a given time period.
    :param fname: Path to netCDF file or directory containing netCDF files. If directory, must use wildcard (/path/*.nc)
    :param save_dir: Directory to save averaged surface current netCDF file
    :param average: Average by time. Default: time.month; See: http://xarray.pydata.org/en/latest/time-series.html#datetime-components
    :param g_u: Error in Eastward seawater velocity (u). Default: .6
    :param g_v: Error in Northward seawater velocity (v). Default: .6
    """
    # open netcdf file or files
    ds = xr.open_mfdataset(files)

    # filter both u and v by .6. anything over will be removed
    ds = filter_by_error(ds, error[0], error[1])

    # group dataset by
    new_ds = ds.groupby(average).apply(custom_mean)

    # Create save directory if it doesn't exist
    create_dir(save_dir)

    # Output dataset containing grouped averages to netcdf
    new_ds.to_netcdf(os.path.join(save_dir, 'codar_monthly_average.nc'))


if __name__ == '__main__':
    files = '../data/totals/nc/hourly/*.nc'
    save_dir = '../data/totals/nc/averaged'
    average = 'time.month'
    filter_error_u_v = [.6, .6]
    main(files, save_dir, average, filter_error_u_v)