#!/usr/bin/env python
import os
import xarray as xr
import xarray.ufuncs as xu
from codar_processing.common import create_dir

coords = ('z', 'lat', 'lon')


def custom_mean(group, min_coverage=0.5):
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


def filter_by_gdosa(ds, g_u=1, g_v=1):
    """
    Filter out data where both U and V contain errors above a threshold
    :param ds: netcdf file in the form of an xarray dataset
    :param g_u: gdosa u threshold
    :param g_v: dgosa v threshold
    :return:
    """
    ds = ds.where((ds.u_err < g_u) & (ds.v_err < g_v))
    return ds


def main(fname, save_dir, average_by, g_u=.6, g_v=.6):
    # open netcdf file or files
    ds = xr.open_mfdataset(fname)

    # filter both u and v by .6. anything over will be removed
    ds = filter_by_gdosa(ds, g_u, g_v)

    # group dataset by
    new_ds = ds.groupby(average_by).apply(custom_mean)

    # Create save directory if it doesn't exist
    create_dir(save_dir)

    # Output dataset containing grouped averages to netcdf
    new_ds.to_netcdf(os.path.join(save_dir, 'codar_monthly_average.nc'))


if __name__ == '__main__':
    files = '../data/totals/nc/*.nc'
    save_dir = '../data/totals/nc_averaged'
    average_by = 'time.month'
    main(files, save_dir, average_by)