#!/usr/bin/env python
import xarray as xr

buoy_ncml = 'http://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/44009/44009.ncml'
start_time = '2006-12-01'
end_time = '2007-01-01'
# end_time = '2016-12-31 23:00:00'

met_data = xr.open_dataset(buoy_ncml)
subset = met_data.sel(time=slice(start_time, end_time))
wind_seasons = dict(subset.groupby('time.season'))

wind_angles = dict(DJF=[285, 345], MAM=[165, 225], JJA=[165, 225], SON=[15, 75])

for key in wind_seasons.keys():
    tds = wind_seasons[key]
    tds = tds['wind_dir'].where((wind_angles[key][0] < tds['wind_dir']) & (tds['wind_dir'] < wind_angles[key][-1]))
    print


# print
# times = ['2010-01-01 00:00:00',
#          '2010-02-01 00:00:00',
#          '2010-03-01 00:00:00',
#          '2010-04-01 00:00:00',
#          '2010-05-01 00:00:00',
#          '2010-06-01 00:00:00',
#          '2010-07-01 00:00:00',
#          '2010-08-01 00:00:00',
#          '2010-09-01 00:00:00',
#          '2010-10-01 00:00:00',
#          '2010-11-01 00:00:00',
#          '2010-12-01 00:00:00',]
#
# # def filter_wind(season, data):
#
# ds = xr.open_mfdataset('/Volumes/boardwalk/codar/codar_means/xarray_aggregated/2010*.nc')
#
# for t in times:
#     arr = xr.concat([arr, ds.sel(time=t)], dim='time')