import functions.common as cf
import functions.radials as rf
import numpy as np
import os
import re
import pandas as pd
import xarray as xr
from functions.calc import gridded_index, reckon


def main(radial):
    radial_file_data = cf.parse_lluv(radial)

    # with open(radial, 'r') as rdl_file:
    #     header_dict, dist_dict, columns, units = cf.parse_header(rdl_file)
    #     header_dict = cf.clean_radial_header(header_dict)

    df = radial_file_data['tables']['1']['data']
    header_dict = rf.clean_radial_header(radial_file_data['header'])

    # range cells in this radial
    range_dim = np.arange(df.RNGE.min(), df.RNGE.max()+float(header_dict['RangeResolutionKMeters']), float(header_dict['RangeResolutionKMeters']))

    # bearing_dim
    bearing_dim = np.arange(0, 360+header_dict['AngularResolution'], header_dict['AngularResolution'])

    # create radial grid
    [bearing, ranges] = np.meshgrid(bearing_dim, range_dim)

    # calculate lat/lons from origin, bearing, and ranges
    lonlat = re.findall(r"[-+]?\d*\.\d+|\d+", header_dict['Origin'])

    new_lat, new_lon = reckon(lonlat[1], lonlat[0], bearing, ranges)

    # create dictionary containing variables from dataframe in the shape of radial grid
    d = {key: np.tile(np.nan, bearing.shape) for key in df.keys()}

    # find grid indices from radial grid (bearing, ranges)
    x_ind, y_ind = gridded_index(bearing, ranges, df['BEAR'], df['RNGE'])

    for k,v in d.iteritems():
        v[(y_ind,x_ind)] = df[k]
        k = v

    ds = xr.Dataset()
    ds['lon'] = (('range', 'bearing'), new_lon)
    ds['lat'] = (('range', 'bearing'), new_lat)

    coords = ('time', 'range', 'bearing')
    # expand dimensions for time
    u = np.expand_dims(np.float32(d['VELU']), axis=0)
    ds['u'] = (coords, u)

    v = np.expand_dims(np.float32(d['VELV']), axis=0)
    ds['v'] = (coords, v)

    ds.coords['time'] = pd.date_range(header_dict['TimeStamp'], periods=1)
    ds.coords['bearing'] = bearing_dim
    ds.coords['range'] = range_dim


if __name__ == '__main__':
    radial_path = 'data/radials/AMAG/RDLi_AMAG_2017_05_02_0700.ruv'
    radial_file = os.path.join(os.path.dirname(__file__), '..', radial_path)
    main(radial_file)