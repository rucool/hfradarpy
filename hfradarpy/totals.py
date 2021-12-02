import numpy as np
import xarray as xr
import pandas as pd
from hfradarpy.common import CTFParser
from hfradarpy.calc import gridded_index

import logging
logger = logging.getLogger(__name__)

# def concat_totals(radial_list):
#     """
#     This function takes a list of radial files. Loads them all separately using the Radial object and then combines
#     them along the time dimension using xarrays built-in concatenation routines.
#     :param radial_list: list of radial files that you want to concatenate
#     :return: radial files concatenated into an xarray dataset by range, bearing, and time
#     """
#
#     totals_dict = {}
#     for each in sorted(radial_list):
#         total = Total(each, multi_dimensional=True)
#         totals_dict[total.file_name] = total.ds
#
#     ds = xr.concat(totals_dict.values(), 'time')
#     return ds


class Totals(CTFParser):
    """
    Totals Subclass.

    This class should be used when loading a CODAR radial (.ruv) file. This class utilizes the generic LLUV class from
    ~/hfradar/ctf.py in order to load CODAR Radial files
    """
    def __init__(self, fname, replace_invalid=True, grid=False):

        CTFParser.__init__(self, fname)
        for key in self._tables.keys():
            table = self._tables[key]
            if 'LLUV' in table['TableType']:
                self.data = table['data']
            elif 'src' in table['TableType']:
                self.diagnostics_source = table['data']

        if replace_invalid:
            self.replace_invalid_values()

        # print()

        # if mask_over_land:
        #     land = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        #     land = land[land['continent'] == 'North America']
        #     # ocean = gpd.read_file('/Users/mikesmith/Downloads/ne_10m_ocean')
        #     self.data = gpd.GeoDataFrame(self.data, crs={'init': 'epsg:4326'}, geometry=[Point(xy) for xy in zip(self.data.LOND.values, self.data.LATD.values)])
        #
        #     # Join the geodataframe containing radial points with geodataframe containing leasing areas
        #     self.data = gpd.tools.sjoin(self.data, land, how='left')
        #
        #     # All data in the continent column that lies over water should be nan.
        #     self.data = self.data[keep][self.data['continent'].isna()]
        #     self.data = self.data.reset_index()

        if grid:
            self.to_multi_dimensional(grid)

    def to_multi_dimensional(self, grid_file):
        try:
            # load grid file
            grid = pd.read_csv(grid_file, sep=',', header=None, names=['lon', 'lat'], delim_whitespace=True)
            logging.debug('{} - Grid file successfully loaded '.format(grid_file))
        except Exception as err:
            logging.error('{} - {}. Grid file could not be loaded.'.format(grid_file, err))
            return

        lon = np.unique(grid['lon'].values.astype(np.float32))
        lat = np.unique(grid['lat'].values.astype(np.float32))
        [x, y] = np.meshgrid(lon, lat)

        logging.debug('{} - Gridding data to 2d grid'.format(grid_file))

        # convert 1d data into 2d gridded form. data_dict must be a dictionary.
        x_ind, y_ind = gridded_index(x, y, self.data.lon, self.data.lat)

        coords = ('time', 'range', 'bearing')

        # Intitialize empty xarray dataset
        ds = xr.Dataset()

    def file_type(self):
        """Return a string representing the type of file this is."""
        return 'totals'