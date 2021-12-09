import numba
import numpy as np
import xarray as xr
from scipy import spatial


@numba.jit
def bool2int(x):
    """

    :param x:
    :return:
    """
    y = 0
    for i,j in enumerate(x):
        y += j << i
    return y


@numba.jit(nopython=True)
def calc_index(x_ind, y_ind, X, Y, x, y):
    for i, line in enumerate(x):
        x_ind[i] = np.argmin(np.abs(X[1, :] - x[i]))
        y_ind[i] = np.argmin(np.abs(Y[:, 1] - y[i]))
    return x_ind, y_ind


def gridded_index(X, Y, x, y, flag=np.nan):
    """
    This function gets the multidimensional index of 1d grid onto a 2d grid without interpolation. It calculates the
    index based on the difference between two points.
    :param X: x grid of values (M x N). Must be a numpy.ndarray
    :param Y: y grid of values (M x N). Must be a numpy.ndarray
    :param x: n vector of x values. Must be a numpy.ndarray
    :param y: n vector of y values. Must be a numpy.ndarray
    :param flag: value to use for missing data values of grid. Default is np.nan
    :return:
    """
    # get mapping index
    x_ind = np.tile(flag, x.size).astype(np.int)
    y_ind = np.tile(flag, y.size).astype(np.int)

    # Roll index calculation into other function so we can use numba for speedups
    x_ind, y_ind = calc_index(x_ind, y_ind, X, Y, x, y)
    return x_ind, y_ind


@numba.jit
def lond_jit(lon, lat, dist, bearing, EARTH_RADIUS):
    next_longitude = lon + (np.arctan2(np.sin(bearing) * np.sin(dist/EARTH_RADIUS) * np.cos(lat),
                                             np.cos(dist/EARTH_RADIUS) - np.sin(lat) * np.sin(lat)))
    return next_longitude


def reckon(lat, lon, bearing, distance):
    """
    Point at specified azimuth, range on sphere or ellipsoid
    :param lat: origin latitude (decimal degrees)
    :param lon: origin longitude (decimal degrees)
    :param bearing: np.array of bearings
    :param distance: np.array of ranges
    :return: latitude, longitude of calculated coordinates from o
    """
    EARTH_RADIUS = 6371.00

    # convert origin lat/lons into radians
    latitude = np.radians(lat)
    longitute = np.radians(lon)
    bearing = np.radians(bearing)

    # calculate distance latitude
    next_latitude = np.arcsin(np.sin(latitude) *
                    np.cos(distance/EARTH_RADIUS) +
                    np.cos(latitude) *
                    np.sin(distance/EARTH_RADIUS) *
                    np.cos(bearing))

    # calculate distance longitude. For some reason, lon is faster when calculated with numba function compared to lat. So we split it out into a new function.....
    next_longitude = lond_jit(longitute, latitude, distance, bearing, EARTH_RADIUS)

    # convert points into decimal degrees
    new_lat = np.degrees(next_latitude)
    new_lon = np.degrees(next_longitude)

    return new_lat, new_lon


class KDTreeIndex():
    """ A KD-tree implementation for fast point lookup on a 2D grid

    Keyword arguments:
    dataset -- a xarray DataArray containing lat/lon coordinates
               (named 'lat' and 'lon' respectively)

    from: https://notes.stefanomattia.net/2017/12/12/The-quest-to-find-the-closest-ground-pixel/
    """

    def transform_coordinates(self, coords):
        """ Transform coordinates from geodetic to cartesian

        Keyword arguments:
        coords - a set of lan/lon coordinates (e.g. a tuple or
                 an array of tuples)
        """
        # WGS 84 reference coordinate system parameters
        A = 6378.137  # major axis [km]
        E2 = 6.69437999014e-3  # eccentricity squared

        coords = np.asarray(coords).astype(np.float)

        # is coords a tuple? Convert it to an one-element array of tuples
        if coords.ndim == 1:
            coords = np.array([coords])

        # convert to radiants
        lat_rad = np.radians(coords[:, 0])
        lon_rad = np.radians(coords[:, 1])

        # convert to cartesian coordinates
        r_n = A / (np.sqrt(1 - E2 * (np.sin(lat_rad) ** 2)))
        x = r_n * np.cos(lat_rad) * np.cos(lon_rad)
        y = r_n * np.cos(lat_rad) * np.sin(lon_rad)
        z = r_n * (1 - E2) * np.sin(lat_rad)

        return np.column_stack((x, y, z))

    def __init__(self, dataset):
        # store original dataset shape
        self.shape = dataset.shape

        # reshape and stack coordinates
        coords = np.column_stack((dataset.lat.values.ravel(),
                                  dataset.lon.values.ravel()))

        # construct KD-tree
        self.tree = spatial.cKDTree(self.transform_coordinates(coords))

    def query(self, point):
        """ Query the kd-tree for nearest neighbour.

        Keyword arguments:
        point -- a (lat, lon) tuple or array of tuples
        """
        _, index = self.tree.query(self.transform_coordinates(point))

        # regrid to 2D grid
        index = np.unravel_index(index, self.shape)

        # return DataArray indexers
        return xr.DataArray(index[0], dims='location'), xr.DataArray(index[1], dims='location')

    def query_ball_point(self, point, radius):
        """ Query the kd-tree for all point within distance
        radius of point(s) x

        Keyword arguments:
        point -- a (lat, lon) tuple or array of tuples
        radius -- the search radius (km)
        """

        index = self.tree.query_ball_point(self.transform_coordinates(point),
                                           radius)

        # regrid to 2D grid
        index = np.unravel_index(index[0], self.shape)

        # return DataArray indexers
        return xr.DataArray(index[0], dims='location'), xr.DataArray(index[1], dims='location')