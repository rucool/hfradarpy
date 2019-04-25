import numpy as np
import numba


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