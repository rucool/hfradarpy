import numpy as np
from numpy.matlib import tile


def bool2int(x):
    y = 0
    for i,j in enumerate(x):
        y += j << i
    return y


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
    x_ind = tile(flag, x.size).astype(np.int)
    y_ind = tile(flag, y.size).astype(np.int)

    for i, line in enumerate(x):
        x_ind[i] = np.argmin(np.abs(X[1, :] - x[i]))
        y_ind[i] = np.argmin(np.abs(Y[:, 1] - y[i]))

    return x_ind, y_ind


def reckon(lat, lon, bearing, distance):
    """
    Point at specified azimuth, range on sphere or ellipsoid
    :param lat: origin latitude (decimal degrees)
    :param lon: origin longitude (decimal degrees)
    :param bearing: np.array of bearings
    :param distance: np.array of ranges
    :return: latitude, longitude of calculated coordinates from o
    """
    # standard earth radius

    EARTH_RADIUS = 6371.00

    # convert Latitude and Longitude
    # into radians for calculation
    latitude = np.radians(lat)
    longitute = np.radians(lon)

    # calculate next latitude
    next_latitude = np.arcsin(np.sin(latitude) *
                    np.cos(distance/EARTH_RADIUS) +
                    np.cos(latitude) *
                    np.sin(distance/EARTH_RADIUS) *
                    np.cos(np.radians(bearing)))

    # calculate next longitude
    next_longitude = longitute + (np.arctan2(np.sin(np.radians(bearing)) * np.sin(distance/EARTH_RADIUS) * np.cos(latitude),
                                             np.cos(distance/EARTH_RADIUS) - np.sin(latitude) * np.sin(next_latitude)))

    # convert points into decimal degrees
    new_lat = np.degrees(next_latitude)
    new_lon = np.degrees(next_longitude)

    return new_lat, new_lon


