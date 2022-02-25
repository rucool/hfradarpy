import numpy as np
from pyproj import Geod

import logging
logger = logging.getLogger(__name__)

geodesic = Geod(ellps='WGS84') #define the coordinate system. WGS84 is the standard used by GPS.

def reckon(origin_lon, origin_lat, forward_azimuth, range):
    """ Calculate lon, lat of a point from a specified azimuth, range on sphere or ellipsoid
    Helper function for pyproj.Geod forward transformation

    Args:
        origin_lon (array, numpy.ndarray, list, tuple, or scalar): Longitude(s) of initial point(s)
        origin_lat (array, numpy.ndarray, list, tuple, or scalar): Latitude(s) of initial point(s)
        forward_azimuth (array, numpy.ndarray, list, tuple, or scalar): Azimuth/bearing(s) of the terminus point relative to the initial point(s)
        range (array, numpy.ndarray, list, tuple, or scalar): Distance(s) between initial and terminus point(s) in kilometers

    Returns:
        array, numpy.ndarray, list, tuple, or scalar: Longitude(s) of terminus point(s)
        array, numpy.ndarray, list, tuple, or scalar: Latitude(s) of terminus point(s)
        array, numpy.ndarray, list, tuple, or scalar: Backwards azimuth(s) of terminus point(s)
    """
    terminus_lon, terminus_lat, _ = geodesic.fwd(origin_lon, origin_lat, forward_azimuth, range*1000)
    logging.info(f'Calculating longitude and latitude of terminus points from origin based off of a bearing of {forward_azimuth} (degrees) and range {range} (km).') 
    logging.info(f'Origin - Lon (degrees): {origin_lon} degrees, Lat (degrees): {origin_lat}') 
    logging.info(f'Terminus - Lon(s) (degrees): {terminus_lon} degrees, Lat(s) (degrees): {terminus_lat}') 
    return terminus_lon, terminus_lat

def inverse_transformation(lons1, lats1, lons2, lats2):
    """Inverse computation of bearing and distance given the latitudes and longitudes of an initial and terminus point.

    Args:
        lons1 (array, numpy.ndarray, list, tuple, or scalar): Longitude(s) of initial point(s)
        lats1 (array, numpy.ndarray, list, tuple, or scalar): Latitude(s) of initial point(s)
        lons2 (array, numpy.ndarray, list, tuple, or scalar): Longitude(s) of terminus point(s)
        lats2 (array, numpy.ndarray, list, tuple, or scalar): Latitude(s) of terminus point(s)

    Returns:
       array, numpy.ndarray, list, tuple, or scalar: Forward azimuth(s)
       array, numpy.ndarray, list, tuple, or scalar: Back azimuth(s)
       array, numpy.ndarray, list, tuple, or scalar: Distance(s) between initial and terminus point(s) in kilometersmeters
    """
    # Inverse transformation using pyproj
    # Determine forward and back azimuths, plus distances between initial points and terminus points.
    forward_azimuth, back_azimuth, distance = geodesic.inv(lons1, lats1, lons2, lats2)

    distance = distance/1000 # Lets stick with kilometers as the output since RUV ranges are in kilometers
        
    logging.info(f'Inversely calculating azimuth and range of initial point(s) to terminus point(s)')
    logging.info(f'Origin - Lon (degrees): {lons1} degrees, Lat (degrees): {lats1}') 
    logging.info(f'Terminus - Lon(s) (degrees): {lons2} degrees, Lat(s) (degrees): {lats2}') 
    logging.info(f'Forward azimuth (CWN): {forward_azimuth} degrees, Back bearing (CWN): {back_azimuth} degrees, Range: {distance} kilometers') 
    return forward_azimuth, back_azimuth, distance

def spdir2uv(speed, direction, degrees=False):
    """Calculate u and v velocities from speed and direction
    
    Args:
        speed (numpy.ndarray): Speed (meters/second)
        direction (nump.ndarray): Direction (degrees)
        degrees (bool, optional): True if data is in degrees. Defaults to False.
        
    Returns:
        (numpy.ndarray): u - Eastward Seawater Velocity (meters/second)
        (numpy.ndarray): v - Northward Seawater Velocity (meters/second)
    """
    if degrees:
        # Calculating u and v requires the data to be in radians
        direction = np.deg2rad(direction)

    u = speed * np.sin(direction) # Calculate east-west component
    v = speed * np.cos(direction) # Calculate north-south component
    return u, v

def uv2spdir(u, v, mag=0, rot=0):
    """Calculate speed (m/s) and direction (degrees) from u (eastward) and v (northward) velocity (meters/second) components
    Converts rectangular to polar coordinate, geographic convention
    Allows for rotation and magnetic declination correction.

    Args:
        u (numpy.ndarray): u - Eastward Seawater Velocity (meters/second)
        v (numpy.ndarray): v - Northward Seawater Velocity (meters/second)
        mag (int, optional): Magnetic Correction (degrees). Defaults to 0.
        rot (int, optional): Angle for rotation (degrees). Defaults to 0.

    Returns:
        numpy.ndarray: Speed (meters/second)
        numpy.ndarray: Direction (degrees) is traveling towards
    """
    u, v, mag, rot = list(map(np.asarray, (u, v, mag, rot)))

    vector = u + 1j * v
    speed = np.abs(vector) #np.hypot(u, v) also works
    direction = np.angle(vector, deg=True)
    direction = direction - mag + rot
    direction = np.mod(90.0 - direction, 360.0)  # Zero is North.

    return speed, direction

# @numba.jit
# def bool2int(x):
#     """
#     Convert a boolean to an integer
    
#     Keyword arguments:
#     argument -- description
#     Return: return_description
#     """
#     y = 0
#     for i,j in enumerate(x):
#         y += j << i
#     return y

# @numba.jit(nopython=True)
# def calc_index(x_ind, y_ind, X, Y, x, y):
#     """
#     Calculate index

#     Args:
#         x_ind (_type_): x index 
#         y_ind (_type_): y index
#         X (numpy.ndarray): x grid of values (M x N)
#         Y (numpy.ndarray): y grid of values (M x N)
#         x (numpy.ndarray): x grid of values (M x N)
#         y (numpy.ndarray): y grid of values (M x N)

#     Returns:
#         _type_: _description_
#         _type_: _description_
#     """
#     for i, line in enumerate(x):
#         x_ind[i] = np.argmin(np.abs(X[1, :] - x[i]))
#         y_ind[i] = np.argmin(np.abs(Y[:, 1] - y[i]))
#     return x_ind, y_ind

# def gridded_index(X, Y, x, y, flag=np.nan):
#     """
#     This function gets the multidimensional index of 1d grid onto a 2d grid without interpolation. It calculates the
#     index based on the difference between two points.

#     Args:
#         X (numpy.ndarray): x grid of values (M x N)
#         Y (numpy.ndarray): y grid of values (M x N)
#         x (numpy.ndarray): x grid of values (M x N)
#         y (numpy.ndarray): y grid of values (M x N)
#         flag (_type_, optional): value to use for missing data values of grid.. Defaults to np.nan.

#     Returns:
#         _type_: _description_
#         _type_: _description_

#     """
#     # get mapping index
#     x_ind = np.tile(flag, x.size).astype(np.int)
#     y_ind = np.tile(flag, y.size).astype(np.int)

#     # Roll index calculation into other function so we can use numba for speedups
#     x_ind, y_ind = calc_index(x_ind, y_ind, X, Y, x, y)
#     return x_ind, y_ind