import numpy as np
import xarray as xr

ds = xr.open_dataset('/Users/mikesmith/Documents/git/rucool/codar_processing/data/totals/nc/RU_MARA_20180301T000000Z.nc')


def curly_vectors(x, y, u, v, maxmag=None, minmag=None, numpoints=50, thin=1, linewidth=1, color='jet'):
    dx = np.abs(np.diff())