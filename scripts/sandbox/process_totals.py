from codar_processing.src.radials import Radial
import numpy as np
import pandas as pd

sites = ['SEAB', 'BRAD', 'SPRK', 'BRNT', 'BRMR', 'RATH']
time = '2018_01_01_0000'
pattern_type = 'RDLi'

radials = []
for site in sites:
    radials.append(f'/Volumes/home/codaradm/data/radials/{site}/2018_01/{pattern_type}_{site}_{time}.ruv')

loaded = {}
for radial in radials:
    loaded[radial] = Radial(radial, to_xarray=False, mask_over_land=False)

grid_file = '../totals/grid_files/maracoos_grid_2km.txt'

# load csv file containing the grid
grid = pd.read_csv(grid_file, sep=',', header=None, names=['lon', 'lat'], delim_whitespace=True)
lon = np.unique(grid['lon'].values.astype(np.float32))
lat = np.unique(grid['lat'].values.astype(np.float32))
[x, y] = np.meshgrid(lon, lat)

x = x.ravel()
y = y.ravel()

sds = {}
for radial in loaded:
    print(radial)
#     for k = 1:size(TUV.LonLat, 1):
#     s = lonlat2dist(TUV.LonLat(k,:), R(j).LonLat' );
#
# sds
# {j, k} = find(s < sw(k));
# end
# end