import os
import xarray as xr
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
from oceans.ocfis import uv2spdir, spdir2uv
from mpl_toolkits.axes_grid1 import make_axes_locatable
from hfradarpy.common import create_dir
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

f = '../data/totals/oi/nc/hourly/RU_MARA_20190101T000000Z.nc'
save_dir = '../data/plots/totals/'
title_str = 'MARACOOS'
velocity_min = 0
velocity_max = 15
resoluton = 150  # plot resolution in DPI
sub = 2  # subset data by every n sub

ds = xr.open_dataset(f)

LAND = cfeature.NaturalEarthFeature(
    'physical', 'land', '10m',
    edgecolor='face',
    facecolor='tan'
)

state_lines = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')

fig = plt.figure()
tds = ds.squeeze()

u = tds['u'].data
v = tds['v'].data

lon = tds.coords['lon'].data
lat = tds.coords['lat'].data
# time = tds.coords['time'].data

u = ma.masked_invalid(u)
v = ma.masked_invalid(v)

angle, speed = uv2spdir(u, v)
us, vs = spdir2uv(np.ones_like(speed), angle, deg=True)

lons, lats = np.meshgrid(lon, lat)

speed_clipped = np.clip(speed[::sub, ::sub], velocity_min, velocity_max).squeeze()

fig, ax = plt.subplots(figsize=(11, 8),
                         subplot_kw=dict(projection=ccrs.PlateCarree()))

# Plot title
plt.title('{}\n{}'.format(title_str, str(ds.time.values[0])))

# plot arrows over pcolor
h = ax.quiver(lons[::sub, ::sub], lats[::sub, ::sub],
              us[::sub, ::sub], vs[::sub, ::sub],
              speed_clipped,
              cmap='jet',
              scale=60)

divider = make_axes_locatable(ax)
cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
fig.add_axes(cax)

# generate colorbar
cb = plt.colorbar(h, cax=cax, ticks=[0,5,10,15])
cb.set_label('cm/s')
# cb.ax.set_yticklabels(['0', '5', '10', '15'])  # vertically oriented colorbar


ax.plot(-74.666667, 38.533333, marker='o', markersize=8, color='crimson')

# Gridlines and grid labels
gl = ax.gridlines(draw_labels=True,
                  linewidth=1,
                  color='black',
                  alpha=0.5, linestyle='--')
gl.xlabels_top = gl.ylabels_right = False
gl.xlabel_style = {'size': 15, 'color': 'gray'}
gl.ylabel_style = {'size': 15, 'color': 'gray'}
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

# Axes properties and features
ax.set_extent([-76.5, -68.5, 35, 42.75])
ax.add_feature(LAND, zorder=0, edgecolor='black')
ax.add_feature(cfeature.LAKES)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(state_lines, edgecolor='black')

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 8.5
plt.rcParams["figure.figsize"] = fig_size

sname = '{}.png'.format(os.path.basename(f))
save_name = os.path.join(save_dir, sname)
create_dir(save_dir)

plt.savefig(save_name, dpi=resoluton)
plt.close('all')
