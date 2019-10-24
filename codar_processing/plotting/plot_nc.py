import xarray as xr
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from oceans.ocfis import uv2spdir, spdir2uv
from mpl_toolkits.axes_grid1 import make_axes_locatable
from codar_processing.src.common import create_dir

LAND = cfeature.NaturalEarthFeature(
    'physical', 'land', '10m',
    edgecolor='face',
    facecolor='tan'
)

state_lines = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none'
)


def plot_totals(nc_file, save_dir, *,
                sub=2, velocity_min=None, velocity_max=None, markers=None, title='HF Radar'):
    """
    param markers:  a list of 3-tuple/lists containng [lon, lat, marker kwargs] as should be
                    passed into ax.plot()
                    eg. [
                            [-74.6, 38.5, dict(marker='o', markersize=8, color='r')],
                            [-70.1, 35.2, dict(marker='o', markersize=8, color='b')]
                        ]
    """
    markers = markers or []

    with xr.open_dataset(str(nc_file)) as ds:
        fig = plt.figure()

        tds = ds.squeeze()

        u = tds['u'].data
        v = tds['v'].data

        time = str(ds.time.values[0])
        lon = tds.coords['lon'].data
        lat = tds.coords['lat'].data

    u = ma.masked_invalid(u)
    v = ma.masked_invalid(v)

    angle, speed = uv2spdir(u, v)
    us, vs = spdir2uv(
        np.ones_like(speed),
        angle,
        deg=True
    )

    lons, lats = np.meshgrid(lon, lat)

    velocity_min = velocity_min or 0
    velocity_max = velocity_max or np.nanmax(speed) or 15
    
    speed_clipped = np.clip(
        speed[::sub, ::sub],
        velocity_min,
        velocity_max
    ).squeeze()

    fig, ax = plt.subplots(
        figsize=(11, 8),
        subplot_kw=dict(projection=ccrs.PlateCarree())
    )

    # Plot title
    plt.title('{}\n{}'.format(title, time))

    # plot arrows over pcolor
    h = ax.quiver(
        lons[::sub, ::sub],
        lats[::sub, ::sub],
        us[::sub, ::sub],
        vs[::sub, ::sub],
        speed_clipped,
        cmap='jet',
        scale=60
    )

    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
    fig.add_axes(cax)

    # generate colorbar
    ticks = np.linspace(velocity_min, velocity_max, 5)
    cb = plt.colorbar(h, cax=cax, ticks=ticks)
    cb.ax.set_yticklabels([ f'{s:.2f}' for s in ticks ])
    cb.set_label('cm/s')

    for m in markers:
        ax.plot(m[0], m[1], **m[2])

    # Gridlines and grid labels
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=1,
        color='black',
        alpha=0.5,
        linestyle='--'
    )
    gl.xlabels_top = gl.ylabels_right = False
    gl.xlabel_style = {'size': 10, 'color': 'gray'}
    gl.ylabel_style = {'size': 10, 'color': 'gray'}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Axes properties and features
    ax.set_extent([
        lon.min() - 1,
        lon.max() + 1,
        lat.min() - 1,
        lat.max() + 1
    ])
    ax.add_feature(LAND, zorder=0, edgecolor='black')
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(state_lines, edgecolor='black')

    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 12
    fig_size[1] = 8.5
    plt.rcParams["figure.figsize"] = fig_size

    sname = '{}.png'.format(os.path.basename(str(nc_file)))
    save_name = os.path.join(save_dir, sname)
    create_dir(save_dir)

    resoluton = 150  # plot resolution in DPI
    plt.savefig(save_name, dpi=resoluton)
    plt.close('all')
