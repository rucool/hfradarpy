from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from hfradarpy.common import create_dir
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature

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


def plot_common(time, lon, lat, u, v, *,
                output_file=None,
                color_clipped=None,
                meshgrid=True,
                sub=2,
                markers=None,
                offset=None,
                extent=None,
                lon_ticks=None,
                lat_ticks=None,
                cmap='jet',
                colorbar=True,
                ticks=None,
                title='HF Radar',
                scale=120, headwidth=2.5, headlength=4, headaxislength=4):
    """
    param markers:  a list of 3-tuple/lists containng [lon, lat, marker kwargs] as should be
                    passed into ax.plot()
                    eg. [
                            [-74.6, 38.5, dict(marker='o', markersize=8, color='r')],
                            [-70.1, 35.2, dict(marker='o', markersize=8, color='b')]
                        ]
    """
    markers = markers or []

    if meshgrid is True:
        lons, lats = np.meshgrid(lon, lat)
    else:
        lons, lats = lon, lat

    extent = extent or [
        lon.min() - 1,
        lon.max() + 1,
        lat.min() - 1,
        lat.max() + 1
    ]

    fig, ax = plt.subplots(
        figsize=(11, 8),
        subplot_kw=dict(projection=ccrs.Mercator())
    )

    # Plot title
    plt.title('{} - {}'.format(title, pd.to_datetime(time).strftime('%Y-%m-%d %H:%M:%S GMT')))

    qargs = dict(cmap=cmap, scale=scale, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength)
    qargs['transform'] = ccrs.PlateCarree()
    qargs['norm'] = offset

    # plot arrows over pcolor
    h = ax.quiver(
        lons[::sub],
        lats[::sub],
        u[::sub],
        v[::sub],
        color_clipped,
        **qargs
    )

    # generate colorbar
    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
        fig.add_axes(cax)

        cb = plt.colorbar(h, cax=cax, ticks=ticks)
        cb.ax.set_yticklabels([f'{s:d}' for s in ticks])
        cb.set_label('cm/s')

    for m in markers:
        ax.plot(m[0], m[1], transform=qargs['transform'], **m[2])

    # Gridlines and grid labels
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=.5,
        color='black',
        alpha=0.25,
        linestyle='--'
    )
    gl.xlabels_top = gl.ylabels_right = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

    if lon_ticks:
        gl.xlocator = mticker.FixedLocator(lon_ticks)
    if lat_ticks:
        gl.ylocator = mticker.FixedLocator(lat_ticks)

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Axes properties and features
    ax.set_extent(extent)
    ax.add_feature(LAND, zorder=0, edgecolor='black')
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(state_lines, edgecolor='black')

    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 12
    fig_size[1] = 8.5
    plt.rcParams["figure.figsize"] = fig_size

    if output_file is not None:
        create_dir(str(Path(output_file).parent))
        resoluton = 300  # plot resolution in DPI
        plt.savefig(output_file, dpi=resoluton, bbox_inches='tight', pad_inches=0.1)
        plt.close('all')
    else:
        return plt
