from hfradarpy.radials import Radial
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import colors
from matplotlib.colors import Normalize
from matplotlib.colors import TwoSlopeNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from oceans.ocfis import uv2spdir, spdir2uv
import cmocean
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_ruv(radial_file, save_path=None, fname=None, speed_display = 'color', redblue=True, plotflag=None,  scale = 50, vlims=(-100, 100)):
    """
    Main function to plot radial files.

    Args:
        radial_file (str or Path): Path to radial file or a Radial object
        save_path (str or Path): Path to save figures
        fname (str): Output file name. If not specified, the radial object filename is used,  Defaults to None
        speed_display (str, optional): 'color' or 'arrowlength' to specify whether current speed is depicted by color or arrow length, Defaults to color
        redblue (bool, optional): If True, colorbar scheme is redblue, Defaults to True
        plotflag (str, optional): QARTOD QC test code, fail and suspect flags for that test will be highlighted, Defaults to None
        scale (int, optional): Scaling factor for drawing the vectors, Default = 50
        vlims (tuple, optional): Velocity limits for the colorbar, Default = (-100,100)
    """
    if not isinstance(radial_file, Radial):
        r = Radial(radial_file)
    else:
        r = radial_file

    if not r.is_valid():
        return
    if r._iscorrupt:
        return
    
    if fname == None:
        fname = r.file_name[0:-4]
    
    # Adjust some standard plotting settings to make them the size of a sheet of paper
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 12
    fig_size[1] = 8
    plt.rcParams["figure.figsize"] = fig_size

    # Set colors of the land.
    edgecolor = 'black'
    landcolor = 'tan'

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
    bf = 0.3  #degrees
    extent = [r.data.LOND.min()-bf , r.data.LOND.max()+bf, r.data.LATD.min()-bf, r.data.LATD.max()+bf]

    # Split out everything into seperate variables in order to pass them easier to the plotting functions
    time = r.time
    lon = r.data.LOND.to_numpy()
    lat = r.data.LATD.to_numpy()
    u = r.data.VELU.to_numpy()
    v = r.data.VELV.to_numpy()
    velocity = r.data.VELO.to_numpy()
    sitename = r.metadata['Site'][0:4]
    ptype = r.metadata['PatternType']

    # Mask nans just in case there are any
    u = ma.masked_invalid(u)
    v = ma.masked_invalid(v)

    # convert U and V component velocities to angle and speed
    angle, speed = uv2spdir(u, v)

    # convert angle and speed right back back to U and V component velocities,
    # Passing speed as an array of 1's allows for the normalizing of the arrow sizes
    # if we pass the correct
    u, v = spdir2uv(
        np.ones_like(speed),
        angle,
        deg=True
    )

    # Get the receiver location
    receiver_location = [float(x) for x in r.metadata['Origin'].split('  ')]
    receiver_location.reverse()

    # Intialize an empty subplot using cartopy
    fig, ax = plt.subplots(
        figsize=(11, 8),
        subplot_kw=dict(projection=ccrs.Mercator())
    )
    #plt.quiver(lon, lat, u, v, transform=ccrs.PlateCarree())
    plt.plot(receiver_location[0], receiver_location[1], 'o', markersize=10, markeredgecolor='black', color='red',
             transform=ccrs.PlateCarree())

    map_features(ax,extent,LAND,edgecolor,landcolor,state_lines)

    # The next lines specify the arrow shapes. You can customize this to your preference, usually by trial and error.
    # scale = 50
    headwidth = 2.5
    headlength = 4
    headaxislength = 4
    sub = 1

    # if user requested speed displayed as arrow length
    if speed_display == 'arrowlength':

        scale_units='width'
        width=0.005


        if not plotflag == None:
            fail = r.data[plotflag]==4
            suspect = r.data[plotflag]==3
            noteval = r.data[plotflag] == 2
            away = r.data.VELO > 0
            plt.quiver(lon, lat, u, v, transform=ccrs.PlateCarree(), scale=scale, scale_units=scale_units, width=width, color='lightpink')
            plt.quiver(lon[away], lat[away], u[away], v[away], transform=ccrs.PlateCarree(), scale=scale, scale_units=scale_units, width=width, color='lightblue')
            #plt.quiver(lon, lat, u, v, transform=ccrs.PlateCarree(), scale=scale, color='white')
            plt.quiver(lon[fail], lat[fail], u[fail], v[fail], transform=ccrs.PlateCarree(), scale=scale, scale_units=scale_units, width=width, color='red')
            plt.quiver(lon[suspect], lat[suspect], u[suspect], v[suspect], transform=ccrs.PlateCarree(), scale=scale, scale_units=scale_units, width=width, color='gold')
            plt.quiver(lon[noteval], lat[noteval], u[noteval], v[noteval], transform=ccrs.PlateCarree(), scale=scale, scale_units=scale_units, width=width, color='gray')
            plt.title(f'{sitename} {ptype} {plotflag}\nFail(red) Suspect(yellow) Not Evaluated(grey)\n{time}')
            plt.savefig(save_path + '/' + fname + '_' + plotflag)
            plt.close('all')
        elif redblue:
            away = r.data.VELO > 0
            plt.quiver(lon, lat, u, v, transform=ccrs.PlateCarree(), scale=scale, scale_units=scale_units, width=width, color='red')
            plt.quiver(lon[away], lat[away], u[away], v[away], transform=ccrs.PlateCarree(), scale=scale, scale_units=scale_units, width=width, color='blue')
            plt.title(f'{sitename} {ptype}\n{time}')
            plt.savefig(save_path + '/' + fname + '_rb')
            plt.close('all')
        else:
            plt.quiver(lon, lat, u, v, transform=ccrs.PlateCarree(), scale=scale, scale_units=scale_units, width=width, color='wheat')
            plt.title(f'{sitename} {ptype}\n{time}')
            plt.savefig(save_path + '/' + fname)
            plt.close('all')


    # if user requested speed displayed as color
    else:
        if not plotflag == None:

            test = r.data[plotflag]
            velocity[np.where(test >= 0)] = 1
            velocity[np.where(test == 4)] = -1

            color_clipped = np.clip(
                r.data.VELO[::sub],
                -1,
                1
            ).squeeze()
            offset = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
            cmap = colors.ListedColormap(['red', 'wheat'])

            plt.title(f'{sitename} {ptype} {plotflag} Fail\n{time}')

            qargs = dict(cmap=cmap, scale=scale, headwidth=headwidth, headlength=headlength,
                         headaxislength=headaxislength)
            qargs['transform'] = ccrs.PlateCarree()
            qargs['norm'] = offset

            # plot arrows over pcolor
            h = ax.quiver(
                lon[::sub],
                lat[::sub],
                u[::sub],
                v[::sub],
                color_clipped,
                **qargs
            )

            plt.savefig(save_path + '/' + fname + '_' + plotflag + '_fail')
            plt.close('all')

        elif redblue:

            cmap = 'bwr'
            # velocity_temp = velocity.where(velocity > 0, other=-1)  # Going away from radar
            velocity[np.where(velocity < 0)] = -1
            velocity[np.where(velocity >= 0)] = 1
            # We will create temporary variable of velocities that sets any velocity less than 0 to 1
            # color_clipped = velocity_temp.where(velocity < 0, other=1)  # Going towards radar
            # Define arrow colors. Limited by velocity_min and velocity_max
            color_clipped = np.clip(
                r.data.VELO[::sub],
                -1,
                1
            ).squeeze()

            offset = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

            plt.title(f'{sitename} {ptype}\n{time}')

            qargs = dict(cmap=cmap, scale=scale, headwidth=headwidth, headlength=headlength,
                         headaxislength=headaxislength)
            qargs['transform'] = ccrs.PlateCarree()
            qargs['norm'] = offset

            # plot arrows over pcolor
            h = ax.quiver(
                lon[::sub],
                lat[::sub],
                u[::sub],
                v[::sub],
                color_clipped,
                **qargs
            )
            # map_features(ax, extent, LAND, edgecolor, landcolor, state_lines)
            plt.savefig(save_path + '/' + fname + '_rb')
            plt.close('all')

        else:

            plt.title(f'{sitename} {ptype}\n{time}')

            cmap = cmocean.cm.balance
            # Colorbar options
            velocity_min = vlims[0]  # The minimum speed that should be displayed on the colorbar
            velocity_max = vlims[1]  # The maximum speed that should be displayed on the colorbar
            cbar_step = 10  # The step between each colorbar tick

            offset = Normalize(vmin=velocity_min, vmax=velocity_max, clip=True)

            # Define arrow colors. Limited by velocity_min and velocity_max
            color_clipped = np.clip(
                r.data.VELO[::sub],
                velocity_min,
                velocity_max
            ).squeeze()

            ticks = np.append(np.arange(velocity_min, velocity_max, cbar_step), velocity_max)

            qargs = dict(cmap=cmap, scale=scale, headwidth=headwidth, headlength=headlength,
                         headaxislength=headaxislength)
            qargs['transform'] = ccrs.PlateCarree()
            qargs['norm'] = offset

            # plot arrows over pcolor
            h = ax.quiver(
                lon[::sub],
                lat[::sub],
                u[::sub],
                v[::sub],
                color_clipped,
                **qargs
            )

            # map_features(ax, extent, LAND, edgecolor, landcolor, state_lines)

            # generate colorbar
            divider = make_axes_locatable(ax)
            cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
            fig.add_axes(cax)

            cb = plt.colorbar(h, cax=cax, ticks=ticks)
            cb.ax.set_yticklabels([f'{s:d}' for s in ticks])
            cb.set_label('cm/s')

            plt.savefig(save_path + '/' + fname)
            plt.close('all')







# Create a re-usable function for map features that we can pass an axes to.
def map_features(ax,extent,LAND, edgecolor, landcolor, state_lines):
    # Axes properties and features
    ax.set_extent(extent)
    ax.add_feature(LAND, edgecolor=edgecolor, facecolor=landcolor)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(state_lines, zorder=11, edgecolor=edgecolor)

    # Gridlines and grid labels
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=.5,
        color='black',
        alpha=0.25,
        linestyle='--'
    )

    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    gl.xlocator = mticker.MaxNLocator(integer=True)
    gl.ylocator = mticker.MaxNLocator(integer=True)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ax.tick_params(which='major',
                   direction='out',
                   bottom=True, top=True,
                   labelbottom=True, labeltop=False,
                   left=True, right=True,
                   labelleft=True, labelright=False,
                   length=5, width=2)

    ax.tick_params(which='minor',
                   direction='out',
                   bottom=True, top=True,
                   labelbottom=True, labeltop=False,
                   left=True, right=True,
                   labelleft=True, labelright=False,
                   width=1)
