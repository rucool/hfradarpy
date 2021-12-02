import xarray as xr
import numpy as np
import numpy.ma as ma
from oceans.ocfis import uv2spdir, spdir2uv
from hfradarpy.plot.common import plot_common
import logging
from matplotlib import colors
from matplotlib.colors import TwoSlopeNorm, Normalize
import cmocean

logger = logging.getLogger(__name__)


def plot_radials(dataset, *,
                 plot_type='velocity',
                 output_file=None,
                 extent=None, lon_ticks=None, lat_ticks=None,
                 scale=True,
                 sub=2,
                 cbar_step=10,
                 velocity_min=None, velocity_max=None,
                 markers=None,
                 prim_filter=False,
                 title='HF Radar'):
    """
    param dataset:  a file-path to an xarray compatible object or an xarray Dataset object
    """
    try:
        ds = xr.open_dataset(dataset)
        closing = ds.close
    except AttributeError:
        if isinstance(dataset, xr.Dataset):
            ds = dataset
            closing = int  # dummy func to close nothing
        else:
            raise

    tds = ds.squeeze()

    if prim_filter:
        if 'primary_flag_qc' in list(tds.keys()):
            tds = tds.where(tds.primary_flag_qc == 1).squeeze()
            title = title + ' - QC'
        else:
            logging.warning('PRIM flag not found. Bypassing quality control filters')

    time = ds.time.values[0]
    lon = tds.coords['lon'].data
    lat = tds.coords['lat'].data
    u = tds['u'].data
    v = tds['v'].data

    u = ma.masked_invalid(u)
    v = ma.masked_invalid(v)

    if scale:
        angle, speed = uv2spdir(u, v)  # convert u/v to angle and speed
        u, v = spdir2uv(  # convert angle and speed back to u/v, normalizing the arrow sizes
            np.ones_like(speed),
            angle,
            deg=True
        )

    velocity_min = velocity_min or -40
    velocity_max = velocity_max or 40

    kwargs = dict(extent=extent, lon_ticks=lon_ticks, lat_ticks=lat_ticks,
                  output_file=output_file,
                  sub=sub,
                  markers=markers,
                  title=title,
                  meshgrid=False,
                  scale=120, headwidth=2.5, headlength=4, headaxislength=4)

    if 'velocity' in plot_type:
        """
        Velocity displays the direction and magnitude of the radials
        """
        kwargs['colorbar'] = True
        kwargs['cmap'] = cmocean.cm.balance

        # Define arrow colors. Limited by velocity_min and velocity_max
        kwargs['color_clipped'] = np.clip(
            tds.velocity.data[::sub],
            velocity_min,
            velocity_max
        ).squeeze()
        kwargs['offset'] = Normalize(vmin=velocity_min, vmax=velocity_max, clip=True)
        kwargs['ticks'] = np.append(np.arange(velocity_min, velocity_max, cbar_step), velocity_max)
    elif 'motion' in plot_type:
        """
        Motion displays the direction (towards or away) from radar
        """
        kwargs['colorbar'] = False
        kwargs['cmap'] = 'bwr'

        velocity = tds.velocity
        velocity_temp = velocity.where(velocity > 0, other=-1)  # Going away from radar
        kwargs['color_clipped'] = velocity_temp.where(velocity < 0, other=1).data  # Going towards radar
        kwargs['offset'] = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    elif 'qc_pass_fail' in plot_type:
        kwargs['colorbar'] = False
        kwargs['cmap'] = colors.ListedColormap(['limegreen', 'red'])

        if prim_filter:
            tds = ds.squeeze()
        kwargs['color_clipped'] = tds.primary_flag_qc.where(tds.primary_flag_qc != 1, other=-1).data  # PRIM == 1 where vectors pass qc
        kwargs['offset'] = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

    closing()
    plot_common(time, lon, lat, u, v, **kwargs)


def plot_totals(dataset, *,
                output_file=None,
                extent=None,
                scale=True,
                sub=2,
                cbar_step=20,
                velocity_min=None, velocity_max=None,
                markers=None,
                title='HF Radar'):
    """
    param dataset:  a file-path to an xarray compatible object or an xarray Dataset object
    """
    try:
        ds = xr.open_dataset(dataset)
        closing = ds.close
    except AttributeError:
        if isinstance(dataset, xr.Dataset):
            ds = dataset
            closing = int  # dummy func to close nothing
        else:
            raise

    tds = ds.squeeze()
    u = tds['u'].data
    v = tds['v'].data
    lon = tds.coords['lon'].data
    lat = tds.coords['lat'].data

    try:
        time = str(ds.time.values[0])
    except IndexError:
        time = ds.time.values

    closing()

    if scale:
        angle, speed = uv2spdir(u, v)  # convert u/v to angle and speed
        u, v = spdir2uv(  # convert angle and speed back to u/v, normalizing the arrow sizes
            np.ones_like(speed),
            angle,
            deg=True
        )

    velocity_min = velocity_min or np.int32(np.nanmin(speed)) or 0
    velocity_max = velocity_max or np.int32(np.nanmax(speed)) or 15

    kwargs = dict(output_file=output_file,
                  sub=sub,
                  markers=markers,
                  title=title,
                  meshgrid=True,
                  scale=60,
                  headwidth=3,
                  headlength=5,
                  headaxislength=4.5)

    kwargs['color_clipped'] = np.clip(
        speed[::sub],
        velocity_min,
        velocity_max
    ).squeeze()
    kwargs['offset'] = Normalize(vmin=velocity_min, vmax=velocity_max, clip=True)
    kwargs['ticks'] = np.append(np.arange(velocity_min, velocity_max, cbar_step), velocity_max)
    kwargs['extent'] = extent

    plt = plot_common(
        time, lon, lat, u, v,
        **kwargs
    )

    if output_file is None:
        return plt
