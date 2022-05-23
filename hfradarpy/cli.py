import click
import glob
import os
from hfradarpy import __version__
from hfradarpy.calc import inverse_transformation as inv
from hfradarpy.common import list_to_dataframe
from hfradarpy.radials import concat
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path

fmt = "%Y-%m-%dT%H%M%S"

@click.group()
def main(args=None):
    """Console script for HFRadarPy python package."""
    click.echo("Toolbox to read in High Frequency Radar (HFR) files written in CODAR Tabular Format (CTF).")
    click.echo(f"hfradarpy: {__version__}")

    return 0
    # pass


@click.command()
@click.option("--radial_file", "-r", default=None, type=str, help="Name of a radial file.")
@click.option("--path_data", "-pd", default=os.getcwd(), type=str, help="Path containing radial files.")
@click.option("--path_save", "-ps", default=os.getcwd(), type=str, help="Path to save output file")
@click.option("--time_start", "-t0", type=str, default=None, help="Start Time (yyyy-mm-dd HH:MM:SS)")
@click.option("--time_end", "-t1", type=str, default=None, help="End Time (yyyy-mm-dd HH:MM:SS)")
@click.option("--prefix", "-p", type=str, default=None, help="Filter files based on a prefix such as RDLm or RDLi")
@click.option("--outname", "-f", type=str, help="Export filename.")
@click.option("--speed_display", "-sd", default="color", type=str,
              help="Method of displaying speed: color, arrowlength")
@click.option("--redblue", "-rb", default=True, type=bool, help="If True, use a red/blue color palette.")
@click.option("--plotflag", "-pf", default=None, type=str,
              help="Name of flag column, current vectors with flag values > 1 will be highlighted in the map")
@click.option("--scale", "-sc", default=50, type=float, help="Scale factor for vectors.")
@click.option("--limits", "-lim", nargs=2, type=int, default=None, help="Velocity limits for colorbar. Example --limits -100 100")

# https://click.palletsprojects.com/en/8.1.x/options/#tuples-as-multi-value-options
def plotruv(radial_file, path_data, path_save, time_start, time_end, prefix, outname, speed_display, redblue, plotflag, scale, limits):
    try:
        from hfradarpy.plot.plot_ruv import plot_ruv
    except ImportError:
        click.echo("Matplotlib and cartopy not installed. Please install.")
        return

    click.echo("Executing hfradarpy.cli.plot_ruv")
    click.echo("")

    os.makedirs(path_save, exist_ok=True)

    if radial_file != None:
        plot_ruv(radial_file, save_path=path_save, fname=outname, speed_display=speed_display, redblue=redblue, plotflag=plotflag, scale=scale, vlims=limits)
    else:
        click.echo(f"Checking for radial data in {path_data}")

        files = sorted(glob.glob(os.path.join(path_data, prefix + '*.ruv')))
        click.echo(f"{len(files)} radial files found")
        df = list_to_dataframe(files)
        click.echo(f"Data from {df.iloc[0].name} to {df.iloc[-1].name}")

        if time_start or time_end:
            if time_start:
                time_start = pd.to_datetime(time_start).tz_localize(None)
            else:
                time_start = df.iloc[0].name

            if time_end:
                time_end = pd.to_datetime(time_end).tz_localize(None)
            else:
                time_end = df.iloc[-1].name

            click.echo(f"Subsetting from {time_start} to {time_end}")
            click.echo("")
            df = df[time_start:time_end]

        filelist = df['file'].tolist()

        for f in filelist:
            print(f)
            plot_ruv(f, save_path=path_save, fname=outname, speed_display=speed_display, redblue=redblue,
                     plotflag=plotflag, scale=scale, vlims=limits)


@click.command()
@click.option("--point", "-p", type=str, required=True, multiple=True, help="lon, lat of point(s) you want to extract.")
@click.option("--time_start", "-t0", type=str, default=None, help="Start Time (yyyy-mm-dd HH:MM:SS)")
@click.option("--time_end", "-t1", type=str, default=None, help="End Time (yyyy-mm-dd HH:MM:SS)")
@click.option("--path_data", "-pd", type=str, default=os.getcwd(), help="Path to radial files.")
@click.option("--path_save", "-ps", type=str, default=os.getcwd(), help="Path to save timeseries output")
@click.option("--fname", "-f", type=str, help="Export filename. The time range will be appended to this")
@click.option("--type",  "-o", type=str, default="csv", help="Format to save output: csv, netcdf")
@click.option('--parallel', '-m', is_flag=True, default=False, help="Enable parallel processing.")
def extract_timeseries(point, time_start, time_end, path_data, path_save, fname, type, parallel):
    click.echo("Executing hfradarpy.cli.extract_timeseries")
    click.echo("")
    os.makedirs(path_save, exist_ok=True)
    click.echo(f"Checking for radial data in {path_data}")

    files = sorted(glob.glob(os.path.join(path_data, '*.ruv')))
    click.echo(f"{len(files)} radial files found")
    df = list_to_dataframe(files)
    click.echo(f"Data from {df.iloc[0].name} to {df.iloc[-1].name}")

    if time_start or time_end:
        if time_start:
            time_start = pd.to_datetime(time_start).tz_localize(None)
        else:
            time_start = df.iloc[0].name

        if time_end:
            time_end = pd.to_datetime(time_end).tz_localize(None)
        else:
            time_end = df.iloc[-1].name

        click.echo(f"Subsetting from {time_start} to {time_end}")
        click.echo("")
        df = df[time_start:time_end]
        print(df.iloc[0])
        print(df.iloc[-1])

    try:
        ds = concat(df['file'].tolist(), method='gridded', enhance=True,
                    parallel=parallel)
    except ValueError as e:
        click.echo(f"{e}: no radial data found. Check your file path")
        return
       
    # Get the receiver location
    rxloc = [float(x) for x in ds.Origin.split('  ')]
    rxloc.reverse()

    click.echo("Extracting timeseries between the following points:")
    click.echo(f"Receiver - Lon: {rxloc[0]}, Lat: {rxloc[1]}")

    lons = []
    lats = []
    for p in point:
        pf = [float(x) for x in p.split(',')]
        click.echo(f"Terminus - Lon, {pf[0]} Lat: {pf[1]}")
        lons.append(pf[0])
        lats.append(pf[1])
    click.echo("")
    
    rxlon = np.full_like(lons, rxloc[0])
    rxlat = np.full_like(lats, rxloc[1])

    # Inverse transformation
    # Determine forward and back azimuths, plus distances between initial points and terminus points.
    forward_azimuth, _, distance = inv(rxlon, rxlat, lons, lats)
    click.echo("")
    click.echo(f"Finding nearest gridpoint(s) based on bearing and range.")
    click.echo(f'Bearing(s) (CWN): {forward_azimuth} degrees, Range(s): {distance} kilometers') # degrees, meters
    
    datasets = []
    for i,_ in enumerate(forward_azimuth):
        # Select nearest neighbor based on bearing and range
        datasets.append(ds.sel(bearing=forward_azimuth[i], range=distance[i], method='nearest'))
    merged = xr.concat(datasets, dim='time')
    click.echo(f'{len(np.unique(merged.time))} timesteps extracted from radial data.')

    tmin = pd.to_datetime(merged.time.min().data)
    tmax = pd.to_datetime(merged.time.max().data)

    if fname:
        fname = f'{fname}_extracted_timeseries_{tmin.strftime(fmt)}-{tmax.strftime(fmt)}'
    else:
        fname = f'{merged.Site.split()[0].lower()}_extracted_timeseries_{tmin.strftime(fmt)}-{tmax.strftime(fmt)}'

    if type == 'csv':
        save_file = os.path.join(path_save, fname + '.csv')
        merged.to_dataframe().to_csv(save_file)
    elif type == "netcdf":
        save_file = os.path.join(path_save, fname + '.nc')
        merged.to_netcdf(save_file) 
    click.echo("")
    click.echo(f'Timeseries successfully extracted from radial data and saved to {save_file}')

# Add each nested command to the main function here.
main.add_command(extract_timeseries)
main.add_command(plotruv)

if __name__ == "__main__":

    data_path = (Path(__file__).parent.with_name("examples") / "data").resolve()
    output_path = (Path(__file__).parent.with_name("examples") / "output").resolve()
    radial_dir = data_path / "radials" / "ruv"/ "SEAB"
    plotruv(
        [
            "--path_data", radial_dir,
            "--path_save", output_path,
            ]
        )

    data_path = (Path(__file__).parent.with_name("examples") / "data").resolve()
    output_path = (Path(__file__).parent.with_name("examples") / "output").resolve()
    radial_dir = data_path / "radials" / "ruv"/ "SEAB"
    point1 = "-73.63, 40.29"
    point2 = "-73.63, 40.31"
    time_start = "2019-01-01T00:00:00Z"
    time_end = "2019-01-02T00:00:00Z"
    extract_timeseries(
        [
            "--point", point1,
            "--point", point2,
            "--time_start", time_start,
            "--time_end", time_end,
            "--path_data", radial_dir,
            "--path_save", output_path,
            "--fname", 'test',
            "--type", 'csv',
            "--parallel"
            ]
        )



