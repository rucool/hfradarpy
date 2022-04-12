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

fmt = "%Y-%m-%dT%H%M%S"

@click.group()
def main(args=None):
    """Console script for HFRadarPy python package."""
    click.echo("Toolbox to read in High Frequency Radar (HFR) files written in CODAR Tabular Format (CTF).")
    click.echo(f"hfradarpy: {__version__}")

    return 0
    # pass

@click.command()
@click.option("--lon", type=float, help="Longitude of point(s)", multiple=True)
@click.option("--lat", type=float, help="Latitude of point(s)", multiple=True)
@click.option("--time_start", "-t0", default=None, type=str, help="Start Time (yyyy-mm-dd HH:MM:SS")
@click.option("--time_end", "-t1", default=None, type=str, help="End Time")
@click.option("--path_data", "-pd", default=os.getcwd(), type=str, help="Path containing radial files.")
@click.option("--path_save", "-ps", default=os.getcwd(), type=str,  help="Path to save output file")
@click.option("--fname", "-f", type=str, help="Export filename. The timeseries range is appended to this")
@click.option("--type",  "-o", default="csv", type=str, help="Format to save output: csv, netcdf")
def extract_timeseries(lon, lat, time_start, time_end, path_data, path_save, fname, type):
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
        ds = concat(df['file'].tolist(), method='gridded', enhance=True)
    except ValueError as e:
        click.echo(f"{e}: no radial data found. Check your file path")
        return
       
    # Get the receiver location
    rxloc = [float(x) for x in ds.Origin.split('  ')]
    rxloc.reverse()

    click.echo("Extracting timeseries between the following points:")
    click.echo(f"Receiver - Lon: {rxloc[0]}, Lat: {rxloc[1]}")

    for i,_ in enumerate(lon):
        click.echo(f"Terminus {i} - Lon, {lon[i]} Lat: {lat[i]}")
    click.echo("")
    
    rxlon = np.full_like(lon, rxloc[0])
    rxlat = np.full_like(lat, rxloc[1])

    # Inverse transformation
    # Determine forward and back azimuths, plus distances between initial points and terminus points.
    forward_azimuth, _, distance = inv(rxlon, rxlat, lon, lat)
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

if __name__ == "__main__":
    # -pd /Users/mikesmith/Documents/ugos/new/qc -o csv -ps /Users/mikesmith/Documents/ -f test -t0 2020-01-01 -t1 2021-01-03
    radial_dir = '/Users/mikesmith/Documents/ugos/new/qc/'
    lon = [-80.5, -81.5] 
    lat = [24.5, 24.25]
    time_start = '2019-01-01T00:00:00Z'
    time_end = '2020-01-02T00:00:00Z'
    time_start = '2020-01-01'
    time_end = '2020-01-03'
    path_save = '/Users/mikesmith/Documents/'
    extract_timeseries(
        [
            "--lon", lon[0],
            "--lat", lat[0],
            "--lon", lon[1],
            "--lat", lat[1],
            "--time_start", time_start,
            "--time_end", time_end,
            "--path_data", radial_dir,
            "--path_save", path_save,
            "--fname", 'test',
            "--type", 'csv'
            ]
        )


    # sys.exit(main())  # pragma: no cover
