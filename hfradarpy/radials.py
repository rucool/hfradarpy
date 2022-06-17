import datetime as dt
import numpy as np
import os
import pandas as pd
import re
import copy
import xarray as xr
from hfradarpy.ctf import CTFParser
from hfradarpy.calc import reckon
from hfradarpy.io.nc import make_encoding
from pathlib import Path
from joblib import Parallel, delayed
import multiprocessing
import logging

logger = logging.getLogger(__name__)


def qc_radial_file(radial_file, qc_values=None, export=None, save_path=None, clean=False, clean_path=None):
    """
    Main function to parse and qc radial files.

    Setting clean to True will create two separate quality controlled radial files. Must set clean_path.
    The first radial file with containing flag metadata will be saved to save_path. This file contains data along with flags.
    The second radial file with data that failed qc removed will be saved to clean_path. This file does not contain any flags.


    Args:
        radial_file (str or Path):
            Path to radial file or a Radial object
        qc_values (dict, optional):
            Dictionary containing thresholds for each QC test. Defaults to None.
        export (str, optional):
            None or 'radial' or 'netcdf-tabular' or 'netcdf-multidimensional'. Defaults to None.
        save_path (str or Path, optional):
            Path to save quality controlled radial file. Defaults to None.
        clean (bool, optional):
            Remove any row of data where the primary flag equals 4 (Failure flag). Defaults to False.
        clean_path (str or Path, optional):
            Path to save quality controlled radial file with data that fails qc removed. Defaults to None.

    Returns:
        Radial object: A quality controlled radial file.
    """
    qc_values = qc_values or dict(
        qc_qartod_avg_radial_bearing=dict(reference_bearing=151, warning_threshold=15, failure_threshold=30),
        qc_qartod_radial_count=dict(min_count=75.0, low_count=225.0),
        qc_qartod_maximum_velocity=dict(max_speed=300.0, high_speed=100.0),
        qc_qartod_spatial_median=dict(
            radial_smed_range_cell_limit=2.1, radial_smed_angular_limit=10, radial_smed_current_difference=30
        ),
        qc_qartod_temporal_gradient=dict(gradient_temp_fail=32, gradient_temp_warn=25),
        qc_qartod_primary_flag=dict(
            include=[
                "qc_qartod_syntax",
                "qc_qartod_valid_location",
                "qc_qartod_radial_count",
                "qc_qartod_maximum_velocity",
                "qc_qartod_spatial_median",
            ]
        ),
    )

    if not isinstance(radial_file, Radial):
        r = Radial(radial_file)
    else:
        r = radial_file

    if r.is_valid():
        if clean:
            rclean = copy.deepcopy(r)
        t0 = r.time - dt.timedelta(hours=1)
        previous_radial = "{}_{}{}".format(
            "_".join(r.file_name.split("_")[:2]), t0.strftime("%Y_%m_%d_%H00"), os.path.splitext(r.file_name)[1]
        )
        previous_full_file = os.path.join(os.path.dirname(r.full_file), previous_radial)
        qc_keys = qc_values.keys()

        # run high frequency radar qartod tests on open radial file
        r.initialize_qc()
        r.qc_qartod_syntax()

        if "qc_qartod_maximum_velocity" in qc_keys:
            r.qc_qartod_maximum_velocity(**qc_values["qc_qartod_maximum_velocity"])

        r.qc_qartod_valid_location()

        if "qc_qartod_radial_count" in qc_keys:
            r.qc_qartod_radial_count(**qc_values["qc_qartod_radial_count"])

        if "qc_qartod_spatial_median" in qc_keys:
            r.qc_qartod_spatial_median(**qc_values["qc_qartod_spatial_median"])

        if 'qc_qartod_temporal_gradient' in qc_keys:
            r.qc_qartod_temporal_gradient(previous_full_file,**qc_values['qc_qartod_temporal_gradient'])

        if "qc_qartod_avg_radial_bearing" in qc_keys:
            r.qc_qartod_avg_radial_bearing(**qc_values["qc_qartod_avg_radial_bearing"])

        if "qc_qartod_primary_flag" in qc_keys:
            r.qc_qartod_primary_flag(**qc_values["qc_qartod_primary_flag"])

        if clean:
            d = rclean.data
            dqc = r.data
            if "PRIM" in r.data:
                rt = d[dqc["PRIM"] != 4]
                rclean.data = rt

                for key in rclean._tables.keys():
                    table = rclean._tables[key]
                    if "LLUV" in table["TableType"]:
                        rclean._tables[1]["TableRows"] = rt.shape[0]
                # else:
                #   warning that it didn't update number of table rows
            # else:
            # warning of failure to update file, the original will be exported

        # Export radial file to either a radial or netcdf
        if export:
            if export == "radial":
                r.to_ruv(os.path.join(save_path, r.file_name))
            elif export == "netcdf-tabular":
                r.to_netcdf(os.path.join(save_path, r.file_name), "tabular")
            elif export == "netcdf-gridded":
                r.to_netcdf(os.path.join(save_path, r.file_name), "gridded")

            if clean:
                if export == "radial":
                    rclean.to_ruv(os.path.join(clean_path, rclean.file_name))
                elif export == "netcdf-tabular":
                    rclean.to_netcdf(os.path.join(save_path, rclean.file_name), "tabular", prepend_extension=True)
                elif export == "netcdf-gridded":
                    rclean.to_netcdf(os.path.join(save_path, rclean.file_name), "gridded", prepend_extension=True)
        else:
            return r


def concat(rlist, method="gridded", enhance=False, parallel=False):
    """
    This function takes a list of Radial objects or radial file paths and
    combines them along the time dimension using xarrays built-in concatenation
    routines.

    Args:
        rlist (list):
            list of radial files or Radial objects that you want to concatenate
        method (str, optional):
            'gridded' or 'tabular'. Defaults to 'gridded'
        enhance (bool, optional):
            Change manufacturer variables to cf compliant variable names. Add CF attributes and metadata. Defaults to False.
        parallel (bool, optional): Enables parallel processing. Defaults to False.

    Returns:
        xarray dataset: radials concatenated into an xarray dataset
    """
    def load_radials(radial, method, enhance=False):
        if not isinstance(radial, Radial):
            radial = Radial(radial)
            ds = radial.to_xarray(method, enhance=enhance)
        return (radial.file_name, ds)

    if parallel:
        num_cores = multiprocessing.cpu_count()
        radials = Parallel(n_jobs=num_cores)(
            delayed(load_radials)(radial=r, method=method, enhance=enhance) for r in rlist
        )
        radial_dict = {radial: ds for (radial, ds) in radials}
    else:
        radial_dict = {}
        for radial in rlist:
            if not isinstance(radial, Radial):
                radial = Radial(radial)
            radial_dict[radial.file_name] = radial.to_xarray(method, enhance=enhance)

    ds = xr.concat(radial_dict.values(), "time")
    return ds.sortby("time")


class Radial(CTFParser):
    """
    Radial Subclass.

    This class should be used when loading a CODAR radial (.ruv) file. This class utilizes the generic LLUV class from
    ~/hfradarpy/ctf.py in order to load CODAR Radial files
    """

    def __init__(self, fname, replace_invalid=True, empty_radial=False, vflip=False):
        """
        Initialize the radial object

        Args:
            fname (str or Path): Full file path to the radial .ruv to be loaded
            replace_invalid (bool, optional): Replace invalid/dummy manufacturer fill values with NaN. Defaults to True.
            empty_radial (bool, optional): Returns an empty Radial object. Defaults to False.
            vflip (bool, optional): Flip sign of velocities [MINV, MAXV, VELO] so + is away and - is towards the radar. Defaults to False.
        """
        logging.info("Loading radial file: {}".format(fname))
        super().__init__(fname)

        if self._iscorrupt:
            return

        self.data = pd.DataFrame()
        self.velocity_sign = "+: towards"

        for key in self._tables.keys():
            table = self._tables[key]
            if "LLUV" in table["TableType"]:
                self.data = table["data"]
            elif "rads" in table["TableType"]:
                self.diagnostics_radial = table["data"]
                self.diagnostics_radial["datetime"] = self.diagnostics_radial[
                    ["TYRS", "TMON", "TDAY", "THRS", "TMIN", "TSEC"]
                ].apply(lambda s: dt.datetime(*s), axis=1)
            elif "rcvr" in table["TableType"]:
                self.diagnostics_hardware = table["data"]
                self.diagnostics_hardware["datetime"] = self.diagnostics_hardware[
                    ["TYRS", "TMON", "TDAY", "THRS", "TMIN", "TSEC"]
                ].apply(lambda s: dt.datetime(*s), axis=1)
            elif "RINF" in table["TableType"]:
                self.range_information = table["data"]
            elif "MRGS" in table["TableType"]:
                self.merge_information = table["data"]

        if not self.data.empty:
            if replace_invalid:
                self.replace_invalid_values()

            if empty_radial:
                self.empty_radial()

        if vflip:
            # CODAR standard: + velocities towards radar, - velocities away from radar
            # Science standard: - velocities towards radar, + velocities away from radar
            # Flip sign so positive velocities are away from the radar
            self.flip_velocities()

    def __repr__(self):
        """
        Represent class's object as a string
        """
        return "<Radial: {}>".format(self.file_name)

    def flip_velocities(self):
        # CODAR standard: + velocities towards radar, - velocities away from radar
        # Science standard: - velocities towards radar, + velocities away from radar
        # Flip sign so positive velocities are away from the radar as per cf conventions
        logging.info(f"Flipping the sign of radial velocities")
        flips = ["MINV", "MAXV", "VELO"]
        for f in flips:
            if f in self.data:
                self.data[f] = -self.data[f]
        self.flipped_velocites = True
        
        if self.velocity_sign == "+: towards":
            self.velocity_sign = "-: away"
        elif self.velocity_sign == "-: away":
            self.velocity_sign = "+: towards"

    def empty_radial(self):
        """
        Defines an empty radial
        """

        self.file_path = ""
        self.file_name = ""
        self.full_file = ""
        self.metadata = ""
        self.is_wera = False
        self._iscorrupt = False
        self.time = []

        for key in self._tables.keys():
            table = self._tables[key]
            self._tables[key]["TableRows"] = "0"
            if "LLUV" in table["TableType"]:
                self.data.drop(self.data.index[:], inplace=True)
                self._tables[key]["data"] = self.data
            elif "rads" in table["TableType"]:
                self.diagnostics_radial.drop(self.diagnostics_radial.index[:], inplace=True)
                self._tables[key]["data"] = self.diagnostics_radial
            elif "rcvr" in table["TableType"]:
                self.diagnostics_hardware.drop(self.diagnostics_hardware.index[:], inplace=True)
                self._tables[key]["data"] = self.diagnostics_hardware
            elif "RINF" in table["TableType"]:
                self.range_information.drop(self.range_information.index[:], inplace=True)
                self._tables[key]["data"] = self.range_information

    def mask_over_land(self, subset=False):
        """
        Return mask of radial data that lies over the land.

        If subset is set to True, you will get a pd.Dataframe back with data over land removed

        Args:
            subset (bool, optional): Removes any points not over water from the dataset. Defaults to True.

        Returns:
            self.data (pd.DataFrame): subset = True, a pandas dataframe is returned with data that is over water.
            or
            water_index (pd.Series): subset = False, pandas series containing the index of points of water is returned
        """
        try:
            import geopandas as gpd
            from shapely.geometry import Point
        except ImportError:
            raise NotImplementedError("You need geopandas>=0.10.0 and shapely installed in order to apply land mask")
        else:
            logging.info("Masking radials over land")

            land = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
            land = land[land["continent"] == "North America"]

            geodata = gpd.GeoDataFrame(
                self.data[["LOND", "LATD"]],
                crs="EPSG:4326",
                geometry=[Point(xy) for xy in zip(self.data.LOND.values, self.data.LATD.values)],
            )
            # Join the geodataframe containing radial points with geodataframe containing leasing areas
            geodata = gpd.tools.sjoin(geodata.to_crs(4326), land.to_crs(4326), how="left", predicate="intersects")

            # All data in the continent column that lies over water should be nan.
            water_index = geodata["continent"].isna()

            if subset:
                # Subset the data to water only
                self.data = self.data.loc[water_index].reset_index(drop=True)
            else:
                return water_index

    def to_xarray(self, model="tabular", enhance=False, range_minmax=None, bearing=None):
        """
        Helper function

        Args:
            model (str, optional):
                Create a 'tabular' (time) or 'gridded' (time, range, bearing) xarray dataset. Defaults to 'tabular'.
            enhance (bool, optional):
                Change manufacturer variables to meaningful variable names. Add attributes and other metadata. Defaults to True
            range_minmax (list or tuple, optional):
                For gridded dataset only, Range minimum (km) and range maximum (km) of the radial.
                Defaults to None which automatically infers the range based off of the radial data.
            bearing (list or tuple, optional):
                For gridded dataset only, Antenna bearing (degrees) and angular resolution (degrees) of the radial grid.
                Defaults to None which automatically infers the bearing based off of the radial data.
        """

        if model == "tabular":
            ds = self._to_xarray_tabular(enhance)
        elif model == "gridded":
            ds = self._to_xarray_gridded(range_minmax=range_minmax, bearing=bearing, enhance=enhance)
        else:
            raise ValueError("Please enter a valid data model type. Must be a string 'tabular' or 'gridded'")
        return ds

    def _to_xarray_gridded(self, range_minmax=None, bearing=None, enhance=True):
        """
        Convert radial file to an xarray dataset on a radial grid.
        This dataset has dimensions of time,range, and bearing.
        Add list range_min and range_max if you want to make the file only include range cells between the two.

        Adapted from MATLAB code by Mark Otero
        http://cordc.ucsd.edu/projects/mapping/documents/HFRNet_Radial_NetCDF.pdf

        Args:
            range_minmax (list or tuple, optional): Range minimum (km) and range maximum (km) of the radial.
                Defaults to None which automatically infers the range based off of the radial data.
            bearing (list or tuple, optional): Antenna bearing (degrees) and angular resolution (degrees) of the radial grid.
                Defaults to None which automatically infers the bearing based off of the radial data.
            enhance (bool, optional):
                Change manufacturer variables to cf compliant variable names. Add CF attributes and metadata. Defaults to True.

        Returns:
            xarray.Dataset: radial file converted to an xarray dataset with dimensions of time, range, and bearing.
        """
        logging.info("Converting radial matrix to multidimensional dataset")

        # CF Standard: T, Z, Y, X
        coords = ("time", "range", "bearing")

        # Intitialize empty xarray dataset
        ds = xr.Dataset()

        if range_minmax:
            # range_include = True, because it is not a NoneType
            if isinstance(range_minmax, (list, tuple)):
                # range_include must be a tuple or a list that has data inside
                range_min = np.float32(range_minmax[0])
                range_max = np.float32(range_minmax[1])
            else:
                # range_minmax is not a tuple or a list, so error out and return to the original function
                raise TypeError(
                    "range_minmax must be a list/tuple or None. See function documentation"
                )
        else:
            # range_minmax = False, because range_minmax is a NoneType.
            logging.info("Inferring min/max ranges based off of radial data")
            range_min = self.data.RNGE.min()
            range_max = self.data.RNGE.max()

        # Get the range resolution of the radar site from the metadata
        range_step = np.double(self.metadata["RangeResolutionKMeters"].split()[0])

        # Create an array from range_min to range_max with range_step between each range
        range_dim = np.arange(range_min, np.round(range_max + range_step, 4), range_step, dtype=np.float32)

        # round floats to 4 decimal places due to inaccuracy of np.arange
        range_dim = [round(x, 4) for x in range_dim]

        # If bearing has user input and
        if bearing:
            # bearing = True, because it is not a NoneType
            if isinstance(bearing, (list, tuple)):
                # bearing must be a tuple or a list contains data
                antenna_bearing = np.float32(bearing[0])
                angular_resolution = np.float32(bearing[1])
            else:
                # bearing is not a tuple or a list, so error out and return to the original function
                raise TypeError(
                    "bearing must be a list/tuple or None. See function documentation"
                )
        else:
            # bearing is a NoneType
            try:
                # Infer it from the radial metadata
                angular_resolution = np.double(self.metadata["AngularResolution"].split()[0])
                antenna_bearing = np.double(self.metadata["AntennaBearing"].split()[0])
            except KeyError:
                logging.info(
                    "AngularResolution and/or AntennaBearing not find in radial file metadata. "
                    "Creating a generic 360 degree grid with 1 degree bin resolution."
                )
                # If AngularResolution and/or AntennaBearing aren't in the radial file,
                # let's create a generic 360 degree grid with 1 degree resolution.
                angular_resolution = 1
                antenna_bearing = 0

        # Calculate bearings with the AntennaBearing as the min
        # and AntennaBearing+360 as the max
        # with angular_resolution as the spacing between each bearing.
        #  We use np.mod(bearings, 360) to convert back to 0 to 360
        bearing_dim = sorted(
            np.mod(np.arange(antenna_bearing, antenna_bearing + 360, angular_resolution, dtype=np.short), 360)
        )

        # create radial grid from bearing and range
        [bearings, ranges] = np.meshgrid(bearing_dim, range_dim)

        # calculate lat/lons from origin, bearing, and ranges
        latlon = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", self.metadata["Origin"])]
        lons = np.full_like(bearings, latlon[1], dtype=float)
        lats = np.full_like(bearings, latlon[0], dtype=float)
        lond, latd = reckon(lons, lats, bearings, ranges)

        # create dictionary containing variables from dataframe in the shape of radial grid
        d = {key: np.tile(np.nan, bearings.shape) for key in self.data.keys()}

        # find grid indices from radial grid (bearing, ranges)
        range_map_idx = np.tile(np.nan, self.data["RNGE"].shape)
        bearing_map_idx = np.tile(np.nan, self.data["BEAR"].shape)

        for i, line in enumerate(self.data["RNGE"]):
            range_map_idx[i] = np.argmin(np.abs(range_dim - self.data.RNGE[i]))
            bearing_map_idx[i] = np.argmin(np.abs(bearing_dim - self.data.BEAR[i]))

        for k, v in d.items():
            v[range_map_idx.astype(int), bearing_map_idx.astype(int)] = self.data[k]
            d[k] = v

        # Add extra dimension for time
        d = {k: np.expand_dims(np.float32(v), axis=0) for (k, v) in d.items()}

        # Add coordinate variables to dataset
        timestamp = dt.datetime(*[int(s) for s in self.metadata["TimeStamp"].split()])
        ds.coords["bearing"] = bearing_dim
        ds.coords["range"] = range_dim
        ds.coords["time"] = pd.date_range(timestamp, periods=1)
        ds.coords["lon"] = (("range", "bearing"), lond.round(4))
        ds.coords["lat"] = (("range", "bearing"), latd.round(4))

        # Add all variables to dataset
        for k, v in d.items():
            ds[k] = (coords, v)

        # Check if calculated longitudes and latitudes align with given longitudes and latitudes
        # plt.plot(ds.lon, ds.lat, 'bo', ds.LOND.squeeze(), ds.LATD.squeeze(), 'rx')

        # Drop extraneous variables
        ds = ds.drop_vars(["LOND", "LATD", "BEAR", "RNGE"])

        # Assign header data to global attributes
        ds = ds.assign_attrs(self.metadata)

        if enhance is True:
            ds = self.enhance_xarray(ds)
            ds = xr.decode_cf(ds)
        return ds

    def _to_xarray_tabular(self, enhance=False):
        """
        Convert radial file to an xarray dataset on a radial grid.
        This dataset has one dimension, time.

        Args:
            enhance (bool, optional): Rename variables to something meaningful and add cf attributes. Defaults to False.

        Returns:
            xarray.Dataset: Enhanced dataset with renamed variables and useful attributes.
        """
        logging.info("Converting radial matrix to tabular dataset")

        # Clean radial header
        # self.clean_header()

        # get timestamp from radial metadata
        timestamp = dt.datetime(*[int(s) for s in self.metadata["TimeStamp"].split()])

        self.data["time"] = timestamp
        self.data.set_index("time", inplace=True)

        # Intitialize xarray dataset
        ds = self.data.to_xarray()
        # ds = self.data.to_xarray().expand_dims('time').assign_coords(time=('time', [timestamp]))
        # ds.coords['time'] = pd.date_range(timestamp, periods=1)
        # ds.expand_dims('time').assign_coords(time=('time', [timestamp]))

        # Check if calculated longitudes and latitudes align with given longitudes and latitudes
        # plt.plot(ds.lon, ds.lat, 'bo', ds.LOND.squeeze(), ds.LATD.squeeze(), 'rx')

        # Assign header data to global attributes
        ds = ds.assign_attrs(self.metadata)

        if enhance is True:
            ds = self.enhance_xarray(ds)
            ds = xr.decode_cf(ds)

        return ds

    @staticmethod
    def enhance_xarray(xds):
        """
        Change manufacturer variables to meaningful variable names.
        Add attributes and other metadata.

        Args:
            xds (xarray.Dataset): radial xarray Dataset to enhance

        Returns:
            xds (xarray.Dataset): enhanced xarray dataset
        """
        rename = dict(
            VELU="u",
            VELV="v",
            VFLG="vector_flag",
            ESPC="spatial_quality",
            ETMP="temporal_quality",
            MAXV="velocity_max",
            MINV="velocity_min",
            ERSC="spatial_count",
            ERTC="temporal_count",
            XDST="dist_east_from_origin",
            YDST="dist_north_from_origin",
            VELO="velocity",
            HEAD="heading",
            SPRC="range_cell",
            EACC="accuracy",  # WERA specific
            LOND="lon",
            LATD="lat",
            BEAR="bearing",
            RNGE="range",
        )

        rename_qc = dict()

        # rename variables to something meaningful if they existin
        # in the xarray dataset
        existing_renames = {k: v for k, v in rename.items() if k in xds}
        xds = xds.rename(existing_renames)

        # set time attribute
        xds["time"].attrs["standard_name"] = "time"

        # Set lon attributes
        xds["lon"].attrs["long_name"] = "Longitude"
        xds["lon"].attrs["standard_name"] = "longitude"
        xds["lon"].attrs["short_name"] = "lon"
        xds["lon"].attrs["units"] = "degrees_east"
        xds["lon"].attrs["axis"] = "X"
        xds["lon"].attrs["valid_min"] = np.float32(-180.0)
        xds["lon"].attrs["valid_max"] = np.float32(180.0)

        # Set lat attributes
        xds["lat"].attrs["long_name"] = "Latitude"
        xds["lat"].attrs["standard_name"] = "latitude"
        xds["lat"].attrs["short_name"] = "lat"
        xds["lat"].attrs["units"] = "degrees_north"
        xds["lat"].attrs["axis"] = "Y"
        xds["lat"].attrs["valid_min"] = np.float32(-90.0)
        xds["lat"].attrs["valid_max"] = np.float32(90.0)

        # Set u attributes
        xds["u"].attrs["long_name"] = "Eastward Surface Current (cm/s)"
        xds["u"].attrs["standard_name"] = "surface_eastward_sea_water_velocity"
        xds["u"].attrs["short_name"] = "u"
        xds["u"].attrs["units"] = "cm s-1"
        xds["u"].attrs["valid_min"] = np.float32(-300)
        xds["u"].attrs["valid_max"] = np.float32(300)
        xds["u"].attrs["coordinates"] = "lon lat"
        xds["u"].attrs["grid_mapping"] = "crs"

        # Set v attributes
        xds["v"].attrs["long_name"] = "Northward Surface Current (cm/s)"
        xds["v"].attrs["standard_name"] = "surface_northward_sea_water_velocity"
        xds["v"].attrs["short_name"] = "v"
        xds["v"].attrs["units"] = "cm s-1"
        xds["v"].attrs["valid_min"] = np.float32(-300)
        xds["v"].attrs["valid_max"] = np.float32(300)
        xds["v"].attrs["coordinates"] = "lon lat"
        xds["v"].attrs["grid_mapping"] = "crs"

        # Set bearing attributes
        xds["bearing"].attrs["long_name"] = "Bearing from origin (away from instrument)"
        xds["bearing"].attrs["short_name"] = "bearing"
        xds["bearing"].attrs["units"] = "degrees"
        xds["bearing"].attrs["valid_min"] = np.float32(0)
        xds["bearing"].attrs["valid_max"] = np.float32(360)
        xds["bearing"].attrs["grid_mapping"] = "crs"
        xds["bearing"].attrs["axis"] = "Y"

        # Set range attributes
        xds["range"].attrs["long_name"] = "Range from origin (away from instrument)"
        xds["range"].attrs["short_name"] = "range"
        xds["range"].attrs["units"] = "km"
        xds["range"].attrs["valid_min"] = np.float32(0)
        xds["range"].attrs["valid_max"] = np.float32(1000)
        xds["range"].attrs["grid_mapping"] = "crs"
        xds["range"].attrs["axis"] = "X"

        # velocity
        xds["velocity"].attrs["valid_range"] = [-1000, 1000]
        xds["velocity"].attrs["standard_name"] = "radial_sea_water_velocity_away_from_instrument"
        xds["velocity"].attrs["units"] = "cm s-1"
        xds["velocity"].attrs["coordinates"] = "lon lat"
        xds["velocity"].attrs["grid_mapping"] = "crs"

        # heading
        if "heading" in xds:
            xds["heading"].attrs["valid_range"] = [0, 3600]
            xds["heading"].attrs["standard_name"] = "direction_of_radial_vector_away_from_instrument"
            xds["heading"].attrs["units"] = "degrees"
            xds["heading"].attrs["coordinates"] = "lon lat"
            xds["heading"].attrs["scale_factor"] = 0.1
            xds["heading"].attrs["grid_mapping"] = "crs"

        # vector_flag
        if "vector_flag" in xds:
            xds["vector_flag"].attrs["long_name"] = "Vector Flag Masks"
            xds["vector_flag"].attrs["valid_range"] = [0, 2048]
            xds["vector_flag"].attrs["flag_masks"] = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
            xds["vector_flag"].attrs[
                "flag_meanings"
            ] = "grid_point_deleted grid_point_near_coast point_measurement no_radial_solution baseline_interpolation"\
                "exceeds_max_speed invalid_solution solution_beyond_valid_spatial_domain insufficient_angular_resolution"\
                "reserved reserved"
            xds["vector_flag"].attrs["coordinates"] = "lon lat"
            xds["vector_flag"].attrs["grid_mapping"] = "crs"

        # spatial_quality
        if "spatial_quality" in xds:
            xds["spatial_quality"].attrs["long_name"] = "Spatial Quality of radial sea water velocity"
            xds["spatial_quality"].attrs["units"] = "cm s-1"
            xds["spatial_quality"].attrs["coordinates"] = "lon lat"
            xds["spatial_quality"].attrs["grid_mapping"] = "crs"

        # temporal_quality
        if "temporal_quality" in xds:
            xds["temporal_quality"].attrs["long_name"] = "Temporal Quality of radial sea water velocity"
            xds["temporal_quality"].attrs["units"] = "cm s-1"
            xds["temporal_quality"].attrs["coordinates"] = "lon lat"
            xds["temporal_quality"].attrs["grid_mapping"] = "crs"

        # velocity_max
        if "velocity_max" in xds:
            xds["velocity_max"].attrs["long_name"] = "Maximum Velocity of sea water (away from instrument)"
            xds["velocity_max"].attrs["units"] = "cm s-1"
            xds["velocity_max"].attrs["coordinates"] = "lon lat"
            xds["velocity_max"].attrs["grid_mapping"] = "crs"

        # velocity_min
        if "velocity_min" in xds:
            xds["velocity_min"].attrs["long_name"] = "Minimum Velocity of sea water (away from instrument)"
            xds["velocity_min"].attrs["units"] = "cm s-1"
            xds["velocity_min"].attrs["coordinates"] = "lon lat"
            xds["velocity_min"].attrs["grid_mapping"] = "crs"

        # spatial_count
        if "spatial_count" in xds:
            xds["spatial_count"].attrs["long_name"] = "Spatial count of sea water velocity (away from instrument)"
            xds["spatial_count"].attrs["coordinates"] = "lon lat"
            xds["spatial_count"].attrs["grid_mapping"] = "crs"

        # temporal_count
        if "temporal_count" in xds:
            xds["temporal_count"].attrs["long_name"] = "Temporal count of sea water velocity (away from instrument)"
            xds["temporal_count"].attrs["coordinates"] = "lon lat"
            xds["temporal_count"].attrs["grid_mapping"] = "crs"

        # east_dist_from_origin
        if "dist_east_from_origin" in xds:
            xds["dist_east_from_origin"].attrs["long_name"] = "Eastward distance from instrument"
            xds["dist_east_from_origin"].attrs["units"] = "km"
            xds["dist_east_from_origin"].attrs["coordinates"] = "lon lat"
            xds["dist_east_from_origin"].attrs["grid_mapping"] = "crs"

        # north_dist_from_origin
        if "dist_north_from_origin" in xds:
            xds["dist_north_from_origin"].attrs["long_name"] = "Northward distance from instrument"
            xds["dist_north_from_origin"].attrs["units"] = "km"
            xds["dist_north_from_origin"].attrs["coordinates"] = "lon lat"
            xds["dist_north_from_origin"].attrs["grid_mapping"] = "crs"

        # range_cell
        if "range_cell" in xds:
            xds["range_cell"].attrs["long_name"] = "Cross Spectra Range Cell  of sea water velocity (away from instrument)"
            xds["range_cell"].attrs["coordinates"] = "lon lat"
            xds["range_cell"].attrs["grid_mapping"] = "crs"

        # range_cell
        if "accuracy" in xds:
            xds["accuracy"].attrs["long_name"] = "Accuracy"
            xds["accuracy"].attrs["coordinates"] = "lon lat"
            xds["accuracy"].attrs["grid_mapping"] = "crs"
            xds["accuracy"].attrs["units"] = "cm s-1"

        # QC06
        if "QC06" in xds:
            xds["QC06"].attrs["long_name"] = "Syntax (QARTOD Test 06) Flag Masks"
            xds["QC06"].attrs["valid_range"] = [1, 9]
            xds["QC06"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["QC06"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["QC06"].attrs["coordinates"] = "lon lat"
            xds["QC06"].attrs["grid_mapping"] = "crs"
            rename_qc["QC06"] = "syntax_qc"

        # QC07
        if "QC07" in xds:
            xds["QC07"].attrs["long_name"] = "Maximum Velocity Threshold (QARTOD Test 07) Flag Masks"
            xds["QC07"].attrs["valid_range"] = [1, 9]
            xds["QC07"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["QC07"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["QC07"].attrs["coordinates"] = "lon lat"
            xds["QC07"].attrs["grid_mapping"] = "crs"
            rename_qc["QC07"] = "max_threshold_qc"

        # QC08
        if "QC08" in xds:
            xds["QC08"].attrs["long_name"] = "Valid Location (QARTOD Test 08) Flag Masks"
            xds["QC08"].attrs["valid_range"] = [1, 9]
            xds["QC08"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["QC08"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["QC08"].attrs["coordinates"] = "lon lat"
            xds["QC08"].attrs["grid_mapping"] = "crs"
            rename_qc["QC08"] = "valid_location_qc"

        # QC09
        if "QC09" in xds:
            xds["QC09"].attrs["long_name"] = "Radial Count (QARTOD Test 09) Flag Masks"
            xds["QC09"].attrs["valid_range"] = [1, 9]
            xds["QC09"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["QC09"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["QC09"].attrs["coordinates"] = "lon lat"
            xds["QC09"].attrs["grid_mapping"] = "crs"
            rename_qc["QC09"] = "radial_count_qc"

        # QC10
        if "QC10" in xds:
            xds["QC10"].attrs["long_name"] = "Spatial Median Filter (QARTOD Test 10) Flag Masks"
            xds["QC10"].attrs["valid_range"] = [1, 9]
            xds["QC10"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["QC10"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["QC10"].attrs["coordinates"] = "lon lat"
            xds["QC10"].attrs["grid_mapping"] = "crs"
            rename_qc["QC10"] = "spatial_median_filter_qc"

        # QC11
        if "QC11" in xds:
            xds["QC11"].attrs["long_name"] = "Temporal Gradient (QARTOD Test 11) Flag Masks"
            xds["QC11"].attrs["valid_range"] = [1, 9]
            xds["QC11"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["QC11"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["QC11"].attrs["coordinates"] = "lon lat"
            xds["QC11"].attrs["grid_mapping"] = "crs"
            rename_qc["QC11"] = "temporal_gradient_qc"

        # QC12
        if "QC12" in xds:
            xds["QC12"].attrs["long_name"] = "Average Radial Bearing (QARTOD Test 12) Flag Masks"
            xds["QC12"].attrs["valid_range"] = [1, 9]
            xds["QC12"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["QC12"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["QC12"].attrs["coordinates"] = "lon lat"
            xds["QC12"].attrs["grid_mapping"] = "crs"
            rename_qc["QC12"] = "average_radial_bearing_qc"

        # QC12
        if "PRIM" in xds:
            xds["PRIM"].attrs["long_name"] = "Primary Flag Masks"
            xds["PRIM"].attrs["valid_range"] = [1, 9]
            xds["PRIM"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["PRIM"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["PRIM"].attrs["coordinates"] = "lon lat"
            xds["PRIM"].attrs["grid_mapping"] = "crs"
            rename_qc["PRIM"] = "primary_flag_qc"

        if "QCOP" in xds:
            xds["QCOP"].attrs["long_name"] = "Operator Flag Masks"
            xds["QCOP"].attrs["valid_range"] = [1, 9]
            xds["QCOP"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["QCOP"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["QCOP"].attrs["coordinates"] = "lon lat"
            xds["QCOP"].attrs["grid_mapping"] = "crs"
            xds["QCOP"].attrs[
                "comment"
            ] = "Flag that is manually set by operator. Flag vectors that are not detected by QC tests but are clearly wrong."
            rename_qc["QCOP"] = "operator_flag_qc"

        # rename variables to something meaningful if they exist in the xarray dataset
        xds = xds.rename({k: v for k, v in rename_qc.items() if k in xds})
        # del xds.attrs['TimeStamp']

        return xds

    def file_type(self):
        """
        Return a string representing the type of file this is.
        """
        return "radial"

    def clean_header(self, split_origin=False):
        """
        Clean up the radial header dictionary so that you can upload it to the HFR MySQL Database.

        Args:
            split_origin (bool, optional):
                Split the 'Origin' string for the receiver location into a list of floats [lon, lat]. Defaults to False.
        """
        keep = [
            "CTF",
            "FileType",
            "LLUVSpec",
            "UUID",
            "Manufacturer",
            "Site",
            "TimeStamp",
            "TimeZone",
            "TimeCoverage",
            "Origin",
            "GreatCircle",
            "GeodVersion",
            "LLUVTrustData",
            "RangeStart",
            "RangeEnd",
            "RangeResolutionKMeters",
            "AntennaBearing",
            "ReferenceBearing",
            "AngularResolution",
            "SpatialResolution",
            "PatternType",
            "PatternDate",
            "PatternResolution",
            "TransmitCenterFreqMHz",
            "DopplerResolutionHzPerBin",
            "FirstOrderMethod",
            "BraggSmoothingPoints",
            "CurrentVelocityLimits",
            "BraggHasSecondOrder",
            "RadialBraggPeakDropOff",
            "RadialBraggPeakNull",
            "RadialBraggNoiseThreshold",
            "PatternAmplitudeCorrections",
            "PatternPhaseCorrections",
            "PatternAmplitudeCalculations",
            "PatternPhaseCalculations",
            "RadialMusicParameters",
            "MergedCount",
            "RadialMinimumMergePoints",
            "FirstOrderCalc",
            "MergeMethod",
            "PatternMethod",
            "TransmitSweepRateHz",
            "TransmitBandwidthKHz",
            "SpectraRangeCells",
            "SpectraDopplerCells",
            "TableType",
            "TableColumns",
            "TableColumnTypes",
            "TableRows",
            "TableStart",
            "CurrentVelocityLimit",
        ]

        # TableColumnTypes
        key_list = list(self.metadata.keys())
        for key in key_list:
            if key not in keep:
                del self.metadata[key]

        for k, v in self.metadata.items():
            if "Site" in k:
                # WERA has lines like: '%Site: csw "CSW' and '%Site: gtn "gtn'
                # This should work for both CODAR and WERA files
                split_site = v.split(" ", 1)[0]
                self.metadata[k] = "".join(e for e in split_site if e.isalnum())
            elif k in ("TimeStamp", "PatternDate"):
                try:
                    t_list = [int(s) for s in v.split()]
                    self.metadata[k] = dt.datetime(*t_list)
                except ValueError:
                    # Can't parse a date, set to None
                    self.metadata[k] = None
            elif "TimeZone" in k:
                self.metadata[k] = v.split('"')[1]
            elif "TableColumnTypes" in k:
                self.metadata[k] = " ".join([x.strip() for x in v.strip().split(" ")])
            elif "Origin" in k:
                if split_origin:
                    self.metadata[k] = re.findall(r"[-+]?\d*\.\d+|\d+", v)
                else:
                    self.metadata[k] = v.lstrip()
            elif k in (
                "RangeStart",
                "RangeEnd",
                "AntennaBearing",
                "ReferenceBearing",
                "AngularResolution",
                "SpatialResolution",
                "FirstOrderMethod",
                "BraggSmoothingPoints",
                "BraggHasSecondOrder",
                "MergedCount",
                "RadialMinimumMergePoints",
                "FirstOrderCalc",
                "SpectraRangeCells",
                "SpectraDopplerCells",
                "TableColumns",
                "TableRows",
                "PatternResolution",
                "CurrentVelocityLimit",
                "TimeCoverage",
            ):
                try:
                    self.metadata[k] = int(v)
                except ValueError:
                    temp = v.split(" ")[0]
                    try:
                        self.metadata[k] = int(temp)
                    except ValueError:
                        try:
                            self.metadata[k] = int(temp.split(".")[0])
                        except ValueError:
                            self.metadata[k] = None
            elif k in (
                "RangeResolutionKMeters",
                "CTF",
                "TransmitCenterFreqMHz",
                "DopplerResolutionHzPerBin",
                "RadialBraggPeakDropOff",
                "RadialBraggPeakNull",
                "RadialBraggNoiseThreshold",
                "TransmitSweepRateHz",
                "TransmitBandwidthKHz",
            ):
                try:
                    self.metadata[k] = float(v)
                except ValueError:
                    try:
                        self.metadata[k] = float(v.split(" ")[0])
                    except ValueError:
                        self.metadata[k] = None
            else:
                continue

        required = ["Origin", "TransmitCenterFreqMHz"]
        present_keys = self.metadata.keys()
        for key in required:
            if key not in present_keys:
                self.metadata[key] = None

    def to_netcdf(self, filename, model="tabular", prepend_extension=False, range_minmax=None, bearing=None, enhance=True,
                  user_attributes=None):
        """
        Create a compressed netCDF4 (.nc) file from the Radial object

        Args:
            filename (str or Path): User defined filename of radial file you want to save
            model (str, optional): Create a 'tabular' (time) or 'gridded' (time, range, bearing) NetCDF. Defaults to 'tabular'
            prepend_extension (bool, optional):
                Prepend a descriptive term (bearing or gridded) to the .nc extension. Defaults to False.
            range_minmax (list or tuple, optional):
                For gridded NetCDF only, Range minimum (km) and range maximum (km) of the radial.
                Defaults to None which automatically infers the range based off of the radial data.
            bearing (list or tuple, optional):
                For gridded NetCDF only, Antenna bearing (degrees) and angular resolution (degrees) of the radial grid.
                Defaults to None which automatically infers the bearing based off of the radial data.
            enhance (bool, optional):
                Change manufacturer variables to meaningful variable names. Add attributes and other metadata. Defaults to True
        """
        # Make sure filename is converted into a Path object
        filename = Path(filename)
        os.makedirs(filename.parent.resolve(), exist_ok=True)

        if not self.is_valid():
            raise ValueError("Could not export ASCII data, the input file was invalid.")

        # If the filename does not have a .nc extension, we will add one.
        if ".nc" not in str(filename):
            filename = filename.with_suffix(".nc")

        # If the outputted file exists already, delete the existing file
        if os.path.isfile(filename):
            os.remove(filename)

        xds = self.to_xarray(model, range_minmax=range_minmax, bearing=bearing, enhance=enhance)

        # if 'tabular' in netcdf_type:
        #     xds = self.to_xarray_tabular(enhance=enhance)
        # elif 'gridded' in netcdf_type:
        #     xds = self.to_xarray_gridded(range_minmax=range_minmax, bearing=bearing, enhance=enhance)

        # Check if dataset has distance_from_origin in coordinates. We will prepend the .nc extension
        # with the appropriate name depending on whether the radial file is tabular or gridded.
        if prepend_extension:
            if "bearing" in xds.coords:
                pre_ext = "gridded"
            else:
                pre_ext = "tabular"
            # Change the extension to reflect the type of wave file
            filename = filename.with_suffix(f".{pre_ext}.nc")

        # encoding = make_encoding(xds, comp_level=4, fillvalue=np.nan)

        # # Assign header data to global attributes
        # xds['site'] = self.metadata['Site'].strip('"').strip()
        # xds['site'] = xds['site'].assign_attrs(self.metadata)

        # Grab min and max time in dataset for entry into global attributes for cf compliance
        time_start = xds['time'].min().data
        time_end = xds['time'].max().data

        if user_attributes:
            from hfradarpy.io.nc import required_global_attributes
            global_attributes = required_global_attributes(user_attributes, time_start, time_end)

            global_attributes['geospatial_lat_min'] = np.double(xds.lat.min())
            global_attributes['geospatial_lat_max'] = np.double(xds.lat.max())
            global_attributes['geospatial_lon_min'] = np.double(xds.lon.min())
            global_attributes['geospatial_lon_max'] = np.double(xds.lon.max())

            logging.debug('{} - Assigning global attributes to dataset'.format(self.file_name))
            xds = xds.assign_attrs(global_attributes)

        # Encode and compress variables 
        encoding = make_encoding(xds, comp_level=4, fillvalue=np.nan)
        
        if "gridded" in model:
            encoding["bearing"] = dict(zlib=False, _FillValue=None)
            encoding["range"] = dict(zlib=False, _FillValue=None)
        # encoding["time"] = dict(zlib=False, _FillValue=None)

        xds.to_netcdf(filename, encoding=encoding, format="netCDF4", engine="netcdf4", unlimited_dims=["time"])

    def to_ruv(self, filename, validate=True, overwrite=False):
        """
        Create a CODAR Radial (.ruv) file from radial instance

        Args:
            filename (str or Path): User defined filename of radial file you want to save
            validate (boolean): If False, no validation check will be performed before creating the file.
            overwrite (bool): If True, an exported file can overwrite an existing file with the same name. Defaults to False.

        """
        # Make sure filename is converted into a Path object
        filename = Path(filename)

        if validate:
                if not self.is_valid():
                    raise ValueError("Could not export ASCII data, the input file was invalid.")

        # Ensure that the filename passed into the export function is not the same as the filename that we read in.
        # # We do not want to overwrite the original wave file by accident.
        if not overwrite:
            if self.full_file == str(filename):
                suffix = f'.mod{filename.suffix}'
                filename = filename.with_suffix(suffix)

        if os.path.isfile(filename):
            os.remove(filename)

        os.makedirs(os.path.dirname(filename), exist_ok=True)
        rcopy = copy.deepcopy(self)
        with open(filename, "w") as f:
            # Write header
            for metadata_key, metadata_value in self.metadata.items():
                if "ProcessedTimeStamp" in metadata_key:
                    break
                else:
                    print(metadata_key)
                    f.write("%{}: {}\n".format(metadata_key, metadata_value))

            # Write data tables. Anything beyond the first table is commented out.
            for table in self._tables.keys():

                if "datetime" in self._tables[table]["data"].keys():
                    self._tables[table]["data"] = self._tables[table]["data"].drop(["datetime"], axis=1)

                for table_key, table_value in self._tables[table].items():
                    if table_key != 'data':
                        if (table_key == 'TableType') & (table == '1'):
                            if 'QCD' in self.metadata:
                               for qcd_info in self.metadata['QCD']:
                                    f.write('%{}\n'.format(qcd_info))
                            if 'QCTest' in self.metadata:
                                f.write('%QCFileVersion: 1.0.0\n')
                                f.write('%QCReference: Quality control reference: IOOS QARTOD HF Radar ver 1.0 May 2016\n')
                                f.write('%QCFlagDefinitions: 1=pass 2=not_evaluated 3=suspect 4=fail 9=missing_data\n')
                                f.write('%QCTestFormat: "test_name [qc_thresholds]: test_result"\n')

                                for test in self.metadata["QCTest"]:
                                    f.write("%QCTest: {}\n".format(test))
                            f.write("%{}: {}\n".format(table_key, table_value))
                        elif table_key == "TableColumns":
                            f.write("%TableColumns: {}\n".format(len(self._tables[table]["data"].columns)))
                        elif table_key == "TableColumnTypes":
                            f.write("%TableColumnTypes: {}\n".format(" ".join(self._tables[table]["data"].columns.to_list())))
                        elif table_key == "TableStart":
                            f.write("%{}: {}\n".format(table_key, table_value))
                        elif table_key == "_TableHeader":
                            pass
                        else:
                            f.write("%{}: {}\n".format(table_key, table_value))

                if table == 1:
                    # Fill NaN with 999.000 which is the standard fill value for codar lluv files
                    self.data = self.data.fillna(999.000)

                    try:
                        self.data["LOND"] = self.data["LOND"].apply(lambda x: "{:.7f}".format(x))
                        self.data["LATD"] = self.data["LATD"].apply(lambda x: "{:.7f}".format(x))
                        self.data["ESPC"] = self.data["ESPC"].apply(lambda x: "{:.3f}".format(x))
                        if "ETMP" in self.data.columns:
                            self.data["ETMP"] = self.data["ETMP"].apply(lambda x: "{:.3f}".format(x))
                        self.data["BEAR"] = self.data["BEAR"].apply(lambda x: "{:.1f}".format(x))
                        self.data["HEAD"] = self.data["HEAD"].apply(lambda x: "{:.1f}".format(x))
                    except:
                        self = rcopy
                        print("Unexpected error in formatting one of these columns: LOND LATD ESPC ETMP BEAR HEAD")

                    # Convert _TableHeader to a new dataframe and concatenate to dataframe containing radial data
                    # This allows for the output format to follow CODARS CTF specifications
                    # The below block of code adds the weird header and units format that codar uses in their files
                    row_df = pd.DataFrame([self._tables[1]["_TableHeader"][1]], columns=self._tables[1]["_TableHeader"][0])
                    self.data.columns = self._tables[1]["_TableHeader"][0]
                    self.data = pd.concat([row_df, self.data], ignore_index=True)
                    self.data.insert(0, "%%", np.nan)  # Insert column at the beginning of dataframe of NaNs
                    self.data.iloc[0, self.data.columns.get_loc("%%")] = "%%"  # make the first row in the first column a '%%'

                    # Output data table to string
                    # self.data.to_string(f, index=False, justify='center', header=True, na_rep=' ')
                    self.data.temp = re.sub(
                        " %%", "%%", self.data.to_string(index=False, justify="right", header=True, na_rep=" ")
                    )
                    f.write(self.data.temp)
                else:
                    table_alias = self._tables[table]
                    table_alias["data"].insert(0, "%%", "%")
                    table_alias["data"] = table_alias["data"].fillna(999.000)
                    if table_alias["TableType"] == "rads rad1":
                        f.write(
                            "%%   Time       Calculated      Calculated     Corrected        Noise Floor     SignalToNoise   Diag  Valid  Dual  Radial RadsV Rads   Max   Vel    Vel    Bearing  Radial Spectra Time\n"
                        )
                        f.write(
                            "%% FromStart   Amp1    Amp2   Phase13 Phase23 Phase1 Phase2   NF1   NF2   NF3   SN1  SN2  SN3   Range Dopplr Angle Vector  per  Range Range  Max    Aver   Average   Type   Type   Year Mo Dy  Hr Mn  S\n"
                        )
                        f.write(
                            "%%  Seconds  (1/v^2) (1/v^2)   (deg)   (deg)   (deg)  (deg)  (dBm) (dBm) (dBm)  (dB) (dB) (dB)  Cell  Cells  Prcnt Count  Range Cells  (km) (cm/s) (cm/s) (deg CWN)\n"
                        )
                        table_alias["data"].to_string(f, index=False, justify="center", header=False)
                    elif table_alias["TableType"] == "rcvr rcv4":
                        f.write(
                            "%%   Minutes  Rcvr  Awg XmitTrip  AwgRun   Supply  +5VDC  -5VDC +12VDC XInt XAmp XForw XRefl  Xmit  X+Ampl  X+5V  X2Int X2Amp X2Forw X2Refl  Xmit2 X2+Amp  X2+5V  GpsRcv GpsDsp GpsSat GpsSat   PLL   HiRcvr  Humid Supply Extern Extern CompRunTime    Date                \n"
                        )
                        f.write(
                            "%%  FromStart degC degC  HexCode  Seconds   Volts  Volts  Volts  Volts degC degC Watts Watts  VSWR   Volts  Volts  degC  degC  Watts  Watts  VSWR   Volts  Volts   Mode   Mode   Lock  Unlock Unlock   degC     %    Amps  InputA InputB   Minutes      Year Mo Dy Hr Mn Sec\n"
                        )
                        table_alias["data"].to_string(f, index=False, justify="center", header=False)
                    elif table_alias["TableType"] == "rcvr rcv3":
                        f.write(
                            "%%   Minutes  Rcvr  Awg XmitTrip  AwgRun   Supply  +5VDC  -5VDC +12VDC XInt XAmp XForw XRefl  Xmit  X+Ampl  X+5V   GpsRcv GpsDsp GpsSat GpsSat   PLL   HiRcvr  Humid Supply Extern Extern CompRunTime    Date\n"
                        )
                        f.write(
                            "%%  FromStart degC degC  HexCode  Seconds   Volts  Volts  Volts  Volts degC degC Watts Watts  VSWR   Volts  Volts   Mode   Mode   Lock  Unlock Unlock   degC     %    Amps  InputA InputB   Minutes      Year Mo Dy Hr Mn Sec\n"
                        )
                        table_alias["data"].to_string(f, index=False, justify="center", header=False)
                    elif table_alias["TableType"] == "rcvr rcv2":
                        f.write(
                            "%%   LogTime Rcvr Awg3 XmtTrip Awg3Run Supply +5VDC -5VDC +12VDC XInt XAmp XForw XRefl X+Ampl X+5VDC GpsRcv GpsDsp  PLL   HiRcvr Humid Supply Extern Extern  CompRunTime   Year Mo Dy Hr Mn  S\n"
                        )
                        f.write(
                            "%%   Minutes   °C   °C HexCode Seconds Volts  Volts Volts  Volts   °C   °C Watts Watts  Volts  Volts  Mode   Mode  Unlock    °C     %   Amps  InputA InputB    Minutes                        \n"
                        )
                        table_alias["data"].to_string(f, index=False, justify="center", header=False)
                    else:
                        table_alias["data"].to_string(f, index=False, justify="center", header=True)

                if table > 1:
                    f.write("\n%TableEnd: {}\n".format(table))
                else:
                    f.write("\n%TableEnd: \n")
                f.write("%%\n")

            # Write footer containing processing information
            f.write("%ProcessedTimeStamp: {}\n".format(self.metadata["ProcessedTimeStamp"]))
            for tool in self.metadata["ProcessingTool"]:
                f.write("%ProcessingTool: {}\n".format(tool))
                # f.write('%{}: {}\n'.format(footer_key, footer_value))
            f.write("%End:")

    def initialize_qc(self):
        """
        Initialize an empty list of QC tests
        """
        self.metadata["QCTest"] = []

    # QARTOD QC Tests
    def qc_qartod_avg_radial_bearing(self, reference_bearing, warning_threshold=15, failure_threshold=30):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Valid Location (Test 12)
        Check that the average radial bearing remains relatively constant (Roarty et al. 2012).

        It is expected that the average of all radial velocity bearings AVG_RAD_BEAR obtained during a sample
        interval (e.g., 1 hour) should be close to a reference bearing REF_RAD_BEAR and not vary beyond warning
        or failure thresholds.

        Args:
            reference_bearing (int or float): Reference bearing
            warning_threshold (int, optional): Warning Threshold. Defaults to 15.
            failure_threshold (int, optional): Failure Threshold. Defaults to 30.
        """
        test_str = "QC12"
        # Absolute value of the difference between the bearing mean and reference bearing
        absolute_difference = np.abs(self.data["BEAR"].mean() - reference_bearing)

        if absolute_difference >= failure_threshold:
            flag = 4
        elif (absolute_difference >= warning_threshold) & (absolute_difference < failure_threshold):
            flag = 3
        elif absolute_difference < warning_threshold:
            flag = 1

        self.data[test_str] = flag  # Assign the flags to the column
        self.metadata["QCTest"].append(
            (
                f"qc_qartod_avg_radial_bearing ({test_str}) - Test applies to entire file. Thresholds="
                "[ "
                f"reference_bearing={reference_bearing} (degrees) "
                f"warning={warning_threshold} (degrees) "
                f"failure={failure_threshold} (degrees) "
                f"]: See result in column {test_str} below"
            )
        )
        self.append_to_tableheader(test_str, "(flag)")

    def qc_qartod_valid_location(self, use_mask=False):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Valid Location (Test 8)
        Removes radial vectors placed over land or in other unmeasureable areas

        Radial vector coordinates are checked against a reference file containing information about which locations
        are over land or in an unmeasurable area (for example, behind an island or point of land). Radials in these
        areas will be flagged with a code (FLOC) in the radial file (+128 in CODAR radial files) and are not included
        in total vector calculations.

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            use_mask (bool, optional): Use mask_over_land function in addition to manufacturers flags. Defaults to False.
        """

        test_str = "QC08"
        flag_column = "VFLG"

        if flag_column in self.data:
            self.data[test_str] = 1  # add new column of passing values
            self.data.loc[(self.data[flag_column] == 128), test_str] = 4  # set to 4 where land is flagged (manufacturer)
            if use_mask:
                self.data.loc[~self.mask_over_land(), test_str] = 4  # set to 4 where land is flagged (mask_over_land)
            self.metadata["QCTest"].append(
                (
                    f"qc_qartod_valid_location ({test_str}) - Test applies to each row. Thresholds=[{flag_column}==128]: "
                    f"See results in column {test_str} below"
                )
            )
            self.append_to_tableheader(test_str, "(flag)")

        else:
            logger.warning(f"qc_qartod_valid_location not run, no {flag_column} column")

    def qc_qartod_radial_count(self, min_count=150, low_count=300):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Radial Count (Test 9)
        Rejects radials in files with low radial counts (poor radial map coverage).

        The number of radials (RCNT) in a radial file must be above a threshold value RCNT_MIN to pass the test and
        above a value RC_LOW to not be considered suspect. If the number of radials is below the minimum level,
        it indicates a problem with data collection. In this case, the file should be rejected and none of the
        radials used for total vector processing.

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            min_count (int, optional):
                Minimum radial count threshold (failure) below which the file should be rejected. Defaults to 150.
            low_count (int, optional):
                Low radial count threshold (warning) below which the file should be considered suspect. Defaults to 300.
        """
        test_str = "QC09"
        column_flag = "VFLG"

        # If a vector flag is supplied by the vendor, subset by that first
        if column_flag in self.data:
            num_radials = len(self.data[self.data[column_flag] != 128])
        else:
            num_radials = len(self.data)

        if num_radials < min_count:
            radial_count_flag = 4
        elif (num_radials >= min_count) and (num_radials <= low_count):
            radial_count_flag = 3
        elif num_radials > low_count:
            radial_count_flag = 1

        self.data[test_str] = radial_count_flag
        self.metadata["QCTest"].append(
            (
                f"qc_qartod_radial_count ({test_str}) - Test applies to entire file. Thresholds="
                "[ "
                f"failure={min_count} (radials) "
                f"warning_num={low_count} (radials) "
                f"<valid_radials={num_radials}> "
                f"]:  See results in column {test_str} below"
            )
        )
        self.append_to_tableheader(test_str, "(flag)")

    def qc_qartod_maximum_velocity(self, max_speed=250, high_speed=150):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Max Threshold (Test 7)
        Ensures that a radial current speed is not unrealistically high.

        The maximum radial speed threshold (RSPDMAX) represents the maximum reasonable surface radial velocity
        for the given domain.

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            max_speed (int, optional):
                Maximum Radial Speed (cm/s). Radials beyond this speed will be flagged a failure. Defaults to 250
            high_speed (int, optional):
                High Radial Speed (cm/s). Radials between high and max speed will be flagged suspect. Defaults to 150
        """
        test_str = "QC07"

        self.data["VELO"] = self.data["VELO"].astype(float)  # make sure VELO is a float

        # Add new column to dataframe for test, and set every row as passing, 1, flag
        self.data[test_str] = 1

        # velocity is less than radial_max_speed but greater than radial_high_speed, set that row as a warning, 3, flag
        self.data.loc[
            (self.data["VELO"].abs() < max_speed) & (self.data["VELO"].abs() > high_speed), test_str
        ] = 3

        # if velocity is greater than radial_max_speed, set that row as a fail, 4, flag
        self.data.loc[(self.data["VELO"].abs() > max_speed), test_str] = 4

        self.metadata["QCTest"].append(
            (
                f"qc_qartod_maximum_velocity ({test_str}) - Test applies to each row. Thresholds="
                "[ "
                f"high_vel={str(high_speed)} (cm/s) "
                f"max_vel={str(max_speed)} (cm/s) "
                f"]: See results in column {test_str} below"
            )
        )

        self.append_to_tableheader(test_str, "(flag)")

    def qc_qartod_spatial_median(
        self, smed_range_cell_limit=2.1, smed_angular_limit=10, smed_current_difference=30
    ):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Spatial Median (Test 10)

        Ensures that the radial velocity is not too different from nearby radial velocities.
        RV is the radial velocity
        NV is a set of radial velocities for neighboring radial cells
        (cells within radius of 'radial_smed_range_cell_limit' * Range Step (km)
        and whose vector bearing (angle of arrival at site)
        is also within 'radial_smed_angular_limit' degrees of the source vector's bearing)
        Required to pass the test: |RV - median(NV)| <= radial_smed_current_difference
        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            radial_smed_range_cell_limit (float, optional):
                Multiple of range step which depends on the radar type. Defaults to 2.1.
            radial_smed_angular_limit (int, optional):
                Limit for number of degrees from source radial's bearing (degrees). Defaults to 10.
            radial_smed_current_difference (int, optional):
                Current difference (cm/s). Defaults to 30.
        """
        test_str = "QC10"

        self.data[test_str] = 1
        try:
            Rstep = float(self.metadata["RangeResolutionKMeters"])
            # Rstep = np.floor(min(np.diff(np.unique(self.data['RNGE'])))) #use as backup method if other fails?

            Bstep = [float(s) for s in re.findall(r"-?\d+\.?\d*", self.metadata["AngularResolution"])]
            Bstep = Bstep[0]
            # Bstep = int(min(np.diff(np.unique(self.data['BEAR']))))  #use as backup method if other fails?

            RLim = int(round(smed_range_cell_limit))  # if not an integer will cause an error later on
            BLim = int(smed_angular_limit / Bstep)  # if not an integer will cause an error later on

            # convert bearing into bearing cell numbers
            adj = np.mod(min(self.data["BEAR"]), Bstep)
            Bcell = ((self.data["BEAR"] - adj) / Bstep) - 1
            Bcell = Bcell.astype(int)
            # Btable = np.column_stack((self.data['BEAR'], Bcell))  #only for debugging

            # convert range into range cell numbers
            Rcell = np.floor((self.data["RNGE"] / Rstep) + 0.1)
            Rcell = Rcell - min(Rcell)
            Rcell = Rcell.astype(int)
            # Rtable = np.column_stack((self.data['RNGE'], Rcell))   #only for debugging
            Rcell = self.data["SPRC"]

            # place velocities into a matrix with rows defined as bearing cell# and columns as range cell#
            BRvel = np.zeros((int(360 / Bstep), max(Rcell) + 1), dtype=int) + np.nan
            BRind = np.zeros((int(360 / Bstep), max(Rcell) + 1), dtype=int) + np.nan

            for xx in range(len(self.data["VELO"])):
                BRvel[Bcell[xx]][Rcell[xx]] = self.data["VELO"][xx]
                BRind[Bcell[xx]][Rcell[xx]] = xx  # keep track of indices so easier to return to original format

            # deal with 359 to 0 transition in bearing by
            # repeating first BLim rows at the bottom and last BLim rows at the top
            # also pad ranges with NaNs by adding extra columns on the left and right of the array
            # this keeps the indexing for selecting the neighbors from breaking

            BRtemp = np.append(np.append(BRvel[-BLim:], BRvel, axis=0), BRvel[:BLim], axis=0)
            rangepad = np.zeros((BRtemp.shape[0], RLim), dtype=int) + np.nan
            BRpad = np.append(np.append(rangepad, BRtemp, axis=1), rangepad, axis=1)

            # calculate median of neighbors (neighbors include the point itself)
            BRmed = BRpad + np.nan  # initialize with an array of NaN
            for rr in range(RLim, BRvel.shape[1] + RLim):
                for bb in range(BLim, BRvel.shape[0] + BLim):
                    temp = BRpad[bb - BLim: bb + BLim + 1, rr - RLim: rr + RLim + 1]  # temp is the matrix of neighbors
                    BRmed[bb][rr] = np.nanmedian(temp)

            # now remove the padding from the array containing the median values
            BRmedtrim = BRmed[BLim:-BLim, RLim:-RLim]

            # calculate velocity minus median of neighbors
            # and put back into single column using the indices saved in BRind
            BRdiff = (
                BRvel - BRmedtrim
            )  # velocity minus median of neighbors, test these values against current radial_smed_current_difference
            diffcol = self.data["RNGE"] + np.nan  # initialize a single column for the difference results
            for rr in range(BRdiff.shape[1]):
                for bb in range(BRdiff.shape[0]):
                    if not (np.isnan(BRind[bb][rr])):
                        diffcol[BRind[bb][rr]] = BRdiff[bb][rr]
            boolean = diffcol.abs() > smed_current_difference

        except TypeError:
            diffcol = diffcol.astype(float)
            boolean = diffcol.abs() > smed_current_difference

        self.data[test_str] = self.data[test_str].where(~boolean, other=4)
        self.metadata["QCTest"].append(
            (
                f"qc_qartod_spatial_median ({test_str}) - Test applies to each row. Thresholds="
                "[ "
                f"range_cell_limit={str(smed_range_cell_limit)} (range cells) "
                f"angular_limit={str(smed_angular_limit)} (degrees) "
                f"current_difference={str(smed_current_difference)} (cm/s) "
                f"]: See results in column {test_str} below"
            )
        )
        self.append_to_tableheader(test_str, "(flag)")

    def qc_qartod_syntax(self):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Syntax (Test 6)

        This test is required to be QARTOD compliant.

        A collection of tests ensuring proper formatting and existence of fields within a radial file.

        The radial file may be tested for proper parsing and content, for file format (hfrweralluv1.0, for example),
        site code, appropriate time stamp, site coordinates, antenna pattern type (measured or ideal, for DF
        systems), and internally consistent row/column specifications.

        ----------------------------------------------------------------------------------------------------------------------
        Fail: One or more fields are corrupt or contain invalid data, If “File Format” ≠ “hfrweralluv1.0”, flag = 4

        Pass: Applies for test pass condition.
        ----------------------------------------------------------------------------------------------------------------------
        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/
        :param threshold: Maximum Radial Speed (cm/s)
        """
        test_str = "QC06"

        i = 0

        # check for timestamp in filename
        result = re.search(r"\d{4}_\d{2}_\d{2}_\d{4}", self.file_name)
        if result:
            timestr = dt.datetime.strptime(result.group(), "%Y_%m_%d_%H%M")
            i = i + 1

        # Radial tables must not be empty
        if self.is_valid():
            i = i + 1

        # The following metadata must be defined.
        tmp = self.metadata
        if (
            tmp["FileType"] and tmp["Site"] and tmp["TimeStamp"] and tmp["Origin"] and tmp["PatternType"] and tmp["TimeZone"]
        ):
            filetime = dt.datetime(*map(int, self.metadata["TimeStamp"].split()))
            i = i + 1

        # Filename timestamp must match the timestamp reported within the file.
        if timestr == filetime:
            i = i + 1

        # Radial data table columns stated must match the number of columns reported for each row
        if len(self._tables[1]["TableColumnTypes"].split()) == self.data.shape[1]:
            i = i + 1

        # Make sure site location is within range: -180 <= lon <= 180 & -90 <= lat <= 90
        latlon = re.findall(r"[-+]?\d*\.\d+|\d+", self.metadata["Origin"])
        if (-180 <= float(latlon[1]) <= 180) & (-90 <= float(latlon[0]) <= 90):
            i = i + 1

        if i == 6:
            syntax = 1
        else:
            syntax = 4

        self.data[test_str] = syntax
        self.metadata["QCTest"].append(
            f"qc_qartod_syntax ({test_str}) - Test applies to entire file. Thresholds=[N/A]: See results in column {test_str}"
        )
        self.append_to_tableheader(test_str, "(flag)")

    def qc_qartod_temporal_gradient(self, r0, gradient_temp_fail=54, gradient_temp_warn=36):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Temporal Gradient (Test 11)
        Checks for satisfactory temporal rate of change of radial components

        Test determines whether changes between successive radial velocity measurements at a particular range
        and bearing cell are within an acceptable range. GRADIENT_TEMP = |Rt-1 - Rt|

        Flags Condition Codable Instructions
        Fail = 4 The temporal change between successive radial velocities exceeds the gradient failure threshold.

        If GRADIENT_TEMP ≥ GRADIENT_TEMP_FAIL,
        flag = 4

        Suspect = 3 The temporal change between successive radial velocities is less than the gradient failure threshold but
        exceeds the gradient warn threshold.

        If GRADIENT_TEMP < GRADIENT_TEMP_FAIL & GRADIENT_TEMP ≥ GRADIENT_TEMP_WARN,
        flag = 3

        Pass = 1 The temporal change between successive radial velocities is less than the gradient warn threshold.

        If GRADIENT_TEMP < GRADIENT_TEMP_WARN,
        flag = 1

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            r0 (str or Path): Full path to the filename of the previous hourly radial.
            gradient_temp_fail (int, optional): Maximum Radial Speed (cm/s). Defaults to 54.
            gradient_temp_warn (int, optional): Warning Radial Speed (cm/s). Defaults to 36.
        """
        test_str = "QC11"
        # self.data[test_str] = data
        self.metadata["QCTest"].append(
            (
                f"qc_qartod_temporal_gradient ({test_str}) - Test applies to each row. Thresholds="
                "[ "
                f"gradient_temp_warn={str(gradient_temp_warn)} (cm/s*hr) "
                f"gradient_temp_fail={str(gradient_temp_fail)} (cm/s*hr) "
                f"]: See results in column {test_str} below"
            )
        )
        self.append_to_tableheader(test_str, "(flag)")

        if os.path.exists(r0):
            r0 = Radial(r0)

            if r0.is_valid():
                
                merged = self.data.merge(r0.data, on=["LOND", "LATD"], how="left", suffixes=(None, "_x"), indicator="Exist")
                difference = (merged["VELO"] - merged["VELO_x"]).abs()

                # Add new column to dataframe for test, and set every row as passing, 1, flag
                self.data[test_str] = 1

                # If any point in the recent radial does not exist in the previous radial, set row as a not evaluated, 2, flag
                self.data.loc[merged['Exist'] == 'left_only', test_str] = 2

                # velocity is less than radial_max_speed but greater than radial_high_speed, set row as a warning, 3, flag
                self.data.loc[(difference < gradient_temp_fail) & (difference > gradient_temp_warn), test_str] = 3

                # if velocity is greater than radial_max_speed, set that row as a fail, 4, flag
                self.data.loc[(difference > gradient_temp_fail), test_str] = 4
            else:
                # Add new column to dataframe for test, and set every row as not_evaluated, 2, flag
                self.data[test_str] = 2
                logging.warning('{} is corrupt or contains no data. Setting column {} to not_evaluated flag'.format(r0, test_str))

        else:
            # Add new column to dataframe for test, and set every row as not_evaluated, 2, flag
            self.data[test_str] = 2
            logging.warning(
                "{} does not exist at specified location. Setting column {} to not_evaluated flag".format(r0, test_str)
            )

    def qc_qartod_primary_flag(self, include=None):
        """
        A primary flag is a single flag set to the worst case of all QC flags within the data record.

        Args:
            include (list, optional):
                list of quality control tests which should be included in the primary flag.
                Defaults to None, which includes all tests.
        """
        test_str = "PRIM"

        # Set summary flag column all equal to 1
        self.data[test_str] = 1

        # generate dictionary of executed qc tests found in the header
        executed = dict()
        for b in [x.split("-")[0].strip() for x in self.metadata["QCTest"]]:
            i = b.split(" ")
            executed[i[0]] = re.sub(r"[()]", "", i[1])

        if include:
            # only add qartod tests which were set by user to executed dictionary
            included_tests = list({key: value for key, value in executed.items() if key in include}.values())
        else:
            included_tests = list(executed.values())

        equals_3 = self.data[included_tests].eq(3).any(axis=1)
        self.data[test_str] = self.data[test_str].where(~equals_3, other=3)

        equals_4 = self.data[included_tests].eq(4).any(axis=1)
        self.data[test_str] = self.data[test_str].where(~equals_4, other=4)

        included_test_strs = ", ".join(included_tests)
        self.metadata["QCTest"].append(
            (
                f'qc_qartod_primary_flag ({test_str}) - Primary Flag - Highest flag value of {included_test_strs}' +
                '("not_evaluated" flag results ignored)'
            )
        )
        self.append_to_tableheader(test_str, "(flag)")
        # %QCFlagDefinitions: 1=pass 2=not_evaluated 3=suspect 4=fail 9=missing_data

    def append_to_tableheader(self, test_string, test_unit):
        """
        Append key, value of metadata in the radial header and footer to the _TableHeader

        Args:
            test_string (str): key
            test_unit (str): value
        """
        self._tables[1]["_TableHeader"][0].append(test_string)
        self._tables[1]["_TableHeader"][1].append(test_unit)

    def reset(self):
        """
        Reset data variable of object, r.data, back to original dataset
        """
        logging.info("Resetting instance data variable to original dataset")
        self.data = self._tables[1]["data"]


if __name__ == "__main__":
    f = "/Users/mikesmith/Documents/github/rucool/hfradarpy/examples/data/radials/ruv/SEAB/RDLi_SEAB_2019_01_01_0100.ruv"
    r = Radial(f, replace_invalid=True, vflip=True)
    r.mask_over_land()
    ds = r.to_xarray('gridded')
    print(ds)
