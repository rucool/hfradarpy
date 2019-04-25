import datetime as dt
import geopandas as gpd
import logging
import numpy as np
import os
import pandas as pd
import re
from shapely.geometry import Point
import xarray as xr
from codar_processing.common import CTFParser, create_dir, make_encoding
from codar_processing.calc import reckon


logger = logging.getLogger(__name__)


def concatenate_radials(radial_list):
    """
    This function takes a list of radial files. Loads them all separately using the Radial object and then combines
    them along the time dimension using xarrays built-in concatenation routines.
    :param radial_list: list of radial files that you want to concatenate
    :return: radial files concatenated into an xarray dataset by range, bearing, and time
    """

    radial_dict = {}
    for each in sorted(radial_list):
        radial = Radial(each, multi_dimensional=True)
        radial_dict[radial.file_name] = radial.data

    ds = xr.concat(radial_dict.values(), 'time')
    return ds


class Radial(CTFParser):
    """
    Radial Subclass.

    This class should be used when loading a CODAR radial (.ruv) file. This class utilizes the generic LLUV class from
    ~/codar_processing/common.py in order to load CODAR Radial files
    """
    def __init__(self, fname, replace_invalid=True, multi_dimensional=False, mask_over_land=False):
        keep = ['LOND', 'LATD', 'VELU', 'VELV', 'VFLG', 'ESPC', 'ETMP', 'MAXV', 'MINV', 'ERSC', 'ERTC', 'XDST', 'YDST', 'RNGE', 'BEAR', 'VELO', 'HEAD', 'SPRC']

        CTFParser.__init__(self, fname)
        for key in self._tables.keys():
            table = self._tables[key]
            if 'LLUV' in table['TableType']:
                self.data = table['data']
            elif 'rads' in table['TableType']:
                self.diagnostics_radial = table['data']
                self.diagnostics_radial['datetime'] = self.diagnostics_radial[['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)
            elif 'rcvr' in table['TableType']:
                self.diagnostics_hardware = table['data']
                self.diagnostics_hardware['datetime'] = self.diagnostics_hardware[['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)

        if replace_invalid:
            self.replace_invalid_values()

        if mask_over_land:
            land = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
            land = land[land['continent'] == 'North America']
            # ocean = gpd.read_file('/Users/mikesmith/Downloads/ne_10m_ocean')
            self.data = gpd.GeoDataFrame(self.data, crs={'init': 'epsg:4326'}, geometry=[Point(xy) for xy in zip(self.data.LOND.values, self.data.LATD.values)])

            # Join the geodataframe containing radial points with geodataframe containing leasing areas
            self.data = gpd.tools.sjoin(self.data, land, how='left')

            # All data in the continent column that lies over water should be nan.
            self.data = self.data[keep][self.data['continent'].isna()]
            self.data = self.data.reset_index()

        if multi_dimensional:
            self.to_multi_dimensional()

    def __repr__(self):
        return "<Radial: {}>".format(self.file_name)

    def to_multi_dimensional(self, range_max=66.4466):
        """
        Adapted from MATLAB code from Mark Otero
        http://cordc.ucsd.edu/projects/mapping/documents/HFRNet_Radial_NetCDF.pdf
        :param range_min:
        :param range_max:
        :return:
        """
        # Clean radial header
        # self.clean_header()

        # CF Standard: T, Z, Y, X
        coords = ('time', 'range', 'bearing')

        # Intitialize empty xarray dataset
        ds = xr.Dataset()

        range_step = float(self.metadata['RangeResolutionKMeters'].split()[0])
        # range_dim = np.arange(range_min, range_max + float(self.metadata['RangeResolutionKMeters']), float(self.metadata['RangeResolutionKMeters']))
        range_dim = np.arange(self.data.RNGE.min(), np.round(range_max + range_step, 4), range_step)
        bearing_dim = np.arange(1, 361, 1).astype(np.float)  # Complete 360 degree bearing coordinate allows for better aggregation

        # create radial grid from bearing and range
        [bearing, range] = np.meshgrid(bearing_dim, range_dim)

        # calculate lat/lons from origin, bearing, and ranges
        latlon = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", self.metadata['Origin'])]
        latd, lond = reckon(latlon[0], latlon[1], bearing, range)

        # create dictionary containing variables from dataframe in the shape of radial grid
        d = {key: np.tile(np.nan, bearing.shape) for key in self.data.keys()}

        # find grid indices from radial grid (bearing, ranges)
        range_map_idx = np.tile(np.nan, self.data['RNGE'].shape)
        bearing_map_idx = np.tile(np.nan, self.data['BEAR'].shape)

        for i, line in enumerate(self.data['RNGE']):
            range_map_idx[i] = np.argmin(np.abs(range_dim - self.data.RNGE[i]))
            bearing_map_idx[i] = np.argmin(np.abs(bearing_dim - self.data.BEAR[i]))

        for k, v in d.items():
            v[range_map_idx.astype(int), bearing_map_idx.astype(int)] = self.data[k]
            d[k] = v

        # Add extra dimension for time
        d = {k: np.expand_dims(np.float32(v), axis=0) for (k, v) in d.items()}

        ds['lon'] = (('range', 'bearing'), lond.round(4))
        ds['lat'] = (('range', 'bearing'), latd.round(4))

        # Add coordinate variables to dataset
        timestamp = dt.datetime(*[int(s) for s in self.metadata['TimeStamp'].split()])
        ds.coords['time'] = pd.date_range(timestamp, periods=1)
        ds.coords['bearing'] = bearing_dim
        ds.coords['range'] = range_dim

        # Add all variables to dataset
        for k, v in d.items():
            ds[k] = (coords, v)

        # Check if calculated longitudes and latitudes align with given longitudes and latitudes
        # plt.plot(ds.lon, ds.lat, 'bo', ds.LOND.squeeze(), ds.LATD.squeeze(), 'rx')

        # Drop extraneous variables
        ds = ds.drop(['LOND', 'LATD', 'BEAR', 'RNGE'])

        # Flip sign so positive velocities are away from the radar as per cf conventions
        ds['MINV'] = -ds.MINV
        ds['MAXV'] = -ds.MAXV
        ds['VELO'] = -ds.VELO

        # Assign header data to global attributes
        self.data = ds.assign_attrs(self.metadata)

    def file_type(self):
        """Return a string representing the type of file this is."""
        return 'radial'

    def clean_header(self, split_origin=False):
        """
        Clean up the radial header dictionary so that you can upload it to the HFR MySQL Database.
        :return:
        """
        keep = ['CTF', 'FileType', 'LLUVSpec', 'UUID', 'Manufacturer', 'Site', 'TimeStamp', 'TimeZone', 'TimeCoverage',
                'Origin', 'GreatCircle', 'GeodVersion', 'LLUVTrustData', 'RangeStart', 'RangeEnd', 'RangeResolutionKMeters',
                'AntennaBearing', 'ReferenceBearing', 'AngularResolution', 'SpatialResolution', 'PatternType', 'PatternDate',
                'PatternResolution', 'TransmitCenterFreqMHz', 'DopplerResolutionHzPerBin', 'FirstOrderMethod',
                'BraggSmoothingPoints', 'CurrentVelocityLimits', 'BraggHasSecondOrder', 'RadialBraggPeakDropOff',
                'RadialBraggPeakNull', 'RadialBraggNoiseThreshold', 'PatternAmplitudeCorrections', 'PatternPhaseCorrections',
                'PatternAmplitudeCalculations', 'PatternPhaseCalculations', 'RadialMusicParameters', 'MergedCount',
                'RadialMinimumMergePoints', 'FirstOrderCalc', 'MergeMethod', 'PatternMethod', 'TransmitSweepRateHz',
                'TransmitBandwidthKHz', 'SpectraRangeCells', 'SpectraDopplerCells', 'TableType', 'TableColumns', 'TableColumnTypes',
                'TableRows', 'TableStart', 'CurrentVelocityLimit']

        # TableColumnTypes
        key_list = list(self.metadata.keys())
        for key in key_list:
            if not key in keep:
                del self.metadata[key]

        for k, v in self.metadata.items():
            if 'Site' in k:
                self.metadata[k] = ''.join(e for e in v if e.isalnum())
            elif k in ('TimeStamp', 'PatternDate'):
                t_list = [int(s) for s in v.split()]
                self.metadata[k] = dt.datetime(*t_list)
            elif 'TimeZone' in k:
                self.metadata[k] = v.split('"')[1]
            elif 'TableColumnTypes' in k:
                self.metadata[k] = ' '.join([x.strip() for x in v.strip().split(' ')])
            elif 'Origin' in k:
                if split_origin:
                    self.metadata[k] = re.findall(r"[-+]?\d*\.\d+|\d+", v)
                else:
                    self.metadata[k] = v.lstrip()
            elif k in ('RangeStart', 'RangeEnd', 'AntennaBearing', 'ReferenceBearing', 'AngularResolution', 'SpatialResolution',
                       'FirstOrderMethod', 'BraggSmoothingPoints', 'BraggHasSecondOrder', 'MergedCount',
                       'RadialMinimumMergePoints', 'FirstOrderCalc', 'SpectraRangeCells', 'SpectraDopplerCells',
                       'TableColumns', 'TableRows',  'PatternResolution', 'CurrentVelocityLimit', 'TimeCoverage'):
                try:
                    self.metadata[k] = int(v)
                except ValueError:
                    temp = v.split(' ')[0]
                    try:
                        self.metadata[k] = int(temp)
                    except ValueError:
                        self.metadata[k] = int(temp.split('.')[0])
            elif k in ('RangeResolutionKMeters', 'CTF', 'TransmitCenterFreqMHz', 'DopplerResolutionHzPerBin',
                       'RadialBraggPeakDropOff', 'RadialBraggPeakNull', 'RadialBraggNoiseThreshold', 'TransmitSweepRateHz',
                       'TransmitBandwidthKHz'):
                try:
                    self.metadata[k] = float(v)
                except ValueError:
                    self.metadata[k] = float(v.split(' ')[0])
            else:
                continue

        required = ['Origin', 'TransmitCenterFreqMHz']
        present_keys = self.metadata.keys()
        for key in required:
            if key not in present_keys:
                self.metadata[key] = None

    def create_netcdf(self, filename):
        """
        Create a compressed netCDF4 (.nc) file from the radial instance
        :param filename: User defined filename of radial file you want to save
        :return:
        """
        create_dir(os.path.dirname(filename))
        rename = dict(VELU='u',
                      VELV='v',
                      VFLG='vector_flag',
                      ESPC='spatial_quality',
                      ETMP='temporal_quality',
                      MAXV='velocity_max',
                      MINV='velocity_min',
                      ERSC='spatial_count',
                      ERTC='temporal_count',
                      XDST='dist_east_from_origin',
                      YDST='dist_north_from_origin',
                      VELO='velocity',
                      HEAD='heading',
                      SPRC='range_cell')

        # rename variables to something meaningful
        self.data.rename(rename, inplace=True)

        # set time attribute
        self.data['time'].attrs['standard_name'] = 'time'

        # Set lon attributes
        self.data['lon'].attrs['long_name'] = 'Longitude'
        self.data['lon'].attrs['standard_name'] = 'longitude'
        self.data['lon'].attrs['short_name'] = 'lon'
        self.data['lon'].attrs['units'] = 'degrees_east'
        self.data['lon'].attrs['axis'] = 'X'
        self.data['lon'].attrs['valid_min'] = np.float32(-180.0)
        self.data['lon'].attrs['valid_max'] = np.float32(180.0)

        # Set lat attributes
        self.data['lat'].attrs['long_name'] = 'Latitude'
        self.data['lat'].attrs['standard_name'] = 'latitude'
        self.data['lat'].attrs['short_name'] = 'lat'
        self.data['lat'].attrs['units'] = 'degrees_north'
        self.data['lat'].attrs['axis'] = 'Y'
        self.data['lat'].attrs['valid_min'] = np.float32(-90.0)
        self.data['lat'].attrs['valid_max'] = np.float32(90.0)

        # Set u attributes
        self.data['u'].attrs['long_name'] = 'Eastward Surface Current (cm/s)'
        self.data['u'].attrs['standard_name'] = 'surface_eastward_sea_water_velocity'
        self.data['u'].attrs['short_name'] = 'u'
        self.data['u'].attrs['units'] = 'cm s-1'
        self.data['u'].attrs['valid_min'] = np.float32(-300)
        self.data['u'].attrs['valid_max'] = np.float32(300)
        self.data['u'].attrs['coordinates'] = 'lon lat'
        self.data['u'].attrs['grid_mapping'] = 'crs'

        # Set v attributes
        self.data['v'].attrs['long_name'] = 'Northward Surface Current (cm/s)'
        self.data['v'].attrs['standard_name'] = 'surface_northward_sea_water_velocity'
        self.data['v'].attrs['short_name'] = 'v'
        self.data['v'].attrs['units'] = 'cm s-1'
        self.data['v'].attrs['valid_min'] = np.float32(-300)
        self.data['v'].attrs['valid_max'] = np.float32(300)
        self.data['v'].attrs['coordinates'] = 'lon lat'
        self.data['v'].attrs['grid_mapping'] = 'crs'

        # Set bearing attributes
        self.data['bearing'].attrs['long_name'] = 'Bearing from origin (away from instrument)'
        self.data['bearing'].attrs['short_name'] = 'bearing'
        self.data['bearing'].attrs['units'] = 'degrees'
        self.data['bearing'].attrs['valid_min'] = np.float32(0)
        self.data['bearing'].attrs['valid_max'] = np.float32(360)
        self.data['bearing'].attrs['grid_mapping'] = 'crs'
        self.data['bearing'].attrs['axis'] = 'Y'

        # Set range attributes
        self.data['range'].attrs['long_name'] = 'Range from origin (away from instrument)'
        self.data['range'].attrs['short_name'] = 'range'
        self.data['range'].attrs['units'] = 'km'
        self.data['range'].attrs['valid_min'] = np.float32(0)
        self.data['range'].attrs['valid_max'] = np.float32(1000)
        self.data['range'].attrs['grid_mapping'] = 'crs'
        self.data['range'].attrs['axis'] = 'X'

        # vector_flag
        self.data['vector_flag'].attrs['long_name'] = 'Vector Flag Masks'
        self.data['vector_flag'].attrs['valid_range'] = [0, 2048]
        self.data['vector_flag'].attrs['flag_masks'] = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
        self.data['vector_flag'].attrs['flag_meanings'] = 'grid_point_deleted grid_point_near_coast point_measurement no_radial_solution baseline_interpolation exceeds_max_speed invalid_solution solution_beyond_valid_spatial_domain insufficient_angular_resolution reserved reserved'
        self.data['vector_flag'].attrs['coordinates'] = 'lon lat'
        self.data['vector_flag'].attrs['grid_mapping'] = 'crs'

        # spatial_quality
        self.data['spatial_quality'].attrs['long_name'] = 'Spatial Quality of radial sea water velocity'
        self.data['spatial_quality'].attrs['units'] = 'cm s-1'
        self.data['spatial_quality'].attrs['coordinates'] = 'lon lat'
        self.data['spatial_quality'].attrs['grid_mapping'] = 'crs'

        # temporal_quality
        self.data['temporal_quality'].attrs['long_name'] = 'Temporal Quality of radial sea water velocity'
        self.data['temporal_quality'].attrs['units'] = 'cm s-1'
        self.data['temporal_quality'].attrs['coordinates'] = 'lon lat'
        self.data['temporal_quality'].attrs['grid_mapping'] = 'crs'

        # velocity_max
        self.data['velocity_max'].attrs['long_name'] = 'Maximum Velocity of sea water (away from instrument)'
        self.data['velocity_max'].attrs['units'] = 'cm s-1'
        self.data['velocity_max'].attrs['coordinates'] = 'lon lat'
        self.data['velocity_max'].attrs['grid_mapping'] = 'crs'

        # velocity_min
        self.data['velocity_min'].attrs['long_name'] = 'Minimum Velocity of sea water (away from instrument)'
        self.data['velocity_min'].attrs['units'] = 'cm s-1'
        self.data['velocity_min'].attrs['coordinates'] = 'lon lat'
        self.data['velocity_min'].attrs['grid_mapping'] = 'crs'

        # spatial_count
        self.data['spatial_count'].attrs['long_name'] = 'Spatial count of sea water velocity (away from instrument)'
        self.data['spatial_count'].attrs['coordinates'] = 'lon lat'
        self.data['spatial_count'].attrs['grid_mapping'] = 'crs'

        # temporal_count
        self.data['temporal_count'].attrs['long_name'] = 'Temporal count of sea water velocity (away from instrument)'
        self.data['temporal_count'].attrs['coordinates'] = 'lon lat'
        self.data['temporal_count'].attrs['grid_mapping'] = 'crs'

        # east_dist_from_origin
        self.data['dist_east_from_origin'].attrs['long_name'] = 'Eastward distance from instrument'
        self.data['dist_east_from_origin'].attrs['units'] = 'km'
        self.data['dist_east_from_origin'].attrs['coordinates'] = 'lon lat'
        self.data['dist_east_from_origin'].attrs['grid_mapping'] = 'crs'

        # north_dist_from_origin
        self.data['dist_north_from_origin'].attrs['long_name'] = 'Northward distance from instrument'
        self.data['dist_north_from_origin'].attrs['units'] = 'km'
        self.data['dist_north_from_origin'].attrs['coordinates'] = 'lon lat'
        self.data['dist_north_from_origin'].attrs['grid_mapping'] = 'crs'

        # velocity
        self.data['velocity'].attrs['valid_range'] = [-1000, 1000]
        self.data['velocity'].attrs['standard_name'] = 'radial_sea_water_velocity_away_from_instrument'
        self.data['velocity'].attrs['units'] = 'cm s-1'
        self.data['velocity'].attrs['coordinates'] = 'lon lat'
        self.data['velocity'].attrs['grid_mapping'] = 'crs'

        # heading
        self.data['heading'].attrs['valid_range'] = [0, 3600]
        self.data['heading'].attrs['standard_name'] = 'direction_of_radial_vector_away_from_instrument'
        self.data['heading'].attrs['units'] = 'degrees'
        self.data['heading'].attrs['coordinates'] = 'lon lat'
        self.data['heading'].attrs['scale_factor'] = 0.1
        self.data['heading'].attrs['grid_mapping'] = 'crs'

        # range_cell
        self.data['range_cell'].attrs['long_name'] = 'Cross Spectra Range Cell  of sea water velocity (away from instrument)'
        self.data['range_cell'].attrs['coordinates'] = 'lon lat'
        self.data['range_cell'].attrs['grid_mapping'] = 'crs'

        del self.data.attrs['TimeStamp']
        encoding = make_encoding(self.data, comp_level=4, fillvalue=np.nan)
        encoding['bearing'] = dict(zlib=False, _FillValue=False)
        encoding['range'] = dict(zlib=False, _FillValue=False)
        # encoding['z'] = dict(zlib=False, _FillValue=False)
        self.data.to_netcdf(filename, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims=['time'])

    def create_ruv(self, filename):
        """
        Create a CODAR Radial (.ruv) file from radial instance
        :param filename: User defined filename of radial file you want to save
        :return:
        """
        create_dir(os.path.dirname(filename))

        with open(filename, 'w') as f:
            # Write header
            for metadata_key, metadata_value in self.metadata.items():
                if 'ProcessedTimeStamp' in metadata_key:
                    break
                else:
                    f.write('%{}: {}\n'.format(metadata_key, metadata_value))

            # Write data tables. Anything beyond the first table is commented out.
            for table in self._tables.keys():
                for table_key, table_value in self._tables[table].items():
                    if table_key is not 'data':
                        if (table_key == 'TableType') & (table == '1'):
                            if 'QCTest' in self.metadata:
                                f.write('%QCReference: Quality control reference: IOOS QARTOD HF Radar ver 1.0 May 2016\n')
                                f.write('%QCFlagDefinitions: 1=pass 2=not_evaluated 3=suspect 4=fail 9=missing_data\n')
                                f.write('%QCTestFormat: "test_name [qc_thresholds]: test_result"\n')

                                for test in self.metadata['QCTest']:
                                    f.write('%QCTest: {}\n'.format(test))
                            f.write('%{}: {}\n'.format(table_key, table_value))
                        elif table_key == 'TableColumns':
                            f.write('%TableColumns: {}\n'.format(len(self._tables[table]['data'].columns)))
                        elif table_key == 'TableStart':
                            f.write('%{}: {}\n'.format(table_key, table_value))
                            for line in self._tables[table]['_TableHeader']:
                                f.write(line)
                        elif table_key == '_TableHeader':
                            pass
                        else:
                            f.write('%{}: {}\n'.format(table_key, table_value))

                if 'datetime' in self._tables[table]['data'].keys():
                    self._tables[table]['data'] = self._tables[table]['data'].drop(['datetime'], axis=1)

                if table == '1':
                    # f.write('%{}\n'.format(self._tables[table]['TableColumnTypes']))
                    # Fill NaN with 999.000 which is the standard fill value for codar lluv filesself._tables[table]['TableColumnTypes']
                    self.data = self.data.fillna(999.000)
                    self.data.to_string(f, index=False, justify='center', header=False)
                else:
                    f.write('%%{}\n'.format(self._tables[table]['TableColumnTypes']))
                    self._tables[table]['data'].insert(0, '%%', '%')
                    self._tables[table]['data'] = self._tables[table]['data'].fillna(999.000)
                    self._tables[table]['data'].to_string(f, index=False, justify='center', header=False)

                if int(table) > 1:
                    f.write('\n%TableEnd: {}\n'.format(table))
                else:
                    f.write('\n%TableEnd: \n')
                f.write('%%\n')

            # Write footer containing processing information
            f.write('%ProcessedTimeStamp: {}\n'.format(self.metadata['ProcessedTimeStamp']))
            for tool in self.metadata['ProcessingTool']:
                f.write('%ProcessingTool: {}\n'.format(tool))
                # f.write('%{}: {}\n'.format(footer_key, footer_value))
            f.write('%End:')

    def export(self, filename, file_type='radial'):
        """
        Export radial file as either a codar .ruv file or a netcdf .nc file
        :param filename: User defined filename of radial file you want to save
        :param file_type: Type of file to export radial: radial (default) or netcdf
        :return:
        """

        if not self.is_valid():
            raise ValueError("Could not export ASCII data, the input file was invalid.")

        if os.path.isfile(filename):
            os.remove(filename)

        if file_type == 'radial':
            self.create_ruv(filename)
        elif file_type == 'netcdf':
            self.create_netcdf(filename)

    def initialize_qc(self):
        self.metadata['QCTest'] = []

    # QARTOD QC Tests
    def qc_qartod_avg_radial_bearing(self, reference_bearing, warning_threshold=15, failure_threshold=30):
        """
        Integrated Ocean Observing System (IOOS) Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Valid Location (Test 12)
        Check that the average radial bearing remains relatively constant (Roarty et al. 2012).

        It is expected that the average of all radial velocity bearings AVG_RAD_BEAR obtained during a sample
        interval (e.g., 1 hour) should be close to a reference bearing REF_RAD_BEAR and not vary beyond warning
        or failure thresholds.
        :return:
        """
        # Absolute value of the difference between the bearing mean and reference bearing
        absolute_difference = np.abs(self.data['BEAR'].mean() - reference_bearing)

        if absolute_difference >= failure_threshold:
            flag = 4
        elif absolute_difference >= warning_threshold & absolute_difference < failure_threshold:
            flag = 3
        elif absolute_difference < warning_threshold:
            flag = 1
        self.metadata['QCTest'].append(f'"qc_qartod_avg_radial_bearing (QC12) [reference_bearing={str(reference_bearing)} (degrees) warning={str(warning_threshold)} (degrees) failure={str(failure_threshold)} (degrees)]: {str(flag)}"')

    def qc_qartod_valid_location(self):
        """
        Integrated Ocean Observing System (IOOS) Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Valid Location (Test 8)
        Removes radial vectors placed over land or in other unmeasureable areas

        Radial vector coordinates are checked against a reference file containing information about which locations
        are over land or in an unmeasurable area (for example, behind an island or point of land). Radials in these
        areas will be flagged with a code (FLOC) in the radial file (+128 in CODAR radial files) and are not included
        in total vector calculations.

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/
        :return:
        """
        self.data['QC08'] = 1
        boolean = self.data['VFLG'] == 128

        self.data['QC08'] = self.data['QC08'].where(~boolean, other=4)
        self._tables['1']['TableColumnTypes'] += ' QC08'
        self.metadata['QCTest'].append('"qc_qartod_valid_location (QC08)[VFLG==128]: See results in column QC08 below"')

    def qc_qartod_radial_count(self, min_radials=150, low_radials=300):
        """
        Integrated Ocean Observing System (IOOS) Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Radial Count (Test 9)
        Rejects radials in files with low radial counts (poor radial map coverage).

        The number of radials (RCNT) in a radial file must be above a threshold value RCNT_MIN to pass the test and
        above a value RC_LOW to not be considered suspect. If the number of radials is below the minimum level,
        it indicates a problem with data collection. In this case, the file should be rejected and none of the
        radials used for total vector processing.

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/
        :param min_radials: Minimum radial count threshold below which the file should be rejected. min_radials < low_radials
        :param low_radials: Low radial count threshold below which the file should be considered suspect. low_radials > min_radials
        :return:
        """
        num_radials = self.data.shape[0]
        if num_radials < min_radials:
            radial_count_flag = 4
        elif (num_radials >= min_radials) and (num_radials <= low_radials):
            radial_count_flag = 3
        elif num_radials > low_radials:
            radial_count_flag = 1

        # self.metadata['qc_qartod_radial_count'] = str(radial_count_flag)
        self.metadata['QCTest'].append(f'"qc_qartod_radial_count (QC09) [failure={str(min_radials)} (number of radials) warning_num={str(low_radials)}]: {str(radial_count_flag)} (number of radials)]"')

    def qc_qartod_maximum_velocity(self, threshold=250):
        """
        Integrated Ocean Observing System (IOOS) Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Max Threshold (Test 7)
        Ensures that a radial current speed is not unrealistically high.

        The maximum radial speed threshold (RSPDMAX) represents the maximum reasonable surface radial velocity
        for the given domain.

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/
        :param threshold: Maximum Radial Speed (cm/s)
        :return:
        """
        self.data['QC07'] = 1
        try:
            boolean = self.data['VELO'].abs() > threshold
        except TypeError:
            self.data['VELO'] = self.data['VELO'].astype(float)
            boolean = self.data['VELO'].abs() > threshold

        self.data['QC07'] = self.data['QC07'].where(~boolean, other=4)
        self._tables['1']['TableColumnTypes'] += ' QC07'
        self.metadata['QCTest'].append(f'"qc_qartod_maximum_velocity (QC07) [max_vel={str(threshold)} (cm/s)]: See results in column QC07 below"')

    def qc_qartod_spatial_median(self, rclim=2.1, anglim=10, threshold=30):
        """
        Integrated Ocean Observing System (IOOS) Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Spatial Median (Test 10)
        Ensures that the radial velocity is not too different from nearby radial velocities.
        RV is the radial velocity
        NV is a set of radial velocities for neighboring radial cells (cells within radius of 'rclim' * Range Step (km)
        and whose vector bearing (angle of arrival at site) is also within 'anglim' degrees of the source vector's bearing)
        Required to pass the test: |RV - median(NV)| <= threshold
        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/
        :param RCLim: multiple of range step which depends on the radar type
        :param AngLim: limit for number of degrees from source radial's bearing (degrees)
        :param CurLim: Current difference threshold (cm/s)
        :return:
        """
        self.data['QC10'] = 1
        try:
            Rstep = float(self.metadata['RangeResolutionKMeters'])
            # Rstep = np.floor(min(np.diff(np.unique(self.data['RNGE'])))) #use as backup method if other fails?

            Bstep = [float(s) for s in re.findall(r'-?\d+\.?\d*', self.metadata['AngularResolution'])]
            Bstep = Bstep[0]
            # Bstep = int(min(np.diff(np.unique(self.data['BEAR']))))  #use as backup method if other fails?

            RLim = int(round(rclim))  # if not an integer will cause an error later on
            BLim = int(anglim / Bstep)  # if not an integer will cause an error later on

            # convert bearing into bearing cell numbers
            adj = int(5 - min(self.data['BEAR']))
            Bcell = ((self.data['BEAR'] + adj) / Bstep) - 1
            Bcell = Bcell.astype(int)
            # Btable = np.column_stack((self.data['BEAR'], Bcell))  #only for debugging

            # convert range into range cell numbers
            Rcell = (np.floor((self.data['RNGE'] / Rstep) + 0.1))
            Rcell = Rcell - min(Rcell)
            Rcell = Rcell.astype(int)
            # Rtable = np.column_stack((self.data['RNGE'], Rcell))   #only for debugging
            Rcell = self.data['SPRC']

            # place velocities into a matrix with rows defined as bearing cell# and columns as range cell#
            BRvel = np.zeros((int(360 / Bstep), max(Rcell) + 1), dtype=int) + np.nan
            BRind = np.zeros((int(360 / Bstep), max(Rcell) + 1), dtype=int) + np.nan

            for xx in range(len(self.data['VELO'])):
                BRvel[Bcell[xx]][Rcell[xx]] = self.data['VELO'][xx]
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
                    temp = BRpad[bb - BLim:bb + BLim + 1, rr - RLim:rr + RLim + 1]  # temp is the matrix of neighbors
                    BRmed[bb][rr] = np.nanmedian(temp)

            # now remove the padding from the array containing the median values
            # I should be able to index the cells I want, but I could not
            # figure out how to do this so I did a two step deletion process instead!
            BRmedtrim = np.delete(BRmed, [[0, BLim - 1, 1], [BRmed.shape[0] - BLim, BRmed.shape[0] - 1, 1]], axis=0)
            BRmedtrim = np.delete(BRmedtrim, [[0, RLim - 1, 1], [BRmed.shape[1] - RLim, BRmed.shape[1] - 1, 1]], axis=1)

            # calculate velocity minus median of neighbors
            # and put back into single column using the indices saved in BRind
            BRdiff = BRvel - BRmedtrim  # velocity minus median of neighbors, test these values against current threshold
            diffcol = self.data['RNGE'] + np.nan  # initialize a single column for the difference results
            for rr in range(BRdiff.shape[1]):
                for bb in range(BRdiff.shape[0]):
                    if not (np.isnan(BRind[bb][rr])):
                        diffcol[BRind[bb][rr]] = BRdiff[bb][rr]
            boolean = diffcol.abs() > threshold

            # Another method would take the median of data from any radial cells within a certain
            # distance (radius) of the radial being tested.  This method, as coded below, was very slow!
            # Perhaps there is a better way to write the code.
            # dist contains distance between each location and the other locations
            # for rvi in range(len(self.data['VELO'])):
            #     dist = np.zeros((len(self.data['LATD']), 1))
            #     for i in range(len(self.data['LATD'])):
            #         rvloc = self.data['LATD'][i],self.data['LOND'][i]
            #         dist[i][0] = distance.distance((self.data['LATD'][rvi],self.data['LOND'][rvi]),(self.data['LATD'][i],self.data['LOND'][i])).kilometers

        except TypeError:
            diffcol = diffcol.astype(float)
            boolean = diffcol.abs() > threshold

        self.data['QC10'] = self.data['QC10'].where(~boolean, other=4)
        # self.data['VFLG'] = self.data['VFLG'].where(~boolean, other=4) # for testing only, shows up as "marker" flag in SeaDisplay
        self._tables['1']['TableColumnTypes'] += ' QC10'
        self.metadata['QCTest'].append(f'"qc_qartod_spatial_median (QC10) [range_cell_limit={str(rclim)} (range cells) angular_limit={str(anglim)} (degrees) current_difference={str(threshold)} (cm/s)]: See results in column QC10 below"')

    def qc_qartod_syntax(self):
        """
        Integrated Ocean Observing System (IOOS) Quality Assurance of Real-Time Oceanographic Data (QARTOD)
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
        :return:
        """
        i = 0

        # check for timestamp in filename
        result = re.search('\d{4}_\d{2}_\d{2}_\d{4}', self.file_name)
        if result:
            timestr = dt.datetime.strptime(result.group(), '%Y_%m_%d_%H%M')
            i = i + 1

        # Radial tables must not be empty
        if self.is_valid():
            i = i + 1

        # The following metadata must be defined.
        if self.metadata['FileType'] and self.metadata['Site'] and self.metadata['TimeStamp'] and self.metadata['Origin'] and self.metadata['PatternType'] and self.metadata['TimeZone']:
            filetime = dt.datetime(*map(int, self.metadata['TimeStamp'].split()))
            i = i + 1

        # Filename timestamp must match the timestamp reported within the file.
        if timestr == filetime:
            i = i + 1

        # Radial data table columns stated must match the number of columns reported for each row
        if len(self._tables['1']['TableColumnTypes'].split()) == self.data.shape[1]:
            i = i + 1

        # Make sure site location is within range: -180 <= lon <= 180 & -90 <= lat <= 90
        latlon = re.findall(r"[-+]?\d*\.\d+|\d+", self.metadata['Origin'])
        if (-180 <= float(latlon[1]) <= 180) & (-90 <= float(latlon[0]) <= 90):
            i = i + 1

        if i == 6:
            syntax = 1
        else:
            syntax = 4
        self.metadata['QCTest'].append(f'"qc_qartod_syntax (QC06) [N/A]: {syntax}"')

    def reset(self):
        logging.info('Resetting instance data variable to original dataset')
        self._tables['1']
        self.data = self._data_backup

