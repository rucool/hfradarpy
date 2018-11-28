import datetime as dt
import geopandas as gpd
import logging
import numpy as np
import os
import pandas as pd
import re
from shapely.geometry import Point
import xarray as xr
from codar_processing.common import LLUVParser, create_dir
from codar_processing.calc import reckon, gridded_index


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
        radial_dict[radial.file_name] = radial.ds

    ds = xr.concat(radial_dict.values(), 'time')
    return ds


class Radial(LLUVParser):
    """
    Radial Subclass.

    This class should be used when loading a CODAR radial (.ruv) file. This class utilizes the generic LLUV class from
    ~/codar_processing/common.py in order to load CODAR Radial files
    """
    def __init__(self, fname, replace_invalid=True, multi_dimensional=False, mask_over_land=False):
        keep = ['LOND', 'LATD', 'VELU', 'VELV', 'VFLG', 'ESPC', 'ETMP', 'MAXV', 'MINV', 'ERSC', 'ERTC', 'XDST', 'YDST', 'RNGE', 'BEAR', 'VELO', 'HEAD', 'SPRC']

        LLUVParser.__init__(self, fname)
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

    def to_multi_dimensional(self, range_min=0, range_max=150):
        rename = dict(VELU='u',
                      VELV='v',
                      VFLG='vector_indicator_flag',
                      ESPC='spatial_standard_deviation',
                      ETMP='temporal_standard_deviation',
                      MAXV='maximum_current_velocity',
                      MINV='minimum_current_velocity',
                      ERSC='spatial_num_velocities',
                      ERTC='temporal_num_velocities',
                      XDST='east_dist_from_origin',
                      YDST='north_dist_from_origin',
                      RNGE='range_old',
                      BEAR='bearing_old',
                      VELO='velocity',
                      HEAD='heading',
                      SPRC='range_cell')
        coords = ('time', 'range', 'bearing')

        # Intitialize empty xarray dataset
        ds = xr.Dataset()
        range_dim = np.arange(range_min, range_max + float(self.metadata['RangeResolutionKMeters']), float(self.metadata['RangeResolutionKMeters']))

        # Clean radial header
        self.clean_header()

        # bearing_dim - 360 degrees
        bearing_dim = np.arange(0, 360 + float(self.metadata['AngularResolution']), self.metadata['AngularResolution'])

        # create radial grid from bearing and range
        [bearing, ranges] = np.meshgrid(bearing_dim, range_dim)

        # calculate lat/lons from origin, bearing, and ranges
        latlon = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", self.metadata['Origin'])]
        new_lat, new_lon = reckon(latlon[0], latlon[1], bearing, ranges)

        # create dictionary containing variables from dataframe in the shape of radial grid
        d = {key: np.tile(np.nan, bearing.shape) for key in self.data.keys()}

        # find grid indices from radial grid (bearing, ranges)
        x_ind, y_ind = gridded_index(bearing, ranges, self.data['BEAR'].values, self.data['RNGE'].values)

        for k, v in d.items():
            v[(y_ind, x_ind)] = self.data[k]
            d[k] = v

        # Add extra dimension for time
        d = {k: np.expand_dims(np.float32(v), axis=0) for (k, v) in d.items()}

        ds['lon'] = (('range', 'bearing'), new_lon.round(4))
        ds['lat'] = (('range', 'bearing'), new_lat.round(4))

        # Add all variables to dataset
        for k, v in d.items():
            ds[k] = (coords, v)

        # Add coordinate variables to dataset
        ds.coords['time'] = pd.date_range(self.metadata['TimeStamp'], periods=1)
        ds.coords['bearing'] = bearing_dim
        ds.coords['range'] = range_dim

        # rename variables to something meaningful
        ds.rename(rename, inplace=True)

        # Assign header data to global attributes
        self.data = ds.assign_attrs(self.metadata)

    def file_type(self):
        """Return a string representing the type of file this is."""
        return 'radial'

    def clean_header(self):
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

    # def create_nc(self, filename):
    #     """
    #     Create a compressed netCDF4 (.nc) file from the radial instance
    #     :param filename: User defined filename of radial file you want to save
    #     :return:
    #     """

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
                        f.write('%{}: {}\n'.format(table_key, table_value))

                if 'datetime' in self._tables[table]['data'].keys():
                    self._tables[table]['data'] = self._tables[table]['data'].drop(['datetime'], axis=1)

                if table == '1':
                    f.write('%{}\n'.format(self._tables[table]['TableColumnTypes']))
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
        # elif file_type == 'netcdf':
        #     self.create_netcdf(filename)

    # QARTOD QC Tests
    def qc_qartod_location(self):
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
        self.data['VLOC'] = 1
        boolean = self.data['VFLG'] == 128

        self.data['VLOC'] = self.data['VLOC'].where(~boolean, other=4)
        self._tables['1']['TableColumnTypes'] += ' VLOC'
        self.metadata['ProcessingTool'].append('"hfr_processing/Radial.qc_qartod_location"')

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

        self.metadata['qc_qartod_radial_count'] = str(radial_count_flag)
        self.metadata['ProcessingTool'].append('"hfr_processing/Radial.qc_qartod_radial_count"')

    def qc_qartod_speed(self, threshold=250):
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
        self.data['MVEL'] = 1
        try:
            boolean = self.data['VELO'].abs() > threshold
        except TypeError:
            self.data['VELO'] = self.data['VELO'].astype(float)
            boolean = self.data['VELO'].abs() > threshold

        self.data['MVEL'] = self.data['MVEL'].where(~boolean, other=4)
        self._tables['1']['TableColumnTypes'] += ' MVEL'
        self.metadata['ProcessingTool'].append('"hfr_processing/Radial.qc_qartod_speed"')

    def qc_qartod_spatialmedian(self, rclim=2.1, anglim=10, threshold = 30):
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
            self.data['SMED'] = 1
            try:
                Rstep = float(self.header['RangeResolutionKMeters'])
                #Rstep = np.floor(min(np.diff(np.unique(self.data['RNGE'])))) #use as backup method if other fails?

                Bstep = [float(s) for s in re.findall(r'-?\d+\.?\d*', self.header['AngularResolution'])]
                Bstep = Bstep[0]
                #Bstep = int(min(np.diff(np.unique(self.data['BEAR']))))  #use as backup method if other fails?

                RLim = int(round(rclim)) #if not an integer will cause an error later on
                BLim = int(anglim/Bstep) #if not an integer will cause an error later on

                #convert bearing into bearing cell numbers
                adj = int(5 - min(self.data['BEAR']))
                Bcell = ((self.data['BEAR'] + adj) / Bstep) - 1
                Bcell = Bcell.astype(int)
                #Btable = np.column_stack((self.data['BEAR'], Bcell))  #only for debugging

                #convert range into range cell numbers
                Rcell = (np.floor((self.data['RNGE'] / Rstep) + 0.1))
                Rcell = Rcell - min(Rcell)
                Rcell = Rcell.astype(int)
                #Rtable = np.column_stack((self.data['RNGE'], Rcell))   #only for debugging
                Rcell = self.data['SPRC']

                # place velocities into a matrix with rows defined as bearing cell# and columns as range cell#
                BRvel = np.zeros((int(360 / Bstep), max(Rcell)+1), dtype=int) + np.nan
                BRind = np.zeros((int(360 / Bstep), max(Rcell)+1), dtype=int) + np.nan

                for xx in range(len(self.data['VELO'])):
                   BRvel[Bcell[xx]][Rcell[xx]] = self.data['VELO'][xx]
                   BRind[Bcell[xx]][Rcell[xx]] = xx   #keep track of indices so easier to return to original format

                # deal with 359 to 0 transition in bearing by
                # repeating first BLim rows at the bottom and last BLim rows at the top
                # also pad ranges with NaNs by adding extra columns on the left and right of the array
                # this keeps the indexing for selecting the neighbors from breaking

                BRtemp = np.append(np.append(BRvel[-BLim:], BRvel, axis=0), BRvel[:BLim],axis=0)
                rangepad = np.zeros((BRtemp.shape[0], RLim), dtype=int) + np.nan
                BRpad = np.append(np.append(rangepad, BRtemp,axis=1), rangepad,axis=1)

                #calculate median of neighbors (neighbors include the point itself)
                BRmed = BRpad + np.nan  # initialize with an array of NaN
                for rr in range(RLim,BRvel.shape[1]+RLim):
                    for bb in range(BLim,BRvel.shape[0]+BLim):
                        temp = BRpad[bb-BLim:bb+BLim+1,rr-RLim:rr+RLim+1]  #temp is the matrix of neighbors
                        BRmed[bb][rr] = np.nanmedian(temp)

                # now remove the padding from the array containing the median values
                # I should be able to index the cells I want, but I could not
                # figure out how to do this so I did a two step deletion process instead!
                BRmedtrim = np.delete(BRmed, [[0, BLim - 1, 1], [BRmed.shape[0] - BLim, BRmed.shape[0] - 1, 1]], axis=0)
                BRmedtrim = np.delete(BRmedtrim, [[0, RLim - 1, 1], [BRmed.shape[1] - RLim, BRmed.shape[1] - 1, 1]], axis=1)

                # calculate velocity minus median of neighbors
                # and put back into single column using the indices saved in BRind
                BRdiff = BRvel - BRmedtrim # velocity minus median of neighbors, test these values against current threshold
                diffcol = self.data['RNGE'] + np.nan  #initialize a single column for the difference results
                for rr in range(BRdiff.shape[1]):
                    for bb in range(BRdiff.shape[0]):
                        if not(np.isnan(BRind[bb][rr])):
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


            self.data['SMED'] = self.data['SMED'].where(~boolean, other=4)
            #self.data['VFLG'] = self.data['VFLG'].where(~boolean, other=4) # for testing only, shows up as "marker" flag in SeaDisplay
            self.tables['1']['TableColumnTypes'] += ' SMED'
            self.footer['ProcessingTool'].append('"hfr_processing/Radial.qc_qartod_spatialmedian"')


    def reset(self):
        logging.info('Resetting instance data variable to original dataset')
        self._tables['1']
        self.data = self._data_backup
