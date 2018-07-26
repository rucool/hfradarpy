import datetime as dt
import logging
import os
from codar_processing.common import LLUVParser, create_dir

logger = logging.getLogger(__name__)


class Radial(LLUVParser):
    """
    Radial Subclass.

    This class should be used when loading a CODAR radial (.ruv) file. This class utilizes the generic LLUV class from
    ~/codar_processing/common.py in order to load CODAR Radial files
    """

    def __init__(self, fname):
        LLUVParser.__init__(self, fname)

        for key in self.tables.keys():
            if self.is_valid(key):
                table = self.tables[key]
                if 'LLUV' in table['TableType']:
                    table['data'].insert(0, '%%', '')
                    self.data = table['data']
                elif 'rads' in table['TableType']:
                    table['data'].insert(0, '%%', '%')
                    self.diags_radial = table['data']
                    self.diags_radial['datetime'] = self.diags_radial[['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)
                elif 'rcvr' in table['TableType']:
                    table['data'].insert(0, '%%', '%')
                    self.diags_hardware = table['data']
                    self.diags_hardware['datetime'] = self.diags_hardware[['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)
        # print('test')

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
        key_list = list(self.header.keys())
        for key in key_list:
            if not key in keep:
                del self.header[key]

        for k, v in self.header.items():
            if 'Site' in k:
                self.header[k] = ''.join(e for e in v if e.isalnum())
            elif k in ('TimeStamp', 'PatternDate'):
                t_list = [int(s) for s in v.split()]
                self.header[k] = dt.datetime(*t_list)
            elif 'TimeZone' in k:
                self.header[k] = v.split('"')[1]
            elif 'TableColumnTypes' in k:
                self.header[k] = ' '.join([x.strip() for x in v.strip().split(' ')])
            elif 'Origin' in k:
                self.header[k] = v.lstrip()
            elif k in ('RangeStart', 'RangeEnd', 'AntennaBearing', 'ReferenceBearing', 'AngularResolution', 'SpatialResolution',
                       'FirstOrderMethod', 'BraggSmoothingPoints', 'BraggHasSecondOrder', 'MergedCount',
                       'RadialMinimumMergePoints', 'FirstOrderCalc', 'SpectraRangeCells', 'SpectraDopplerCells',
                       'TableColumns', 'TableRows',  'PatternResolution', 'CurrentVelocityLimit', 'TimeCoverage'):
                try:
                    self.header[k] = int(v)
                except ValueError:
                    temp = v.split(' ')[0]
                    try:
                        self.header[k] = int(temp)
                    except ValueError:
                        self.header[k] = int(temp.split('.')[0])
            elif k in ('RangeResolutionKMeters', 'CTF', 'TransmitCenterFreqMHz', 'DopplerResolutionHzPerBin',
                       'RadialBraggPeakDropOff', 'RadialBraggPeakNull', 'RadialBraggNoiseThreshold', 'TransmitSweepRateHz',
                       'TransmitBandwidthKHz'):
                try:
                    self.header[k] = float(v)
                except ValueError:
                    self.header[k] = float(v.split(' ')[0])
            else:
                continue

        required = ['Origin', 'TransmitCenterFreqMHz']
        present_keys = self.header.keys()
        for key in required:
            if key not in present_keys:
                self.header[key] = None

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
            for header_key, header_value in self.header.items():
                f.write('%{}: {}\n'.format(header_key, header_value))

            # Write data tables. Anything beyond the first table is commented out.
            for table in self.tables.keys():
                for table_key, table_value in self.tables[table].items():
                    if table_key is not 'data':
                        f.write('%{}: {}\n'.format(table_key, table_value))

                if 'datetime' in self.tables[table]['data'].keys():
                    self.tables[table]['data'] = self.tables[table]['data'].drop(['datetime'], axis=1)

                # Fill NaN with 999.000 which is the standard fill value for codar lluv files
                if table == '1':
                    self.data = self.data.fillna(999.000)
                    self.data.to_string(f, index=False, justify='center')
                else:
                    self.tables[table]['data'] = self.tables[table]['data'].fillna(999.000)
                    self.tables[table]['data'].to_string(f, index=False, justify='center')

                if int(table) > 1:
                    f.write('\n%TableEnd: {}\n'.format(table))
                else:
                    f.write('\n%TableEnd: \n')
                f.write('%%\n')

            # Write footer containing processing information
            for footer_key, footer_value in self.footer.items():
                if footer_key == 'ProcessingTool':
                    for tool in self.footer['ProcessingTool']:
                        f.write('%ProcessingTool: {}\n'.format(tool))
                else:
                    f.write('%{}: {}\n'.format(footer_key, footer_value))

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
        self.tables['1']['TableColumnTypes'] += ' VLOC'
        self.footer['ProcessingTool'].append('"hfr_processing/Radial.qc_qartod_location"')

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

        self.header['qc_qartod_radial_count'] = str(radial_count_flag)
        self.footer['ProcessingTool'].append('"hfr_processing/Radial.qc_qartod_radial_count"')

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
        self.tables['1']['TableColumnTypes'] += ' MVEL'
        self.footer['ProcessingTool'].append('"hfr_processing/Radial.qc_qartod_speed"')

    def reset(self):
        logging.info('Resetting instance data variable to original dataset')
        self.tables['1']
        self.data = self._data_backup
    # # Modify any tableheader information that needs to be updated
    # radial_file_data['tables']['1']['TableColumns'] = radial_data.shape[1]
    # header_list = radial_data.columns.tolist()
    # header_list.remove('%%')
    # radial_file_data['tables']['1']['TableColumnTypes'] = ' '.join(header_list)

    # def qc_average_radial_bearing(df):
    #     return df

    # def qc_spatial_median_filter():
    # """
    # """
    # return def

    # def qc_temporal_gradient(df, deriv_limit):
    #     """
    #     Remove data points that have a backward derivative greater than the deriv_limit that is passed to the function
    #     :param df:
    #     :param deriv_limit:
    #     :return:
    #     """
    #     # function[DATA, nout, ind] = backwardsTemporalDerivative(rad_vel, deriv_limit)
    #     # deriv_limit is expressed in velocity change(cm / s) over the course of an hour
    #
    #     # Convert the derivative limit from cm / s * hour to cm / s * s
    #     deriv_limit = deriv_limit / 3600  # units cm / s ^ 2
    #
    #     # use the diff function to calculate the derivative
    #     deriv = diff(rad_vel(:, 2))./ (24 * 60 * 60 * diff(rad_vel(:, 1)))
    #
    #     # the backward derivative is the output of the diff calculation with the last value removed
    #     bckwd_deriv = abs(deriv);
    #     bckwd_deriv(length(bckwd_deriv)) = []
    #
    #     # Find the indices of the backward derivatives that are greater than the derivative limit
    #     ind = find(bckwd_deriv > deriv_limit)
    #
    #     # find the number of outliers
    #     nout = length(ind);
    #
    #     # add 1 to the indices so that they match the indices of the vectors that you are screening
    #     ind = ind + 1;
    #
    #     # remove the outliers
    #     rad_vel(ind,:)=[];
    #
    #     # assign the new matrix to the variable DATA
    #     DATA = rad_vel;
    #     return df