import datetime as dt
import numpy as np
import pandas as pd
import re
# import xarray as xr
# from codar_processing.calc import reckon
from codar_processing.common import LLUVParser


class Waves(LLUVParser):
    """
    Waves Subclass.

    This class should be used when loading a CODAR wave (.wls) file. This class utilizes the generic LLUV class from
    ~/codar_processing/common.py in order to load CODAR wave files
    """

    def __init__(self, fname, replace_invalid=True, n_dimensional=True):
        rename = dict(datetime='time',
                      MWHT='wave_height',
                      MWPD='wave_period',
                      WAVB='wave_bearing',
                      WNDB='wind_bearing',
                      PMWH='maximum_observable_wave_height',
                      ACNT='cross_spectra_averaged_count',
                      DIST='distance_from_origin',
                      RCLL='range_cell_result',
                      WDPT='doppler_points_used',
                      MTHD='wave_method',
                      FLAG='vector_flag',
                      WHNM='num_valid_source_wave_vectors',
                      WHSD='standard_deviation_of_valid_wave_heights')

        LLUVParser.__init__(self, fname)
        if self._tables['1']['data']['DIST'].isnull().all():
            self.data = self._tables['1']['data']
            # new = True
        else:
            data_tables = []
            for key in self._tables.keys():
                data_tables.append(self._tables[key]['data'])
            self.data = pd.concat(data_tables, axis=0)
            # new = False

        # Cleanup pandas dataframe
        self.data['datetime'] = self.data[['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)
        self.data = self.data.set_index('datetime')

        # Drop extraneous columns
        self.data = self.data.drop(['TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'], axis=1)

        if replace_invalid:
            self.replace_invalid_values()

        if n_dimensional:
            # Convert from pandas dataframe into an xarray dataset
            self.data = self.data.to_xarray()

            # rename variables to something meaningful
            self.data.rename(rename, inplace=True)

            # Clean up wave header and assign header data to global attributes
            self.clean_wave_header()
            self.data = self.data.assign_attrs(self._metadata)

    def clean_wave_header(self):
        """
        Cleans the header data from the wave data for proper input into MySQL database
        :param head_dict: dictionary containing the header data
        :return: dictionary containing the cleaned header data
        :rtype:
        """

        keep = ['TimeCoverage', 'WaveMinDopplerPoints', 'AntennaBearing', 'DopplerCells', 'TransmitCenterFreqMHz',
                'CTF', 'TableColumnTypes', 'TimeZone', 'WaveBraggPeakDropOff', 'RangeResolutionKMeters',
                'CoastlineSector', 'WaveMergeMethod', 'RangeCells', 'WaveBraggPeakNull', 'WaveUseInnerBragg',
                'BraggSmoothingPoints', 'Manufacturer', 'TimeStamp', 'FileType', 'TableRows', 'BraggHasSecondOrder',
                'Origin', 'MaximumWavePeriod', 'UUID', 'WaveBraggNoiseThreshold', 'TransmitBandwidthKHz', 'Site',
                'TransmitSweepRateHz', 'WaveBearingLimits', 'WavesFollowTheWind', 'CurrentVelocityLimit']

        key_list = list(self._metadata.keys())
        for key in key_list:
            if not key in keep:
                del self._metadata[key]

        for k, v in self._metadata.items():
            if 'Site' in k:
                self._metadata[k] = ''.join(e for e in v if e.isalnum())
            elif 'TimeStamp' in k:
                t_list = v.split()
                t_list = [int(s) for s in t_list]
                self._metadata[k] = dt.datetime(t_list[0], t_list[1], t_list[2], t_list[3], t_list[4], t_list[5]).strftime(
                    '%Y-%m-%d %H:%M:%S')
            elif k in ('TimeCoverage', 'RangeResolutionKMeters'):
                self._metadata[k] = re.findall("\d+\.\d+", v)[0]
            elif k in ('WaveMergeMethod', 'WaveUseInnerBragg', 'WavesFollowTheWind'):
                self._metadata[k] = re.search(r'\d+', v).group()
            elif 'TimeZone' in k:
                self._metadata[k] = re.search('"(.*)"', v).group(1)
            elif k in ('WaveBearingLimits', 'CoastlineSector'):
                bearings = re.findall(r"[-+]?\d*\.\d+|\d+", v)
                self._metadata[k] = ', '.join(e for e in bearings)
            else:
                self._metadata

    def file_type(self):
        """Return a string representing the type of file this is."""
        return 'wave'

    # def flag_wave_heights(self, wave_min=0.2, wave_max=5):
    #     """
    #     Flag bad wave heights in Wave instance. This method labels wave heights between wave_min and wave_max good,
    #     while labeling anything else bad
    #     :param wave_min: Minimum Wave Height - Waves above this will be considered good
    #     :param wave_max: Maximum Wave Height - Waves less than this will be considered good
    #     :return:
    #     """
    #     # self.data['mwht_flag'] = 1
    #     # boolean = self.data['MWHT'].between(wave_min, wave_max, inclusive=True)
    #     # self.data['mwht_flag'] = self.data['mwht_flag'].where(boolean, other=4)
    #     # test = self.data['wave_height'].where((wave_min < self.data['wave_height']) & (self.data['wave_height'] < wave_max))

    def is_valid(self):
        if self.data.empty:
            return False
        else:
            return True

    def remove_bad_wave_heights(self, wave_min=0.2, wave_max=5):
        """
        Remove any data that lies outside of wave heights at wave_min and wave_max
        :param wave_min: Minimum Wave Height - Waves above this will be considered good
        :param wave_max: Maximum Wave Height - Waves less than this will be considered good
        :return:
        """
        self.flag_wave_heights(wave_min, wave_max)
        self.data = self.data[self.data['mwht_flag'] == 1]
        self.data = self.data.drop(['mwht_flag'], axis=1)

    def replace_invalid_values(self, values=[999.00, 1080.0]):
        # Convert 999.00 and 1080.0 CODAR 'uncalculable values' to NaN
        self.data.replace(values, np.nan, inplace=True)