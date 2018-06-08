import datetime as dt
import pandas as pd
import re
from codar_processing.common import LLUVParser


class Waves(LLUVParser):
    """
    Waves Subclass.

    This class should be used when loading a CODAR wave (.wls) file. This class utilizes the generic LLUV class from
    ~/codar_processing/common.py in order to load CODAR wave files
    """

    def __init__(self, fname):
        LLUVParser.__init__(self, fname)
        if self.tables['1']['data']['DIST'].isnull().all():
            self.data = self.tables['1']['data']
        else:
            data_tables = []
            for key in self.tables.keys():
                data_tables.append(self.tables[key]['data'])
            self.data = pd.concat(data_tables, axis=0)
        self.data['datetime'] = self.data[['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)

    def file_type(self):
        """Return a string representing the type of file this is."""
        return 'wave'

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

        key_list = list(self.header.keys())
        for key in key_list:
            if not key in keep:
                del self.header[key]

        for k, v in self.header.items():
            if 'Site' in k:
                self.header[k] = ''.join(e for e in v if e.isalnum())
            elif 'TimeStamp' in k:
                t_list = v.split()
                t_list = [int(s) for s in t_list]
                self.header[k] = dt.datetime(t_list[0], t_list[1], t_list[2], t_list[3], t_list[4], t_list[5]).strftime(
                    '%Y-%m-%d %H:%M:%S')
            elif k in ('TimeCoverage', 'RangeResolutionKMeters'):
                self.header[k] = re.findall("\d+\.\d+", v)[0]
            elif k in ('WaveMergeMethod', 'WaveUseInnerBragg', 'WavesFollowTheWind'):
                self.header[k] = re.search(r'\d+', v).group()
            elif 'TimeZone' in k:
                self.header[k] = re.search('"(.*)"', v).group(1)
            elif k in ('WaveBearingLimits', 'CoastlineSector'):
                bearings = re.findall(r"[-+]?\d*\.\d+|\d+", v)
                self.header[k] = ', '.join(e for e in bearings)
            else:
                self.header

    def flag_wave_heights(self, wave_min=0.2, wave_max=5):
        """
        Flag bad wave heights in Wave instance. This method labels wave heights between wave_min and wave_max good,
        while labeling anything else bad
        :param wave_min: Minimum Wave Height - Waves above this will be considered good
        :param wave_max: Maximum Wave Height - Waves less than this will be considered good
        :return:
        """
        self.data['mwht_flag'] = 1
        boolean = self.data['MWHT'].between(wave_min, wave_max, inclusive=True)
        self.data['mwht_flag'] = self.data['mwht_flag'].where(boolean, other=4)

    def remove_bad_data(self):
        """
        Remove any data that contains NaN. In the original RUV file, these usually refer to data with 999 or 1080 values
        :return:
        """
        self.data = self.data.dropna()

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