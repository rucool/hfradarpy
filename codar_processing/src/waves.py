import datetime as dt
import pandas as pd
import re
import xarray as xr
from codar_processing.src.common import CTFParser


def concatenate_waves(wave_list):
    """
    This function takes a list of radial files. Loads them all separately using the Wave object and then combines
    them along the time dimension using xarrays built-in concatenation routines.
    :param wave_list: list of radial files that you want to concatenate
    :return: wave files concatenated into an xarray dataset by range, bearing, and time
    """

    wave_dict = {}
    for each in sorted(wave_list):
        wave = Waves(each, multi_dimensional=True)
        wave_dict[wave.file_name] = wave.ds

    ds = xr.concat(wave_dict.values(), 'time')
    return ds


class Waves(CTFParser):
    """
    Waves Subclass.

    This class should be used when loading a CODAR wave (.wls) file. This class utilizes the generic LLUV class from
    ~/codar_processing/common.py in order to load CODAR wave files
    """

    def __init__(self, fname, replace_invalid=True, multi_dimensional=True):
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
                      WHSD='standard_deviation_of_wave_heights')

        CTFParser.__init__(self, fname)
        if self._tables['1']['data']['DIST'].isnull().all():
            df = self._tables['1']['data']
            self.data = df
            index = 'datetime'  # define index so pd.to_xarray function will automatically assign dimension and coordinates
        else:
            data_tables = []
            for key in self._tables.keys():
                df = self._tables[key]['data']
                data_tables.append(df)
            self.data = pd.concat(data_tables, axis=0)
            index = ['datetime', 'DIST']  # define two indices for multidimensional indexing.

        # Use separate date and time columns to create datetime column and drop those columns.
        self.data['datetime'] = self.data[['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)

        if replace_invalid:
            self.replace_invalid_values()

        if multi_dimensional:
            # Set index of dataframe and also drop columns that we don't need to see anymore
            self.data = self.data.set_index(index).drop(['TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'], axis=1)

            # Convert from pandas dataframe into an xarray dataset
            self.ds = self.data.to_xarray()

            # rename variables to something meaningful
            self.ds.rename(rename, inplace=True)

            # Clean up wave header and assign header data to global attributes
            self.clean_wave_header()
            self.data = self.ds.assign_attrs(self.metadata)

    def clean_wave_header(self):
        """
        Cleans the header data from the wave data for proper input into MySQL database
        """

        keep = ['TimeCoverage', 'WaveMinDopplerPoints', 'AntennaBearing', 'DopplerCells', 'TransmitCenterFreqMHz',
                'CTF', 'TableColumnTypes', 'TimeZone', 'WaveBraggPeakDropOff', 'RangeResolutionKMeters',
                'CoastlineSector', 'WaveMergeMethod', 'RangeCells', 'WaveBraggPeakNull', 'WaveUseInnerBragg',
                'BraggSmoothingPoints', 'Manufacturer', 'TimeStamp', 'FileType', 'TableRows', 'BraggHasSecondOrder',
                'Origin', 'MaximumWavePeriod', 'UUID', 'WaveBraggNoiseThreshold', 'TransmitBandwidthKHz', 'Site',
                'TransmitSweepRateHz', 'WaveBearingLimits', 'WavesFollowTheWind', 'CurrentVelocityLimit']

        key_list = list(self.metadata.keys())
        for key in key_list:
            if not key in keep:
                del self.metadata[key]

        for k, v in self.metadata.items():
            if 'Site' in k:
                self.metadata[k] = ''.join(e for e in v if e.isalnum())
            elif 'TimeStamp' in k:
                t_list = v.split()
                t_list = [int(s) for s in t_list]
                self.metadata[k] = dt.datetime(t_list[0], t_list[1], t_list[2], t_list[3], t_list[4], t_list[5]).strftime(
                    '%Y-%m-%d %H:%M:%S')
            elif k in ('TimeCoverage', 'RangeResolutionKMeters'):
                self.metadata[k] = re.findall("\d+\.\d+", v)[0]
            elif k in ('WaveMergeMethod', 'WaveUseInnerBragg', 'WavesFollowTheWind'):
                self.metadata[k] = re.search(r'\d+', v).group()
            elif 'TimeZone' in k:
                self.metadata[k] = re.search('"(.*)"', v).group(1)
            elif k in ('WaveBearingLimits', 'CoastlineSector'):
                bearings = re.findall(r"[-+]?\d*\.\d+|\d+", v)
                self.metadata[k] = ', '.join(e for e in bearings)
            else:
                continue

    def file_type(self):
        """Return a string representing the type of file this is."""
        return 'wave'

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

    def is_valid(self):
        if self.data.empty:
            return False
        else:
            return True