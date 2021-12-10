import datetime as dt
import pandas as pd
import re
import xarray as xr
from hfradarpy.ctf import CTFParser
from hfradarpy.io.nc import make_encoding
from hfradarpy.common import create_dir
import numpy as np
import os
from pathlib import Path

import logging
logger = logging.getLogger(__name__)


def concat(wave_list, enhance=False):
    """
    This function takes a list of Wave objects or wave file paths and
    combines them along the time dimension using xarrays built-in concatenation
    routines.
    :param wave_list: list of wave files or Wave objects that you want to concatenate
    :return: waves concatenated into an xarray dataset by (time) or (time, dist)
    """
    wave_dict = {}
    for wave in wave_list:
        if not isinstance(wave, Waves):
            wave = Waves(wave)
        wave_dict[wave.file_name] = wave.to_xarray(enhance=enhance)

    ds = xr.concat(wave_dict.values(), 'time')
    return ds.sortby('time')


class Waves(CTFParser):
    """
    Waves Subclass.

    This class should be used when loading a CODAR wave (.wls) file. This class utilizes the generic LLUV class from
    ~/hfradarpy/ctf.py in order to load CODAR wave files
    """

    def __init__(self, fname, replace_invalid=True):
        logging.info('Loading wave file: {}'.format(fname))
        super().__init__(fname)

        if self._iscorrupt:
            return

        if self._tables['1']['data']['DIST'].isnull().all():
            df = self._tables['1']['data']
            self.data = df
            self.df_index = 'time'  # define index so pd.to_xarray function will automatically assign dimension and coordinates
        else:
            data_tables = []
            for key in self._tables.keys():
                df = self._tables[key]['data']
                data_tables.append(df)
            self.data = pd.concat(data_tables, axis=0)
            self.df_index = ['time', 'DIST']  # define two indices for multidimensional indexing.

        # Use separate date and time columns to create atetime column and drop those columns.
        self.data['time'] = self.data[['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)

        if not self.data.empty:
            if replace_invalid:
                self.replace_invalid_values()

    def __repr__(self):
        return "<Wave: {}>".format(self.file_name)

    def file_type(self):
        """Return a string representing the type of file this is."""
        return 'wave'

    def clean_header(self):
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


    def flag_wave_heights(self, min=0.2, max=5, remove=False):
        """
        Flag bad wave heights in Wave instance. This method labels wave heights between wave_min and wave_max as good,
        while labeling anything else bad
        :param wave_min: Minimum Wave Height - Waves above this will be considered good
        :param wave_max: Maximum Wave Height - Waves less than this will be considered good
        :param remove: Remove bad wave heights. Defaults to False
        :return:
        """
        boolean = self.data['MWHT'].between(min, max, inclusive='both')
        
        if not remove:
            self.data['mwht_flag'] = 1
            self.data['mwht_flag'] = self.data['mwht_flag'].where(boolean, other=4)
        elif remove:
            self.data = self.data[boolean]


    def is_valid(self):
        if self.data.empty:
            return False
        else:
            return True
    

    def to_xarray(self, enhance=False):
        """
        :param range_min:
        :param range_max:
        :return:
        """
        logging.info('Converting wave data to xarray dataset')

        # Set dataframe to indexes defined during class initialization
        tdf = self.data.set_index(self.df_index).drop(['TIME', 'TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'], axis=1)

        # Intitialize xarray dataset
        ds = tdf.to_xarray()

        # Assign header data to global attributes
        ds = ds.assign_attrs(self.metadata)
        
        if enhance is True:
            # global_attr = required_global_attributes(required_attributes, time_start, time_end)
            ds = self.enhance_xarray(ds)
            ds = xr.decode_cf(ds)

        return ds
    
    def enhance_xarray(self, xds):
        rename = dict()
        rename['MWHT'] = 'wave_height'
        rename['MWPD'] = 'wave_period'
        rename['WAVB'] = 'wave_bearing'
        rename['WNDB'] = 'wind_bearing'
        rename['ACNT'] = 'cross_spectra_averaged_count'
        rename['DIST'] = 'distance_from_origin'
        rename['RCLL'] = 'range_cell_result'
        rename['WDPT'] = 'doppler_points_used'
        rename['MTHD'] = 'wave_method'
        rename['FLAG'] = 'vector_flag'

        if 'PMWH' in self.data.keys():
            rename['PMWH'] = 'maximum_observable_wave_height'

        if 'WHNM' in self.data.keys():
            rename['WHNM'] = 'num_valid_source_wave_vectors'

        if 'WHSD' in self.data.keys():
            rename['WHSD'] = 'standard_deviation_of_wave_heights'

        # rename variables to something meaningful if they existin
        # in the xarray dataset
        existing_renames = { k: v for k, v in rename.items() if k in xds }
        xds = xds.rename(existing_renames)

        length = len(xds.time)
        lonlat = [float(x) for x in self.metadata['Origin'].split()]
        xds['lon'] = xr.DataArray(np.full(length, lonlat[1]), dims=('time'))
        xds['lat'] = xr.DataArray(np.full(length, lonlat[0]), dims=('time'))

        # set time attribute
        xds['time'].attrs['standard_name'] = 'time'
        xds['time'].attrs['long_name'] = 'Universal Time Coordinated (UTC) Time'

        # Set wave_height attributes
        xds['wave_height'].attrs['long_name'] = 'wave model height in meters'
        xds['wave_height'].attrs['standard_name'] = 'sea_surface_wave_significant_height'
        xds['wave_height'].attrs['units'] = 'm'
        xds['wave_height'].attrs['comment'] = 'wave model height in meters for every one of three waves'
        xds['wave_height'].attrs['valid_min'] = np.double(0)
        xds['wave_height'].attrs['valid_max'] = np.double(100)
        xds['wave_height'].attrs['coordinates'] = 'time'
        xds['wave_height'].attrs['grid_mapping'] = 'crs'
        xds['wave_height'].attrs['coverage_content_type'] = 'physicalMeasurement'

        # Set wave_period attributes
        xds['wave_period'].attrs['long_name'] = 'wave spectra period in seconds'
        xds['wave_period'].attrs['standard_name'] = 'sea_surface_wave_mean_period'
        xds['wave_period'].attrs['units'] = 's'
        xds['wave_period'].attrs['comment'] = 'wave spectra period in seconds'
        xds['wave_period'].attrs['valid_min'] = np.double(0)
        xds['wave_period'].attrs['valid_max'] = np.double(100)
        xds['wave_period'].attrs['coordinates'] = 'time'
        xds['wave_period'].attrs['grid_mapping'] = 'crs'
        xds['wave_period'].attrs['coverage_content_type'] = 'physicalMeasurement'

        # Set wave_bearing attributes
        xds['wave_bearing'].attrs['long_name'] = 'wave from direction in degrees'
        xds['wave_bearing'].attrs['standard_name'] = 'sea_surface_wave_from_direction'
        xds['wave_bearing'].attrs['units'] = 'degrees'
        xds['wave_bearing'].attrs['comment'] = 'wave from direction in degrees'
        xds['wave_bearing'].attrs['valid_min'] = np.double(0)
        xds['wave_bearing'].attrs['valid_max'] = np.double(360)
        xds['wave_bearing'].attrs['coordinates'] = 'time'
        xds['wave_bearing'].attrs['grid_mapping'] = 'crs'
        xds['wave_bearing'].attrs['coverage_content_type'] = 'physicalMeasurement'

        # Set wind_bearing attributes
        xds['wind_bearing'].attrs['long_name'] = 'wind from direction in degrees'
        xds['wind_bearing'].attrs['standard_name'] = 'sea_surface_wind_wave_from_direction'
        xds['wind_bearing'].attrs['units'] = 'degrees'
        xds['wind_bearing'].attrs['comment'] = 'wind from direction in degrees'
        xds['wind_bearing'].attrs['valid_min'] = np.double(0)
        xds['wind_bearing'].attrs['valid_max'] = np.double(360)
        xds['wind_bearing'].attrs['coordinates'] = 'time'
        xds['wind_bearing'].attrs['grid_mapping'] = 'crs'
        xds['wind_bearing'].attrs['coverage_content_type'] = 'physicalMeasurement'

        # Set lon attributes
        xds['lon'].attrs['long_name'] = 'Longitude'
        xds['lon'].attrs['standard_name'] = 'longitude'
        xds['lon'].attrs['short_name'] = 'lon'
        xds['lon'].attrs['units'] = 'degrees_east'
        xds['lon'].attrs['axis'] = 'X'
        xds['lon'].attrs['valid_min'] = np.double(-180.0)
        xds['lon'].attrs['valid_max'] = np.double(180.0)
        xds['lon'].attrs['grid_mapping'] = 'crs'

        # Set lat attributes
        xds['lat'].attrs['long_name'] = 'Latitude'
        xds['lat'].attrs['standard_name'] = 'latitude'
        xds['lat'].attrs['short_name'] = 'lat'
        xds['lat'].attrs['units'] = 'degrees_north'
        xds['lat'].attrs['axis'] = 'Y'
        xds['lat'].attrs['valid_min'] = np.double(-90.0)
        xds['lat'].attrs['valid_max'] = np.double(90.0)
        xds['lat'].attrs['grid_mapping'] = 'crs'

        # # add container variables that contain no data
        # xds = xds.assign(**dict(crs=False, instrument=False))
        
        # # Set crs attributes
        # xds['crs'].attrs['grid_mapping_name'] = 'latitude_longitude'
        # xds['crs'].attrs['inverse_flattening'] = 298.257223563
        # xds['crs'].attrs['long_name'] = 'Coordinate Reference System'
        # xds['crs'].attrs['semi_major_axis'] = '6378137.0'
        # xds['crs'].attrs['epsg_code'] = 'EPSG:4326'
        # xds['crs'].attrs['comment'] = 'http://www.opengis.net/def/crs/EPSG/0/4326'
        
        # xds['instrument'].attrs['long_name'] = 'Direction-finding high frequency radar antenna'
        # xds['instrument'].attrs['sensor_type'] = 'Direction-finding high frequency radar antenna'
        # xds['instrument'].attrs['make_model'] = self.metadata['Manufacturer']
        # xds['instrument'].attrs['serial_number'] = 1

        return xds
    
    def create_netcdf(self, filename, prepend_ext=False, enhance=True):
        """
        Create a compressed netCDF4 (.nc) file from the radial instance
        :param filename: User defined filename of radial file you want to save
        :return:
        """

        # If the filename does not have a .nc extension, we will add one.
        if not '.nc' in str(filename):
            filename = filename.with_suffix('.nc')

        # If the outputted file exists already, delete the existing file  
        if os.path.isfile(filename):
            os.remove(filename)

        create_dir(os.path.dirname(filename))

        # Convert pandas dataframe to xarray dataset using built-in pandas to_xarray function
        xds = self.to_xarray(enhance=enhance)

        # Check if dataset has distance_from_origin in coordinates. We will prepend the .nc extension 
        # with the appropriate name depending on whether the wave file is averaged or arranged by distance with manufacturer software
        if prepend_ext:
            if 'distance_from_origin' in xds.coords:
                pre_ext = 'ranged'
            else:
                pre_ext = 'averaged'
            # Change the extension to reflect the type of wave file
            filename = filename.with_suffix(f'.{pre_ext}.nc')

        # Pass through make_encoding function fo automatically 
        encoding = make_encoding(xds, comp_level=4, fillvalue=np.nan)
        encoding['time'] = dict(zlib=False, _FillValue=None)

        # Convert files to netcdf
        xds.to_netcdf(
            filename,
            encoding=encoding,
            format='netCDF4',
            engine='netcdf4',
            unlimited_dims=['time']
        )


    # def create_wls(self, filename):
    #     """
    #     Create a CODAR Wave (.wls) file from wave instance
    #     :param filename: User defined filename of wave file you want to save
    #     :return:
    #     """
    #     # # Ensure that the filename passed into the export function is not the same as the filename that we read in.
    #     # # We do not want to overwrite the original wave file by accident.
    #     # if self.full_file == str(filename):
    #     #     suffix = f'.mod{filename.suffix}'
    #     #     filename = filename.with_suffix(suffix)

    #     # if os.path.isfile(filename):
    #     #     os.remove(filename)

    #     create_dir(os.path.dirname(filename))
    #     rcopy = copy.deepcopy(self)
    #     with open(filename, 'w') as f:
    #         # Write header
    #         for metadata_key, metadata_value in self.metadata.items():
    #             if 'ProcessedTimeStamp' in metadata_key:
    #                 break
    #             else:
    #                 f.write('%{}: {}\n'.format(metadata_key, metadata_value))

    #         # Write data tables. Anything beyond the first table is commented out.
    #         for table in self._tables.keys():
    #             for table_key, table_value in self._tables[table].items():
    #                 if table_key != 'data':
    #                     if (table_key == 'TableType') & (table == '1'):
    #                         if 'QCTest' in self.metadata:
    #                             f.write('%QCFileVersion: 1.0.0\n')
    #                             f.write('%QCReference: Quality control reference: IOOS QARTOD HF Radar ver 1.0 May 2016\n')
    #                             f.write('%QCFlagDefinitions: 1=pass 2=not_evaluated 3=suspect 4=fail 9=missing_data\n')
    #                             f.write('%QCTestFormat: "test_name [qc_thresholds]: test_result"\n')

    #                             for test in self.metadata['QCTest']:
    #                                 f.write('%QCTest: {}\n'.format(test))
    #                         f.write('%{}: {}\n'.format(table_key, table_value))
    #                     elif table_key == 'TableColumns':
    #                         f.write('%TableColumns: {}\n'.format(len(self._tables[table]['data'].columns)))
    #                     elif table_key == 'TableColumnTypes':
    #                         f.write('%TableColumnTypes: {}\n'.format(' '.join(self._tables[table]['data'].columns.to_list())))
    #                     elif table_key == 'TableStart':
    #                         f.write('%{}: {}\n'.format(table_key, table_value))
    #                     elif table_key == '_TableHeader':
    #                         pass
    #                     else:
    #                         f.write('%{}: {}\n'.format(table_key, table_value))

    #             if 'datetime' in self._tables[table]['data'].keys():
    #                 self._tables[table]['data'] = self._tables[table]['data'].drop(['datetime'], axis=1)

    #             if table == '1':
    #                 # Fill NaN with 999.000 which is the standard fill value for codar lluv filesself._tables[table]['TableColumnTypes']
    #                 self.data = self.data.fillna(999.000)

    #                 try:
    #                     self.data['LOND'] = self.data['LOND'].apply(lambda x: "{:.7f}".format(x))
    #                     self.data['LATD'] = self.data['LATD'].apply(lambda x: "{:.7f}".format(x))
    #                     self.data['ESPC'] = self.data['ESPC'].apply(lambda x: "{:.3f}".format(x))
    #                     if 'ETMP' in self.data.columns:
    #                         self.data['ETMP'] = self.data['ETMP'].apply(lambda x: "{:.3f}".format(x))
    #                     self.data['BEAR'] = self.data['BEAR'].apply(lambda x: "{:.1f}".format(x))
    #                     self.data['HEAD'] = self.data['HEAD'].apply(lambda x: "{:.1f}".format(x))
    #                 except:
    #                     self = rcopy
    #                     print("Unexpected error in formatting one of these columns: LOND LATD ESPC ETMP BEAR HEAD")

    #                 # Convert _TableHeader to a new dataframe and concatenate to dataframe containing radial data
    #                 # This allows for the output format to follow CODARS CTF specifications
    #                 row_df = pd.DataFrame([self._tables['1']['_TableHeader'][1]], columns=self._tables['1']['_TableHeader'][0])
    #                 self.data.columns = self._tables['1']['_TableHeader'][0]
    #                 self.data = pd.concat([row_df, self.data], ignore_index=True)
    #                 self.data.insert(0, '%%', np.nan)  # Insert column at the beginning of dataframe of NaNs
    #                 self.data.iloc[0, self.data.columns.get_loc('%%')] = '%%'  # make the first row in the first column a '%%'

    #                 # Output data table to string
    #                 #self.data.to_string(f, index=False, justify='center', header=True, na_rep=' ')
    #                 self.data.temp = re.sub(' %%', '%%', self.data.to_string(index=False, justify='right', header=True, na_rep=' '))
    #                 f.write(self.data.temp)
    #             else:
    #                 self._tables[table]['data'].insert(0, '%%', '%')
    #                 self._tables[table]['data'] = self._tables[table]['data'].fillna(999.000)
    #                 self._tables[table]['data'].to_string(f, index=False, justify='center', header=True)

    #             if int(table) > 1:
    #                 f.write('\n%TableEnd: {}\n'.format(table))
    #             else:
    #                 f.write('\n%TableEnd: \n')
    #             f.write('%%\n')

    #         # Write footer containing processing information
    #         f.write('%ProcessedTimeStamp: {}\n'.format(self.metadata['ProcessedTimeStamp']))
    #         for tool in self.metadata['ProcessingTool']:
    #             f.write('%ProcessingTool: {}\n'.format(tool))
    #             # f.write('%{}: {}\n'.format(footer_key, footer_value))
    #         f.write('%End:')

    def export(self, filename, file_type='netcdf', prepend_ext=False):
        """
        Export wave file as either a codar .wls file or a netcdf .nc file
        :param filename: User defined filename of wave file you want to save
        :param file_type: Type of file to export wave: wave (default) or netcdf
        :param prepend_ext: <False>. Detect the distance method of wave file, and add 'ranged' or 'averaged' distance before .nc file extension (for NetCDF only)
        :return:
        """
        # Make sure filename is converted into a Path object 
        filename = Path(filename)

        if not self.is_valid():
            raise ValueError("Could not export ASCII data, the input file was invalid.")
        
        if file_type == 'wave':
            logging.info('Cannot create .wls file. Writing of.wls CTF files not implemented yet.')
            # self.create_wave(filename)
        elif file_type == 'netcdf':         
            self.create_netcdf(filename, prepend_ext)
