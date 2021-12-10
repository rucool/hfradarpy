import unittest
from pathlib import Path
import numpy as np
import xarray as xr
from hfradarpy.waves import Waves, concat

data_root = (Path(__file__).parent.with_name('examples') / 'data').resolve()
output_path = (Path(__file__).parent.with_name('examples') / 'output').resolve()


def test_codar_averaged_wave_to_netcdf():
    wave_file = data_root / 'waves' / 'wls' / 'SEAB' / 'WVLM_SEAB_2019_01_01_0000.wls'
    nc_file = output_path / 'waves' / 'nc' / 'SEAB' / 'WVLM_SEAB_2019_01_01_0000.nc'

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    wav1 = Waves(wave_file)
    wav1.export(str(nc_file), file_type='netcdf')

    # Convert it to an xarray Dataset with no variable
    # or attribte enhancements
    xds2 = wav1.to_xarray(enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = wav1.to_xarray(enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_codar_ranged_wave_to_netcdf():
    wave_file = data_root / 'waves' / 'wls' / 'SEAB' / 'WVLR_SEAB_2019_01_01_0000.wls'
    nc_file = output_path / 'waves' / 'nc' / 'SEAB' / 'WVLR_SEAB_2019_01_01_0000.nc'

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    wav1 = Waves(wave_file)
    wav1.export(str(nc_file), file_type='netcdf')

    # Convert it to an xarray Dataset with no variable
    # or attribte enhancements
    xds2 = wav1.to_xarray(enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = wav1.to_xarray(enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_wave_filter():
    wave_file = data_root / 'waves' / 'wls' / 'SEAB' / 'WVLR_SEAB_2019_01_01_0000.wls'
    
    # Load the wave file
    wav1 = Waves(wave_file, replace_invalid=False)
    assert len(wav1.data) == 12645

    # Check that flag_wave_heights is removing the bad data (and not creating an extra column of flags) 
    wav1.flag_wave_heights(2, 5, remove=True)
    assert len(wav1.data) < 12645

    # Check that flag_wave_heights is adding an extra column of flags
    wav1.flag_wave_heights(2, 5)
    assert 'mwht_flag' in wav1.data



class TestCombineWaves(unittest.TestCase):

    def setUp(self):
        self.wvlm_paths = list(
            (data_root / 'waves' / 'wls' / 'SEAB').glob('WVLM*.wls')
        )

        self.wvlm_files = [
            str(r) for r in self.wvlm_paths
        ]

        self.wvlm_objects = [
            Waves(str(r)) for r in self.wvlm_files
        ]

        self.wvlr_paths = list(
            (data_root / 'waves' / 'wls' / 'SEAB').glob('WVLR*.wls')
        )

        self.wvlr_files = [
            str(r) for r in self.wvlr_paths
        ]

        self.wvlr_objects = [
            Waves(str(r)) for r in self.wvlr_files
        ]

        # Select even indexed file_paths and odd indexed radial objects
        # into one array of mixed content types for concating
        self.wvlm_mixed = self.wvlm_files[::2] + self.wvlm_objects[1:][::2]
        self.wvlr_mixed = self.wvlr_files[::2] + self.wvlr_objects[1:][::2]

    def test_concat_wave_wvlm_objects(self):
        combined = concat(self.wvlm_objects)
        assert combined.time.size == 8524
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )

    def test_concat_wave_wvlm_files(self):
        combined = concat(self.wvlm_files)
        assert combined.time.size == 8524
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )

    def test_concat_mixed_wvlm_waves(self):
        combined = concat(self.wvlm_mixed)
        assert combined.time.size == 8524
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )

    def test_concat_mixed_wvlm_waves_enhance(self):
        # Select even indexed file_paths and odd indexed radial objects
        # into one array of mixed content types for concating
        combined = concat(self.wvlm_mixed, enhance=True)
        assert combined.time.size == 8524
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )

    def test_concat_wave_wvlr_objects(self):
        combined = concat(self.wvlr_objects)
        assert combined.time.size == 8522
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )

    def test_concat_wave_wvlr_files(self):
        combined = concat(self.wvlr_files)
        assert combined.time.size == 8522
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )

    def test_concat_mixed_wvlr_waves(self):
        combined = concat(self.wvlr_mixed)
        assert combined.time.size == 8522
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )

    def test_concat_mixed_wvlr_waves_enhance(self):
        # Select even indexed file_paths and odd indexed radial objects
        # into one array of mixed content types for concating
        combined = concat(self.wvlr_mixed, enhance=True)
        assert combined.time.size == 8522
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )
 