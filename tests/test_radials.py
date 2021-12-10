import unittest
from pathlib import Path

import numpy as np
import xarray as xr

from hfradarpy.radials import Radial
from hfradarpy.radials import concat as concatenate_radials

data_path = (Path(__file__).parent.with_name('examples') / 'data').resolve()
output_path = (Path(__file__).parent.with_name('examples') / 'output').resolve()


def test_codar_radial_to_tabular_netcdf():
    radial_file = data_path / 'radials' / 'ruv' / 'SEAB' / 'RDLi_SEAB_2019_01_01_0000.ruv'
    nc_file = output_path / 'radials' / 'nc' / 'tabular' / 'SEAB' / 'RDLi_SEAB_2019_01_01_0000.nc'

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.export(str(nc_file), file_type='netcdf-tabular')

    # Convert it to an xarray Dataset with no variable
    # or attribte enhancements
    xds2 = rad1.to_xarray_tabular(enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray_tabular(enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_codar_radial_to_multidimensional_netcdf():
    radial_file = data_path / 'radials' / 'ruv' / 'SEAB' / 'RDLi_SEAB_2019_01_01_0000.ruv'
    nc_file = output_path / 'radials' / 'nc' / 'multidimensional' / 'SEAB' / 'RDLi_SEAB_2019_01_01_0000.nc'

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.export(str(nc_file), file_type='netcdf-multidimensional')

    # Convert it to an xarray Dataset with no variable
    # or attribte enhancements
    xds2 = rad1.to_xarray_multidimensional(enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray_multidimensional(enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_wera_radial_to_tabular_netcdf():
    radial_file = data_path / 'radials' / 'ruv' / 'WERA' / 'RDL_csw_2019_10_24_162300.ruv'
    nc_file = output_path / 'radials' / 'nc' / 'tabular' / 'WERA' / 'RDL_csw_2019_10_24_162300.nc'

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.export(str(nc_file), file_type='netcdf-tabular')

    # Convert it to an xarray Dataset with no variable
    # or attribte enhancements
    xds2 = rad1.to_xarray_tabular(enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray_tabular(enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_wera_radial_to_multidimensional_netcdf():
    radial_file = data_path / 'radials' / 'ruv' / 'WERA' / 'RDL_csw_2019_10_24_162300.ruv'
    nc_file = output_path / 'radials' / 'nc' / 'multidimensional' / 'WERA' / 'RDL_csw_2019_10_24_162300.nc'

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.export(str(nc_file), file_type='netcdf-multidimensional')

    # Convert it to an xarray Dataset with no variable
    # or attribte enhancements
    xds2 = rad1.to_xarray_multidimensional(enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray_multidimensional(enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_wera_mask():
    radial_file = data_path / 'radials' / 'ruv' / 'WERA' / 'RDL_csw_2019_10_24_162300.ruv'
    rad1 = Radial(radial_file, mask_over_land=False, replace_invalid=False)
    # Total points before masking
    assert len(rad1.data) == 6327
    rad1.mask_over_land()
    # Make sure we subset the land points
    assert len(rad1.data) == 5745


def test_wera_qc():
    radial_file = data_path / 'radials' / 'ruv' / 'WERA' / 'RDL_csw_2019_10_24_162300.ruv'
    rad1 = Radial(radial_file, mask_over_land=False, replace_invalid=False)
    rad1.initialize_qc()
    assert len(rad1.data) == 6327
    rad1.mask_over_land()
    rad1.qc_qartod_radial_count()
    rad1.qc_qartod_valid_location()
    rad1.qc_qartod_maximum_velocity()
    rad1.qc_qartod_spatial_median()
    rad1.qc_qartod_avg_radial_bearing(reference_bearing=180)
    rad1.qc_qartod_primary_flag()
    assert len(rad1.data) == 5745
    assert 'QC07' in rad1.data
    assert 'QC08' not in rad1.data  # no VFLG column so we can't run it
    assert 'QC09' in rad1.data
    assert 'QC10' in rad1.data
    # assert 'QC11' in rad1.data  # temporal gradient test
    assert 'QC12' in rad1.data
    assert 'PRIM' in rad1.data


def test_wera_raw_to_quality_multidimensional_nc():
    radial_file = data_path / 'radials' / 'ruv' / 'WERA' / 'RDL_csw_2019_10_24_162300.ruv'
    nc_file = output_path / 'radials' / 'qc' / 'nc' / 'multidimensional' / 'WERA' / 'RDL_csw_2019_10_24_162300.nc'
    rad1 = Radial(radial_file, mask_over_land=False, replace_invalid=False)
    rad1.initialize_qc()
    rad1.mask_over_land()
    rad1.qc_qartod_radial_count()
    rad1.qc_qartod_valid_location()
    rad1.qc_qartod_maximum_velocity()
    rad1.qc_qartod_spatial_median()
    rad1.export(str(nc_file), file_type='netcdf-multidimensional')

    xds2 = rad1.to_xarray_multidimensional(enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        assert len(xds1.QCTest) == 3  # no VFLG column so one test not run
        # The two enhanced files should be identical
        assert xds1.identical(xds2)


def test_wera_raw_to_quality_tabular_nc():
    radial_file = data_path / 'radials' / 'ruv' / 'WERA' / 'RDL_csw_2019_10_24_162300.ruv'
    nc_file = output_path / 'radials' / 'qc' / 'nc' / 'tabular' / 'WERA' / 'RDL_csw_2019_10_24_162300.nc'
    rad1 = Radial(radial_file, mask_over_land=False, replace_invalid=False)
    rad1.initialize_qc()
    rad1.mask_over_land()
    rad1.qc_qartod_radial_count()
    rad1.qc_qartod_valid_location()
    rad1.qc_qartod_maximum_velocity()
    rad1.qc_qartod_spatial_median()
    rad1.export(str(nc_file), file_type='netcdf-tabular')

    xds2 = rad1.to_xarray_tabular(enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        assert len(xds1.QCTest) == 3  # no VFLG column so one test not run
        # The two enhanced files should be identical
        assert xds1.identical(xds2)


def test_miami_radial_multidimensional_nc():
    radial_file = data_path / 'radials' / 'ruv' / 'WERA' / 'RDL_UMiami_STF_2019_06_01_0000.hfrweralluv1.0'
    nc_file = output_path / 'radials' / 'nc' / 'multidimensional' / 'WERA' / 'RDL_UMiami_STF_2019_06_01_0000.nc'

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.export(str(nc_file), file_type='netcdf-multidimensional')

    # Convert it to an xarray Dataset with no variable
    # or attribte enhancements
    xds2 = rad1.to_xarray_multidimensional(enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray_multidimensional(enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_miami_radial_tabular_nc():
    radial_file = data_path / 'radials' / 'ruv' / 'WERA' / 'RDL_UMiami_STF_2019_06_01_0000.hfrweralluv1.0'
    nc_file = output_path / 'radials' / 'nc' / 'tabular' / 'WERA' / 'RDL_UMiami_STF_2019_06_01_0000.nc'

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.export(str(nc_file), file_type='netcdf-tabular')

    # Convert it to an xarray Dataset with no variable
    # or attribte enhancements
    xds2 = rad1.to_xarray_tabular(enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray_tabular(enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


class TestCombineRadials(unittest.TestCase):

    def setUp(self):
        self.file_paths = list(
            (data_path / 'radials' / 'ruv' / 'SEAB').glob('*.ruv')
        )

        self.radial_files = [
            str(r) for r in self.file_paths
        ]

        self.radial_objects = [
            Radial(str(r)) for r in self.radial_files
        ]

        # Select even indexed file_paths and odd indexed radial objects
        # into one array of mixed content types for concating
        self.radial_mixed = self.radial_files[::2] + self.radial_objects[1:][::2]

    def test_concat_radial_objects(self):
        combined = concatenate_radials(self.radial_objects)
        assert combined.time.size == len(self.file_paths)
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )

    def test_concat_radial_files(self):
        combined = concatenate_radials(self.radial_files)
        assert combined.time.size == len(self.file_paths)
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )

    def test_concat_mixed_radials(self):
        combined = concatenate_radials(self.radial_mixed)
        assert combined.time.size == len(self.file_paths)
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )

    def test_concat_mixed_radials_enhance(self):
        # Select even indexed file_paths and odd indexed radial objects
        # into one array of mixed content types for concating
        combined = concatenate_radials(self.radial_mixed, enhance=True)
        assert combined.time.size == len(self.file_paths)
        # Make sure the dataset was sorted by time
        assert np.array_equal(
            combined.time.values,
            np.sort(combined.time.values)
        )
 