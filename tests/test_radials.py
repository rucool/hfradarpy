from pathlib import Path

import numpy as np
import xarray as xr

from codar_processing.src.radials import Radial, concatenate_radials

data_path = (Path(__file__).parent.with_name('codar_processing') / 'data').resolve()
output_path = (Path(__file__).parent.with_name('output')).resolve()


def test_radial_to_netcdf():
    radial_file = data_path / 'radials' / 'SEAB' / 'RDLi_SEAB_2019_01_01_0000.ruv'
    nc1_file = output_path / 'radials_nc' / 'SEAB' / '1_RDLi_SEAB_2019_01_01_0000.nc'
    #nc2_file = output_path / 'radials_nc' / 'SEAB' / '2_RDLi_SEAB_2019_01_01_0000.nc'

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.export(str(nc1_file), file_type='netcdf')

    # Convert it to an xarray Dataset with no variable
    # or attribte enhancements
    xds2 = rad1.to_xarray(enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray(enhance=True)

    with xr.open_dataset(nc1_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


class TestCombineRadials:

    file_paths = list(
        (data_path / 'radials' / 'SEAB').glob('*.ruv')
    )

    radial_files = [
        str(r) for r in file_paths
    ]

    radial_objects = [
        Radial(str(r)) for r in radial_files
    ]

    # Select even indexed file_paths and odd indexed radial objects
    # into one array of mixed content types for concating
    radial_mixed = radial_files[::2] + radial_objects[1:][::2]

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
