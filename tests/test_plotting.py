from pathlib import Path
from codar_processing.plotting.plot_nc import plot_totals as nctotalsplot

data_path = (Path(__file__).parent.with_name('codar_processing') / 'data').resolve()
output_path = (Path(__file__).parent.with_name('output')).resolve()


def test_nc_plot():
    ncfile = data_path / 'radials_nc' / 'SEAB' / 'RDLi_SEAB_2019_01_01_0000.nc'
    nctotalsplot(ncfile, output_path, sub=10)


def test_another_one():
    ncfile = data_path / 'totals' / 'oi' / 'nc' / 'hourly' / 'RU_MARA_20190101T000000Z.nc'
    nctotalsplot(ncfile, output_path, sub=3)
