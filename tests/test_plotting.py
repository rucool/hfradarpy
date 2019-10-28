from pathlib import Path
from codar_processing.plotting.plot_nc import plot_totals as nctotalsplot
from codar_processing.plotting.plot_nc import plot_radials as ncradialsplot

data_path = (Path(__file__).parent.with_name('codar_processing') / 'data').resolve()
output_path = (Path(__file__).parent.with_name('output')).resolve()


def test_codar_radials_plot():
    ncfile = data_path / 'radials_nc' / 'SEAB' / 'RDLi_SEAB_2019_01_01_0000.nc'
    ncradialsplot(ncfile, output_file=output_path / 'RDLi_SEAB_2019_01_01_0000.png', sub=1)


def test_wera_radials_plot():
    ncfile = data_path / 'radials_nc' / 'WERA' / 'RDL_csw_2019_10_24_162300.nc'
    ncradialsplot(ncfile, output_file=output_path / 'RDL_csw_2019_10_24_162300.png', sub=2)


def test_totals_plot():
    ncfile = data_path / 'totals' / 'oi' / 'nc' / 'hourly' / 'RU_MARA_20190101T000000Z.nc'
    nctotalsplot(ncfile, output_file=output_path / 'RU_MARA_20190101T000000Z.png' , sub=2)
