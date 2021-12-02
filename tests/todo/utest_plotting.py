from pathlib import Path
from hfradarpy.plot.plot_nc import plot_totals as nctotalsplot
from hfradarpy.plot.plot_nc import plot_radials as ncradialsplot

data_path = (Path(__file__).parent.with_name('hfradarpy') / 'data').resolve()
output_path = (Path(__file__).parent.with_name('output')).resolve()
types = ['velocity', 'motion']


def test_codar_radials_nc_tabular_plot():
    ncfile = data_path / 'radials' / 'nc' / 'tabular' / 'SEAB' / 'RDLi_SEAB_2019_01_01_0000.nc'
    for t in types:
        output_file = output_path / 'plots' / 'tabular' / 'RDLi_SEAB_2019_01_01_0000-{}.png'.format(t)
        ncradialsplot(ncfile, output_file=output_file, sub=1, plot_type=t)


def test_codar_radials_nc_multidimensional_plot():
    ncfile = data_path / 'radials' / 'nc' / 'multidimensional' / 'SEAB' / 'RDLi_SEAB_2019_01_01_0000.nc'
    for t in types:
        output_file = output_path / 'plots' / 'multidimensional' / 'RDLi_SEAB_2019_01_01_0000-{}.png'.format(t)
        ncradialsplot(ncfile, output_file=output_file, sub=1, plot_type=t)


def test_wera_radials_plot():
    ncfile = data_path / 'radials' / 'nc' / 'multidimensional' / 'WERA' / 'RDL_csw_2019_10_24_162300.nc'
    for t in types:
        output_file = output_path / 'plots' / 'wera' / 'RDL_csw_2019_10_24_162300-{}.png'.format(t)
        ncradialsplot(ncfile, output_file=output_file, sub=1, plot_type=t)


def test_totals_plot():
    ncfile = data_path / 'totals' / 'oi' / 'nc' / 'hourly' / 'RU_MARA_20190101T000000Z.nc'
    output_file = output_path / 'plots' / 'totals' / 'RU_MARA_20190101T000000Z.png'
    nctotalsplot(ncfile, output_file=output_file, sub=2)
