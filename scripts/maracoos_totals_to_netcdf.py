import pandas as pd
from codar_processing.common import list_files, list_to_dataframe
from functions import mat_to_netcdf4

data_dir = '/home/codaradm/data/totals/maracoos/oi/mat/5MHz'
save_dir = '/home/michaesm/codar/data/nc'
hfr_grid = '../totals/grid_files/OI_6km_Grid_Extend.txt'
start_time = pd.Timestamp(2009, 11, 7, 0, 0, 0)
end_time = pd.Timestamp(2018, 4, 26, 0, 0, 0)

avoid = ('ideal', 'measured')  # avoid these subfolders
types = ['*.mat']  # find only files with these types

user_attributes = dict(title='MARACOOS 6km Sea Surface Currents',
                       naming_authority='edu.rutgers.marine.rucool',
                       comment='Network maintained by MARACOOS. For oi_* global attribute explanations, see references attribute',
                       acknowledgment='This data is provided by the Mid-Atlantic Regional Association Coastal Ocean Observing System (MARACOOS). Funding is provided by the U.S. Integration Ocean Observing System (IOOS).',
                       standard_name_vocabulary='CF Standard Name Table v41',
                       creator_name='Michael Smith',
                       creator_email='michaesm@marine.rutgers.edu',
                       creator_url='rucool.marine.rutgers.edu',
                       institution='Center for Ocean Observing and Leadership, Department of Marine & Coastal Sciences, Rutgers University',
                       project='Mid-Atlantic Regional Association Coastal Ocean Observing System - High Frequency Radar Sea Surface Current Mapping',
                       sea_name='Mid-Atlantic Bight',
                       creator_type='person',
                       creator_institution='Rutgers University',
                       contributor_name='Scott Glenn, Josh Kohut, Hugh Roarty, Ethan Handel, Michael Smith, Laura Nazzaro, Teresa Updyke, Larry Atkinson, Rich Arena, Wendell Brown, Mike Muglia, Harvey Seim, Sara Haines, Tony Whipple',
                       contributor_role='Principal Investigator, Principal Investigator, Principal Investigator, Hardware Maintenance, Data Manager, Data Manager, Principal Investigator, Principal Investigator, Hardware Maintenance, Principal Investigator, Hardware Maintenance, Principal Investigator, Data Manager, Hardware Maintenance',
                       platform='MARACOOS HF Radar 5MHz Network',
                       instrument='Network includes CODAR sites AMAG, ASSA, BLCK, BRIG, CEDR, CORE, DUCK, FARO, HATY, HEMP, HOOK, LISL, LOVE, MABO, MRCH, MVCO, NANT, NAUS, PYFC, and WILD',
                       references='http://maracoos.org/node/146 https://rucool.marine.rutgers.edu/facilities https://rucool.marine.rutgers.edu/data',
                       summary='Optimally Interpolated Total Vectors calculated by HFRProgs toolbox using MATLAB. Mercator lat/lon projection',
                       ncei_template_version='NCEI_NetCDF_Grid_Template_v2.0',
                       history='Hourly codar radial data combined into one hourly file containing vectors.',
                       cdm_data_type='Grid',
                       source='CODAR SeaSonde Surface Current Mapping Device',
                       processing_level='Level 3',
                       keywords='Environmental Advisories > Marine Advisories > Marine Weather/Forecast, Oceans > Coastal Processes, Oceans > Ocean Circulation, Oceans > Ocean Waves, Oceans > Ocean Winds, Oceans > Ocean Tides, Spectral/Engineering > Radar',
                       publisher_name='NOAA National Centers for Environmental Information',
                       publisher_email='ncei.info@noaa.gov',
                       publisher_url='www.ncei.noaa.gov')


# load csv file containing the grid
grid = pd.read_csv(hfr_grid, sep=',', header=None, names=['lon', 'lat'], delim_whitespace=True)

file_list = list_files(types, data_dir, avoid)
df = list_to_dataframe(file_list)

df = df[start_time: end_time]

for row in df.itertuples():
    mat_to_netcdf4.main(grid, row.file, save_dir, user_attributes)