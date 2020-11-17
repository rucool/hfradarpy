import logging
import os
import pandas as pd
import sys
from hfradar.src.common import list_files, list_to_dataframe
from hfradar.methods.totals import hfrprogs_mat_to_netcdf4
from glob import glob
import datetime as dt

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

main_dir = '/Users/mikesmith/Documents/projects/swarm/data/totals/'
hfr_grid = '/Users/mikesmith/Documents/GitHub/rucool/codar_processing/codar_processing/data/grid_files/palmer_deep_grid_1km.txt'

end_time = dt.datetime.utcnow()
start_time = end_time - dt.timedelta(days=7)

process = ('ideal', 'measured', 'bestchoice')  # process totals in these subfolders
types = ['*.mat']  # find only files with these types


user_attributes = dict(title='Palmer Deep Antarctica 1.0 km Sea Surface Currents',
                       naming_authority='edu.rutgers.marine.rucool',
                       comment='Network maintained by University of Alaska - Fairbanks and Rutgers University. For oi_* global attribute explanations, see references attribute',
                       acknowledgment='This data is provided by the Rutger University, Center for Ocean Observing Leadership. Funding is provided by the U.S. National Science Foundation (NSF).',
                       standard_name_vocabulary='CF Standard Name Table v41',
                       creator_name='Michael Smith',
                       creator_email='michaesm@marine.rutgers.edu',
                       creator_url='rucool.marine.rutgers.edu',
                       institution='Center for Ocean Observing and Leadership, Department of Marine & Coastal Sciences, Rutgers University',
                       project='Project SWARM - High Frequency Radar Sea Surface Current Mapping',
                       sea_name='Southern Ocean',
                       creator_type='person',
                       creator_institution='Rutgers University',
                       contributor_name='Josh Kohut, Peter Winsor, Hank Statscewich, Hugh Roarty, Ethan Handel, Michael Smith, Laura Nazzaro',
                       contributor_role='Principal Investigator, Principal Investigator, Hardware Maintenance, Manager, Maintenance, Data Manager, Data Manager',
                       platform='Palmer Deep, Antarctica HF Radar 25MHz Network',
                       instrument='Network includes CODAR sites PALM, JOUB, and WAUW',
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
logging.debug('{} - Reading grid file'.format(hfr_grid))
grid = pd.read_csv(hfr_grid, sep=', ', header=None, names=['lon', 'lat'])
logging.debug('{} - Grid file loaded'.format(hfr_grid))


file_list = glob(os.path.join(main_dir, '**', dt.datetime.utcnow().strftime('%Y_%m'), '*.mat'), recursive=True)
df = list_to_dataframe(file_list)
df = df[start_time: end_time]

for row in df.itertuples():
    logging.info('{} - Converting MAT file to netCDF4 format'.format(os.path.basename(row.file)))
    if '/lsq/' in row.file:
        hfrprogs_mat_to_netcdf4.main(grid, row.file, os.path.join(os.path.dirname(row.file), 'nc'), user_attributes, method='lsq')
    elif '/oi/' in row.file:
        hfrprogs_mat_to_netcdf4.main(grid, row.file, os.path.join(os.path.dirname(row.file), 'nc'), user_attributes)
    else:
        logging.info('{} - cannot determine whether oi or lsq'.format(os.path.basename(row.file)))



