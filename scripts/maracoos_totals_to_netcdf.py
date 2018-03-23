import glob
import os
import pandas as pd
from tools import mat_to_netcdf4

main_dir = '/home/codaradm/data/totals/maracoos/oi/mat/5MHz'
save_dir = '/home/michaesm/codar/nc/maracoos'
hfr_grid = '../totals/grid_files/OI_6km_Grid_Extend.txt'

avoid = ('ideal', 'measured') # avoid these subfolders
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
                       contributor_name='Scott Glenn, Josh Kohut, Hugh Roarty, Ethan Handel, Michael Smith, Laura Nazzaro, Teresa Updyke, Larry Atkinson, Rich Arena, Wendell Brown, Mike Muglia, Harvey Seim',
                       contributor_role='Principal Investigator, Principal Investigator, Principal Investigator, Hardware Maintenance, Data Manager, Data Manager, Hardware Maintenance, Principal Investigator, Hardware Maintenance, Principal Investigator, Hardware Maintenance, Principal Investigator',
                       platform='MARACOOS HF Radar 5MHz Network',
                       instrument='Network includes CODAR sites AMAG, ASSA, BLCK, BRIG, CEDR, CORE, DUCK, FARO, HATY, HEMP, HOOK, LISL, LOVE, MABO, MRCH, MVCO, NANT, NAUS, PYFC, and WILD',
                       references='http://maracoos.org/node/146 https://rucool.marine.rutgers.edu/facilities https://rucool.marine.rutgers.edu/data')

sub_dirs = [os.path.join(main_dir, o) for o in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, o)) and o not in avoid]

# load csv file containing the grid
grid = pd.read_csv(hfr_grid, sep=',', header=None, names=['lon', 'lat'], delim_whitespace=True)

# create empty list for finding files
files_grabbed = []

for sub in sub_dirs:
    for ext in types:
        files_grabbed.extend(glob.glob(os.path.join(sub, ext)))


files = sorted(files_grabbed)
for fname in files:
    mat_to_netcdf4.main(grid, fname, save_dir, user_attributes)