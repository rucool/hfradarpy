import glob
import os
import pandas as pd
from functions import mat_to_netcdf4

main_dir = '/Volumes/home/codaradm/data/totals/pldp/25MHz/0.5km/oi/mat/'
save_dir = '/Users/mikesmith/Documents/data/codar/totals/nc'
hfr_grid = '../totals/grid_files/PLDP_0.5km_grid.txt'

avoid = ('ideal', 'measured') # avoid these subfolders
types = ['*.mat']  # find only files with these types

sub_dirs = [os.path.join(main_dir, o) for o in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, o)) and o not in avoid]

# load csv file containing the grid
grid = pd.read_csv(hfr_grid, sep=',', header=None, names=['lon', 'lat'], delim_whitespace=True)

# create empty list for finding files
files_grabbed = []

for sub in sub_dirs:
    for ext in types:
        files_grabbed.extend(glob.glob(os.path.join(sub, ext)))

user_attributes = dict(title='Palmer Deep Antarctica 0.5 km Sea Surface Currents',
                       naming_authority='edu.rutgers.marine.rucool',
                       comment='Network maintained by University of Alaska - Fairbanks and Rutgers University. For oi_* global attribute explanations, see references attribute',
                       acknowledgment='This data is provided by the Rutger University, Center for Ocean Observing Leadership. Funding is provided by the U.S. National Science Foundation (NSF).',
                       creator_name='Michael Smith',
                       creator_email='michaesm@marine.rutgers.edu',
                       creator_url='rucool.marine.rutgers.edu',
                       institution='Center for Ocean Observing and Leadership, Department of Marine & Coastal Sciences, Rutgers University',
                       project='Project Converge - High Frequency Radar Sea Surface Current Mapping',
                       sea_name='Mid-Atlantic Bight',
                       creator_type='person',
                       creator_institution='Rutgers University',
                       contributor_name='Josh Kohut, Peter Winsor, Hank Statscewich, Hugh Roarty, Ethan Handel, Michael Smith, Laura Nazzaro',
                       contributor_role='Principal Investigator, Principal Investigator, Hardware Maintenance, Manager, Maintenance, Data Manager, Data Manager',
                       platform='Palmer Deep, Antarctica HF Radar 25MHz Network',
                       instrument='Network includes CODAR sites AMAG, ASSA, BLCK, BRIG, CEDR, CORE, DUCK, FARO, HATY, HEMP, HOOK, LISL, LOVE, MABO, MRCH, MVCO, NANT, NAUS, PYFC, and WILD',
                       references='http://maracoos.org/node/146 https://rucool.marine.rutgers.edu/facilities https://rucool.marine.rutgers.edu/data')

files = sorted(files_grabbed)
for fname in files:
    mat_to_netcdf4.main(grid, fname, save_dir, user_attributes)