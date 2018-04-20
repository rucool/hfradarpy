#!/usr/bin/env python
import os
from datetime import timedelta, datetime
import shutil
import glob

# main directory that contains radial files
source_dir = '/home/codaradm/data/radials/'

# File types we want to move
types = ('*.ruv', '*.euv')

# list the subdirectories from source_dir. Only add directories
site_dirs = [os.path.join(source_dir, o) for o in os.listdir(source_dir) if os.path.isdir(os.path.join(source_dir, o))]


# check if the date is greater than 30 days old.
def check_date(file):
    base_date = datetime.utcnow() - timedelta(days=30)
    yy = file.split("_")[2]
    yy = int(yy)
    mm = file.split("_")[3]
    mm = int(mm)
    dd = file.split("_")[4]
    dd = int(dd)
    HH = file.split("_")[5]
    HH = int(HH[:2])
    last_used = datetime(yy, mm, dd, HH, 0, 0, 0)
    if last_used < base_date:
        return 1
    return 0


# iterate through the site directories
for site in site_dirs:
    # create empty list for files_grabbed
    files_grabbed = []

    # grab a list of all files in each site directory that has all extension in the tuple, 'types', created above
    for files in types:
        files_grabbed.extend(glob.glob(os.path.join(site, files)))

    # if files_grabbed is not empty
    if files_grabbed:
        # iterate through each file in the list files_grabbed
        for filename in files_grabbed:
            # grab the filename from the path/filename combo, filena,e
            base_name = os.path.basename(filename)
            yy = base_name.split("_")[2]
            mm = base_name.split("_")[3]

            # create directory name in the format yyyy_mm
            new_dir = os.path.join(site,yy + '_' + mm)

            # try to create the directory. if it exists, pass
            try:
                os.makedirs(new_dir)
            except OSError:
                if os.path.exists(new_dir):
                    pass
                else:
                    raise

            # pass the filename to the check_date function to see if the file is older than 30 days
            found_file = check_date(base_name)

            # if file is older than 30 days
            if found_file == 1:
                # is the file archived already?
                if os.path.isfile(os.path.join(new_dir, base_name)):
                    # if it is, remove the file from the main folder of the site directory
                    print('{} is older than 30 days and already archived in {}. Deleting from source directory.'.format(base_name, new_dir))
                    os.remove(filename)

                # if the file isn't archived already,
                else:
                    # move it to the subdirectory, yyyy_mm
                    print('{} is older than 30 days and is not archived. Moving to {}'.format(base_name,new_dir))
                    shutil.move(filename, new_dir)

            # if the file is newer than 30 days old.
            else:
                # is the file archived already?
                if os.path.isfile(os.path.join(new_dir, base_name)):
                    # if it is, leave it in the main folder of the site directory
                    print('{} already archived.'.format(filename))
                    continue
                # if the file isn't archived already
                else:
                    # copy the file to the yyyy_mm directory, but leave a copy in the main site directory
                    print('Archiving {} to {}'.format(filename, new_dir))
                    shutil.copy(filename, new_dir)