#!/bin/bash
PATH=/bin:/usr/bin;
PATH="/home/codaradm/miniconda3/bin:$PATH"
logdir=$HOME/logs/clean_radial_dirs
log_file_name=clean_radial_dirs-$(date --utc +%Y%m%d).log
logfile=$logdir/${log_file_name}

echo ---------------- Start ---------------------- >> $logfile
date >> $logfile

source activate codar_processing
python /home/codaradm/operational_scripts/codar_processing/utilities/clean_radial_dirs.py >> $logfile
source deactivate

echo ---------------- End ------------------------ >> $logfile
