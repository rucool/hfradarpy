#!/bin/bash
PATH=/bin:/usr/bin;
PATH="/home/codaradm/miniconda3/bin:$PATH"
logdir=$HOME/logs/realtime_waves_to_netcdf
log_file_name=realtime_waves_to_netcdf-$(date --utc +%Y%m%d).log
logfile=$logdir/${log_file_name}

echo ---------------- Start ---------------------- >> $logfile
date >> $logfile

source activate codar_processing
python /home/codaradm/operational_scripts/codar_processing/scripts/realtime_waves_to_netcdf.py >> $logfile
#$HOME/logs/parse_wave_files-`date +\%Y\%m\%d`.log
source deactivate

echo ---------------- End ------------------------ >> $logfile
