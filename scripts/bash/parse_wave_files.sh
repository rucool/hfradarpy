#!/bin/bash
PATH=/bin:/usr/bin;
PATH="/home/codaradm/miniconda3/bin:$PATH"
logdir=$HOME/logs/parse_wave_files
log_file_name=parse_wave_files-$(date --utc +%Y%m%d).log
logfile=$logdir/${log_file_name}

echo ---------------- Start ---------------------- >> $logfile
date >> $logfile

source activate codar_processing 
python /home/codaradm/operational_scripts/codar_processing/codar_processing/methods/waves/parse_waves_to_database.py >> $logfile

#$HOME/logs/parse_wave_files-`date +\%Y\%m\%d`.log
source deactivate

echo ---------------- End ------------------------ >> $logfile
