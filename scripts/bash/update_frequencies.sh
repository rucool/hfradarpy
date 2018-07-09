#!/bin/bash
PATH=/bin:/usr/bin;
PATH="/home/codaradm/miniconda3/bin:$PATH"
logdir=$HOME/logs/update_frequencies
log_file_name=update_frequencies-$(date --utc +%Y%m%d).log
logfile=$logdir/${log_file_name}

echo ---------------- Start ---------------------- >> $logfile
date >> $logfile

source activate codar_processing
python /home/codaradm/operational_scripts/codar_processing/utilities/update_frequencies.py >> $logfile
source deactivate

echo ---------------- End ------------------------ >> $logfile
