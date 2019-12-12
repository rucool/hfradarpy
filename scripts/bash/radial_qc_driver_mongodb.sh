#!/bin/bash
PATH=/bin:/usr/bin;
PATH="/home/codaradm/miniconda3/bin:$PATH"
logdir=$HOME/logs/radial_qc_driver
log_file_name=radial_qc_driver-$(date --utc +%Y%m%d).log
logfile=$logdir/${log_file_name}

echo ---------------- Start ---------------------- >> $logfile
date >> $logfile

source ~/miniconda3/etc/profile.d/conda.sh
conda activate codar_processing
python /home/codaradm/operational_scripts/codar_processing/scripts/realtime/radial_qc_driver_mongodb.py >> $logfile
conda deactivate

echo ---------------- End ------------------------ >> $logfile
