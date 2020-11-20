#!/bin/bash
PATH=/bin:/usr/bin;
PATH="/home/codaradm/miniconda3/bin:$PATH"
logdir=$HOME/logs/clean_radial_dirs
log_file_name=clean_radial_dirs-$(date --utc +%Y%m%d).log
logfile=$logdir/${log_file_name}

echo ---------------- Start ---------------------- >> $logfile
date >> $logfile

source ~/miniconda3/etc/profile.d/conda.sh
conda activate hfradar
python /home/codaradm/operational_scripts/hfradarpy/hfradar/utilities/clean_radial_dirs.py >> $logfile
conda deactivate

echo ---------------- End ------------------------ >> $logfile
