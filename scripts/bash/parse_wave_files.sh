#!/bin/bash
PATH=/bin:/usr/bin;
PATH="/home/codaradm/anaconda2/bin:$PATH"
logdir=$HOME/logs/parse_wave_files
log_file_name=parse_wave_files-$(date --utc +%Y%m%d).log
logfile=$logdir/${log_file_name}

source activate hfr_processing 
python /home/codaradm/operational_scripts/hfr_processing/parse_wave_files.py > $logfile
#$HOME/logs/parse_wave_files-`date +\%Y\%m\%d`.log
source deactivate
