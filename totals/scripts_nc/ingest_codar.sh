#!/bin/bash
. ~/.bashrc
#---------------------------------------
# Rutgers Codar Processing Script
# # Main call script for BPU matlab processing to create Netcdf files
# # Edited 10/13/11 by Rivera
# #---------------------------------------
 export MATLAB_SHELL="/bin/sh"

# Setup logfile
 LOGFILE=`date "+/home/codaradm/logs/totals_25MHz_nc/MARASR_NETCDF_%Y_%m_%d_%H%M.txt"`

# Start script
 date >> $LOGFILE

# # Run the Master Matlab file
 cd /home/codaradm/operational_scripts/totals_toolbox/scripts_nc
 /usr/local/MATLAB/R2011a/bin/matlab -nodisplay </home/codaradm/operational_scripts/totals_toolbox/scripts_nc/ingest_codar.m >> $LOGFILE 

 date >> $LOGFILE

 exit 0

