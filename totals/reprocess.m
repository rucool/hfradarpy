% Rutgers HFRadar Processing Toolbox
%
% CheckTotals.m
%
% This is the main input script for total generation. This script works
% for all monostatic SeaSonde system types. Bistatic will be implemented in
% the future.
%
% This script connects to the MySQL database to see if the system type
% (long, mid, standard) exists. If the system type exists, this script
% calls radials2totals.m, which by default, determines the times that 
% totals need to be generated. 
%
% The only inputs you will most likely need to change are "system_type" on
% line 28 and the 'pattType' input for radials2totals.m. You may also 
% change the optional inputs for radials2totals.m. In order to see the 
% availabe options type "help radials2totals.m" If the argument inputs for 
% radials2totals.m are not entered by the user, the script will use the 
% default settings which are preset in order to generate totals for the 
% past twelce hours on arctic.
% 
%
% Created by Mike Smith (michaesm@marine.rutgers.edu) on 3/23/2012.
% 2017-05-24 - Changed from camelCase to more standard variable naming
% pattern
% See also radials2totals, CODAR_configuration, CODAR_driver_totals
addpath(genpath('.'))

% Proceed to total generation by calling generate_totals. This also updates
% the codar database
totals_reprocess(2,...
    'start_time', datenum(2016,1,23,16,0,0),...
    'end_time', datenum(2016,1,23,16,0,0),...
    'region', 7,...
    'radial_directory', '/Users/mikesmith/Documents/projects/bpu/radials/fairways_south/',...
    'save_directory', '/Users/mikesmith/Documents/projects/bpu/totals/fairways_south/');

% Close MySQL Connection
% close(database_session);
fprintf(1, 'Total Creation Finished.\n');
fprintf(1, '---------------------------------------------------------------------------\n');

% Display script ending time and elapsed time to file.
% end_time = now;
% fprintf(1, 'reprocess.m End Time: %s \n', datestr(end_time));
% end_time = abs(end_time(floor(datevec(start_time)), floor(datevec(end_time))));
% elapsed_minutes = floor(end_time/60);
% elapsed_seconds = mod(end_time, 60);
% fprintf(1, 'Elapsed time is %s minutes and %s seconds.\n', num2str(elapsed_minutes), num2str(elapsed_seconds));
