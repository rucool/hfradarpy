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
% past twelve hours on arctic.
% 
%
% Created by Mike Smith (michaesm@marine.rutgers.edu) on 3/23/2012.
% 2017-05-24 - Changed from camelCase to more standard variable naming
% pattern
% See also radials2totals, CODAR_configuration, CODAR_driver_totals
addpath(genpath('.'))

start_time = datenum(2016,1,1,0,0,0);
end_time = datenum(2016,2,1,0,0,0);
system_type = 2;
radial_directory = '/Users/mikesmith/Documents/projects/bpu/radials/hudson_south/';
save_directory = '/Users/mikesmith/Documents/projects/bpu/totals/domain/';

% Build hourly timestamps
time_steps = end_time:-1/24:start_time;  
conf = codar_configuration(system_type, radial_directory, save_directory);

for x = 1:1:length(time_steps)
    % Process current timestamp
    fprintf(1, '****************************************\n');
    fprintf(1, '  Current time: %s\n',datestr(now));
    fprintf(1, '  Processing data time: %s\n',datestr(time_steps(x),0));

    % Hourly Total Creation
    try
        fprintf(1, 'Starting Totals_driver\n');
        [procFname] = codar_driver_totals(time_steps(x), conf, system_type);
    catch
        fprintf(1, 'Driver_totals failed because: \n');
        res = lasterror;                
        fprintf(1, '%s\n',res.message);
    end
end

fprintf(1, 'Total Creation Finished.\n');
fprintf(1, '---------------------------------------------------------------------------\n');

Display script ending time and elapsed time to file.
end_time = now;
fprintf(1, 'reprocess.m End Time: %s \n', datestr(end_time));
end_time = abs(end_time(floor(datevec(start_time)), floor(datevec(end_time))));
elapsed_minutes = floor(end_time/60);
elapsed_seconds = mod(end_time, 60);
fprintf(1, 'Elapsed time is %s minutes and %s seconds.\n', num2str(elapsed_minutes), num2str(elapsed_seconds));