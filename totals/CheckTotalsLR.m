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
javaaddpath('./bin/mysql-connector-java-5.1.6-bin.jar');

% System Type: (long, mid, or standard) range?
system_type = 'long';

start_time = now; 
fprintf('CheckTotals.m Start Time: %s\n', datestr(start_time)); % Display start time of script

database_session = configs_database; % connect to the totals/radial database

% Convert from string to assigned number
range_lookup = struct(...
    'standard', 4,...
    'mid', 3,...
    'long', 2);

vars = fieldnames(range_lookup); % get field names of each variable

for varInd = 1:length(vars)
    variable = vars{varInd};
    
      if strcmpi(system_type, variable)
          system_type = range_lookup.(variable); % get the corresponding #
      else
          continue
      end
end

% Build database execution statement
database_statement = ['select count(distinct(hfrSites.type)) from hfrSites where hfrSites.type = ' num2str(system_type)];

result = fetch_data(database_session, database_statement); % see if hfr system type is in database

% See if the array "result" contains any data. 
if result{1} == 0     %If not, return from function
    fprintf(1, 'System Type %s (%s) does not exist. Exiting.\n', num2str(system_type), system_type);
    return
else    % If it does, continue on past this if statement.
    fprintf(1, 'System Type %s exists. Proceeding to total creation\n', num2str(system_type));
end    
fprintf(1, '---------------------------------------------------------------------------\n');

% Proceed to total generation by calling generate_totals. This also updates
% the codar database
totals_generate(system_type, database_session,...
    'start_time', datenum(2017,8,1,0,0,0),...
    'end_time', datenum(2017,8,1,0,0,0),...
    'pattern_select', 1,...
    'region', 7,...
    'radial_directory', '/Volumes/codaradm/data/',...
    'save_directory', '/Users/mikesmith/Documents/',...
    'redo',true);

% Close MySQL Connection
close(database_session);
fprintf(1, 'Total Creation Finished.\n');
fprintf(1, '---------------------------------------------------------------------------\n');

% Display script ending time and elapsed time to file.
end_time = now;
fprintf(1, 'CheckTotals.m End Time: %s \n', datestr(end_time));
end_time = abs(end_time(floor(datevec(start_time)), floor(datevec(end_time))));
elapsed_minutes = floor(end_time/60);
elapsed_seconds = mod(end_time, 60);
fprintf(1, 'Elapsed time is %s minutes and %s seconds.\n', num2str(elapsed_minutes), num2str(elapsed_seconds));
