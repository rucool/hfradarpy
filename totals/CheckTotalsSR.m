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
% The only inputs you will most likely need to change are "SystemString" on
% line 28 and the 'pattType' input for radials2totals.m. You may also 
% change the optional inputs for radials2totals.m. In order to see the 
% availabe options type "help radials2totals.m" If the argument inputs for 
% radials2totals.m are not entered by the user, the script will use the 
% default settings which are preset in order to generate totals for the 
% past twelce hours on arctic.
% 
%
% Created by Mike Smith (michaesm@marine.rutgers.edu) on 3/23/2012.
%
% See also radials2totals, CODAR_configuration, CODAR_driver_totals

% addpath(genpath('/Users/codar/Documents/MATLAB/totals_toolbox/'))%
% removed hjr 2015.06.02


% System Type: (long, mid, or standard) range?
SystemString = 'mid';

sTime = now; fprintf('CheckTotals.m Start Time: %s\n', datestr(sTime)); % Display start time of script

addToolboxPaths; % Setup toolbox paths

dbConn = connect2TotalsDB;

% Convert from string to assigned number
rangeVars = struct('standard', 4,...
    'mid', 3,...
    'long', 2);

vars = fieldnames(rangeVars);

for varInd = 1:length(vars)
    variable = vars{varInd};
    
      if strcmpi(SystemString, variable)
          systemType = rangeVars.(variable);
      else
          continue
      end
end

% Build database execution statement
dbStatement = ['select count(distinct(hfrSites.type)) from hfrSites where hfrSites.type = ' num2str(systemType)];

result = fetchDatafromDB(dbConn, dbStatement);

% See if the array "result" contains any data. 
if result{1} == 0;     %If not, return from function
    fprintf(1, 'System Type %s (%s) does not exist. Exiting.\n', num2str(systemType), SystemString);
    return
else    % If it does, continue on past this if statement.
    fprintf(1, 'System Type %s exists. Proceeding to total creation\n', num2str(systemType));
end    

fprintf(1, '---------------------------------------------------------------------------\n');
disp('                     ');

% Proceed to total generation by calling radials2totals
radials2totals(systemType, dbConn,...
    'startTime', datenum(2015,6,3,0,0,0),...
    'endTime', datenum(2015,6,4,0,0,0),...
    'pattSelect', 3,...
    'region', 7,...
    'baseDir', '/Volumes/codaradm/data/',...
    'saveDir', '/Users/hroarty/data/realtime/',...
    'reDo',true,...
    'useDb',false);



% Close MySQL Connection
close(dbConn);
fprintf(1, 'Total Creation Finished.\n');
fprintf(1, '---------------------------------------------------------------------------\n');
fprintf(1, '\n');

% Display script ending time and elapsed time to file.
eTime = now;
fprintf(1, 'CheckTotals.m End Time: %s \n', datestr(eTime));
eTime = abs(etime(floor(datevec(sTime)), floor(datevec(eTime))));
eMins = floor(eTime/60);
eSecs = mod(eTime, 60);
fprintf(1, 'Elapsed time is %s minutes and %s seconds.\n', num2str(eMins), num2str(eSecs));
