function totals_generate(system_type, database_session, varargin)
%
% Rutgers HFRadar Processing Toolbox
%
% Usage: generate_totals(system_type, database_session, varargin)
%
% Wrapper function for CODAR_driver_totals that queries the radial
% database and returns the most recent timestamps of all active SeaSonde
% sites. If the optional inputs, startime and end_time, are not entered. The
% function grabs the latest timestamp available for radials and checks if the
% radials were used in the production of totals in the past 12 (default)
% hours.
%
% Setting the parameter 'redo' to true will forcefully process totals even 
% if they already exist. It bypasses checking the database to see if the 
% totals were already processed, creates the totals, and inserts a new
% record into the database.
%
% OPTIONAL (name, value) PARAMETERS
%   'start_time' - start time (datenum)
%   'end_time'   - end time (datenum)
%   'hours'     - hours before end time to check for new radials
%   'radial_directory'   - location of main file folder containing sub-folders
%   for each site
%   'save_directory'   - directory to save totals to
%   'redo'      - Reprocess totals and insert into database
%   'useDb'     - use the database? Enter 'true'  or 'false'
%   'pattern_select' - use the following:  1 = best (chosen by operations), 2 = ideal, 3 = measured, 4 = best, ideal, and measured
%
% Created by Mike Smith (michaesm@marine.rutgers.edu) on 3/23/2012.
%
% See also CheckTotals, CODAR_configuration, CODAR_driver_totals 

% Script name
caller = [mfilename '.m'];

% Default values
hours   = 12;
start_time = [];
end_time   = [];
radial_directory = '/home/codaradm/data/';
save_directory = '/home/codaradm/data/';
redo  = false;
region = 1; %1 = maracoos, 2 = caricoos, 13 = antarctica
pattern_select = 1; % 1 = best; 2 = ideal; 3 = measured; 4 = best, ideal, and measured

% Process optional name/value input parameters -------------------------------
    for x = 1:2:length(varargin)

        name  = varargin{x};
        value = varargin{x+1};

        if isempty(value)
            disp([caller ' - Skipping empty option pairing: ' name]);
            continue;
        elseif ~ischar(name)
            disp([caller ' E: Parameter name must be a string.']);
            return;
        end

        switch lower(name)
            case 'redo'
                if ~isequal(numel(value), 1) || ~islogical(value)
                    fprintf(2,...
                        'Value for %s must be a logical.\n',...
                        name);
                    return;
                end
                redo = value;
            case 'radial_directory'
                if ~ischar(value) || ~isdir(value)
                    fprintf(2,...
                        '%s must be a string specifying a valid directory.\n',...
                        name);
                    return;
                end
                radial_directory = value;
            case 'save_directory'
                if ~ischar(value) || ~isdir(value)
                    fprintf(2,...
                        '%s must be a string specifying a valid directory.\n',...
                        name);
                    return;
                end
                save_directory = value;
            case 'start_time'
                if ~isequal(numel(value),1) || ~isnumeric(value)
                    fprintf(2,...
                        'Invalid start time specified\n');
                    return;
                else
                    try
                        start_time = datenum(value);
                    catch ME
                        fprintf(2,...
                        '%s - %s\n',...
                        ME.identifier,...
                        ME.message);
                    return;
                    end
                end
            case 'end_time'
                if ~isequal(numel(value),1) || ~isnumeric(value)
                    fprintf(2,...
                        'Invalid start time specified\n');
                    return;
                else
                    try
                        end_time = datenum(value);
                    catch ME
                        fprintf(2,...
                        '%s - %s\n',...
                        ME.identifier,...
                        ME.message);
                    return;
                    end
                end
            case 'hours'
                if ~isequal(numel(value),1) || ~isnumeric(value) || value < 0
                    fprintf(2,...
                        'Value for %s must be a positive numeric scalar.\n',...
                        name);
                    return;
                end
                hours = value;
            case 'region'
                if ~isequal(numel(value),1) || ~isnumeric(value) || value < 0
                    fprintf(2,...
                        'Value for %s must be a positive numeric scalar.\n',...
                        name);
                     return;
                 end
                 region = value;
             case 'pattern_select'
                if ~isequal(numel(value),1) || ~isnumeric(value) || value < 0
                    fprintf(2,...
                        'Value for %s must be a positive numeric scalar.\n',...
                        name);
                     return;
                 end
                 pattern_select = value;
            otherwise
                fprintf(2,...
                    'Unknown option: %s.\n',...
                    name);
                return;
        end
    end
    
%Build MySQL query statement
db_query = ['select hfrSites.site, hfrSites.id as ',...
    'siteId, hfrLatestRadials.Timestamp as latestRadial, ',...
    'hfrFileTypes.type as filePrefix from hfrSites ',...
    'inner join hfrFileTypes on ',...
    '(hfrSites.radialProcessingType = hfrFileTypes.id) inner join ',...
    'hfrLatestRadials on (hfrSites.id = hfrLatestRadials.siteId) ',...
    'inner join hfrSystemTypes on ',...
    '(hfrSystemTypes.id = hfrSites.type) ',...
    'where(hfrSites.region = ' num2str(region) ' AND ',...
    'hfrSites.active = 1 ',...
    'AND hfrSites.type = ' num2str(system_type) ' AND ',...
    'hfrSites.numsites = 1)'];

result = fetch_data(database_session, db_query);

% Get site codes
 sites = result(:,1)';

% Get pattern types
pattern_type = result(:,4)';

% Convert Time String to datenum
radial_datestrs = datenum(result(:,3), 'yyyy-mm-dd HH:MM:SS.0');


% Get latest radial time and build past n number of hours
if ~isempty(start_time)
    tNow = end_time;
    tThen = start_time;
else
    tNow = max(unique(radial_datestrs));
    tThen = tNow - hours/24;
end

% Build hourly timestamps 
if system_type == 4
    current_time = tThen:1/24:tNow;
else
    current_time = tNow:-1/24:tThen;
end

    for x = 1:1:length(current_time)
        num_times = {};
        for y = 1:length(sites)
            instance = datestr(current_time(x), 'yyyy_mm_dd_HH00.ruv');
            file_name = [radial_directory 'radials/' sites{y} '/' instance(1:7) '/' pattern_type{y} '_' sites{y} '_' instance];
            if exist(file_name, 'file')
                num_times = [num_times; file_name]; %Append to matrix the times that do exist
            end
        end

        if ~redo

            % Build query to get the number of radials used in the creation
            % of totals with this specific time stamp
            db_query = ['select numRadials, region, pattTypes from hfrTotals where ',...
            'RadialTimestamp = ''' datestr(current_time(x), 'yyyy-mm-dd HH:00:00') '''',...
            ' and systemTypeId = ' num2str(system_type) ' AND region = ' num2str(region),...
            ' ORDER BY dateProcessed DESC LIMIT 1'];

            radials_included = fetch_data(database_session, db_query);

            if strcmpi(radials_included, 'No Data') % Record does not exist. Create Totals and insert record into db.
                fprintf(1, 'A record for %s does not exist. Total Creation Starting. \n', datestr(current_time(x), 'yyyy-mm-dd HH:00:00'));
                sum_of_old_radials = 0;
            else % Record exists, but check if the new number of radials is greater than the old number of radials.
                fprintf(1, 'A record for %s exists. Determining if re-processing is needed. \n', datestr(current_time(x), 'yyyy-mm-dd HH:00:00'));
                sum_of_old_radials = radials_included{1};
            end
        else
            fprintf(1, 'Forcefully reprocessing totals. Bypassing database checking.\n');
            sum_of_old_radials = 0;
        end

        sum_of_recent_radials = length(num_times); % Number of Radials that are available for a specific time period

        metaData.region = region;
        metaData.current_time = current_time;
        metaData.database_session = database_session;
        metaData.sum_of_recent_radials = sum_of_recent_radials;
        metaData.sum_of_old_radials = sum_of_old_radials;
        metaData.x = x;
        metaData.pattern_select = pattern_select;

        % Process current timestamp
        fprintf(1, '****************************************\n');
        fprintf(1, '  Current time: %s\n',datestr(now));
        fprintf(1, '  Processing data time: %s\n',datestr(current_time(x),0));
        fprintf(1,  '****************************************\n');

        if sum_of_recent_radials > sum_of_old_radials
            if pattern_select == 1
                createHourlyTotals(system_type, sites, pattern_type, radial_directory, save_directory, metaData);
            elseif pattern_select == 2
                createHourlyTotals(system_type, sites, {'RDLi'}, radial_directory, save_directory, metaData);
            elseif pattern_select == 3
                createHourlyTotals(system_type, sites, {'RDLm'}, radial_directory, save_directory, metaData);
            elseif pattern_select == 4
                db_query = ['select numRadials, region, pattern_types from hfrTotals where ',...
                    'RadialTimestamp = ''' datestr(current_time(x), 'yyyy-mm-dd HH:00:00') '''',...
                    ' and system_typeId = ' num2str(system_type) ' AND region = ' num2str(region),...
                    ' and pattern_types = 1 ORDER BY dateProcessed DESC LIMIT 1'];
                result = fetch_data(database_session, db_query);

                if strcmpi(result, 'No Data') % Record does not exist. Create Totals and insert record into db.
                    fprintf(1, 'A record for %s does not exist. Total Creation Starting. \n', datestr(current_time(x), 'yyyy-mm-dd HH:00:00'));
                    createHourlyTotals(system_type, sites, pattern_type, radial_directory, save_directory, metaData);
                else % Record exists, but check if the new number of radials is greater than the old number of radials.
                    if sum_of_recent_radials > sum_of_old_radials
                        createHourlyTotals(system_type, sites, pattern_type, radial_directory, save_directory, metaData);
                    else
                        fprintf(1, 'not reprocessing\n');
                    end
                end

                db_query = ['select numRadials, region, pattern_types from hfrTotals where ',...
                    'RadialTimestamp = ''' datestr(current_time(x), 'yyyy-mm-dd HH:00:00') '''',...
                    ' and system_typeId = ' num2str(system_type) ' AND region = ' num2str(region),...
                    ' and pattern_types = 2 ORDER BY dateProcessed DESC LIMIT 1'];
                result = fetch_data(database_session, db_query);

                if strcmpi(result, 'No Data') % Record does not exist. Create Totals and insert record into db.
                    fprintf(1, 'A record for %s does not exist. Total Creation Starting. \n', datestr(current_time(x), 'yyyy-mm-dd HH:00:00'));
                    createHourlyTotals(system_type, sites, {'RDLi'}, radial_directory, save_directory, metaData);
                else % Record exists, but check if the new number of radials is greater than the old number of radials.
                    if sum_of_recent_radials > sum_of_old_radials
                        createHourlyTotals(system_type, sites, {'RDLi'}, radial_directory, save_directory, metaData);
                    else
                        fprintf(1, 'not reprocessing\n');
                    end
                end

                db_query = ['select numRadials, region, pattern_types from hfrTotals where ',...
                    'RadialTimestamp = ''' datestr(current_time(x), 'yyyy-mm-dd HH:00:00') '''',...
                    ' and system_typeId = ' num2str(system_type) ' AND region = ' num2str(region),...
                    ' and pattern_types = 3 ORDER BY dateProcessed DESC LIMIT 1'];
                result = fetch_data(database_session, db_query);

                if strcmpi(result, 'No Data') % Record does not exist. Create Totals and insert record into db.
                    fprintf(1, 'A record for %s does not exist. Total Creation Starting. \n', datestr(current_time(x), 'yyyy-mm-dd HH:00:00'));
                    createHourlyTotals(system_type, sites, {'RDLm'}, radial_directory, save_directory, metaData);
                else % Record exists, but check if the new number of radials is greater than the old number of radials.
                    if sum_of_recent_radials > sum_of_old_radials
                        createHourlyTotals(system_type, sites, {'RDLm'}, radial_directory, save_directory, metaData);
                    else
                        fprintf(1, 'not reprocessing\n');
                    end                    
                end
            end
        elseif sum_of_recent_radials < 2
            fprintf(1, 'Not enough radials to generate totals. Continuing to next timestep. \n');
        else
            fprintf(1, 'Record exists and timestamp does NOT need to be re-processed. \n');
        end
        clear num_times
    end
end
