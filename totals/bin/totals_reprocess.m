function totals_reprocess(systemType, varargin)
%
% Rutgers HFRadar Processing Toolbox
%
% Usage: totals_reprocess(systemType, dbConn, varargin)
%
% Wrapper function for CODAR_driver_totals 
% 
% This script does not require database usage. It bypass the database 
% querying portion of the script and will process any times between 
% start_time and end_time. If these inputs are not given, the script will 
% take the current time in EDT and process the past twelve hours. This is 
% best for reprocessing totals in which we do not plan on displaying to the public.
%
% OPTIONAL (name, value) PARAMETERS
%   'start_time' - start time (datenum)
%   'end_time'   - end time (datenum)
%   'hours'     - hours before end time to check for new radials
%   'radial_directory'   - location of main file folder containing sub-folders
%   for each site
%   'pattern_type'  - RDLi (Ideal) or RDLm (Measured) pattern radials?
%   'save_directory'   - directory to save totals to
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
radial_directory = '/Users/hroarty/data/realtime/';
pattern_type = 'RDLi';
save_directory = '/Users/hroarty/data/realtime/';
region = 1; %1 = maracoos, 2 = caricoos, 13 = antarctica

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
            case 'pattern_type'
                if ~ischar(value)
                    fprintf(2,...
                        'Pattern type (%s) must be a valid pattern type string.\n',...
                        name);
                    return;
                end
                pattern_type = value;
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
            otherwise
                fprintf(2,...
                    'Unknown option: %s.\n',...
                    name);
                return;
        end
    end
    
    %Do not use the database.. forcefully reprocess the files...
        fprintf(1, 'Bypassing database completely. \n');
        if ~isempty(start_time)
            tNow = end_time;
            tThen = start_time;
        else
            tNow = floor(now*24)/24;
            tThen = tNow - hours/24;
        end

        tCurrent = [tNow:-1/24:tThen]; % Build hourly timestamps 

        for x = 1:1:length(tCurrent)
            nTime = {};
            
            % Build configuration file for use in CODAR_driver_totals
            % Leave the second input, Sites, as an empty array with no data
            
            %conf = CODAR_configuration(systemType, [], [], radial_directory, save_directory, region);
            metaData.region = region;
            conf = CODAR_configuration_HJR_laptop(systemType, [], [], radial_directory, save_directory, metaData.region);
            for y = 1:length(conf.Radials.Sites)
                file_name = [conf.Radials.radial_directory conf.Radials.Sites{y} '/' datestr(tCurrent(x), 'yyyy_mm') '/' conf.Radials.Types{y} '_' conf.Radials.Sites{y} '_' datestr(tCurrent(x), 'yyyy_mm_dd_HH00.ruv')];
                if exist(file_name, 'file');
                    nTime = [nTime; file_name]; %Append to matrix the times that do exist
                end
            end

            % Process current timestamp
            fprintf(1, '****************************************\n');
            fprintf(1, '  Current time: %s\n',datestr(now));
            fprintf(1, '  Processing data time: %s\n',datestr(tCurrent(x),0));
            
            % Hourly Total Creation
            try
                fprintf(1, 'Starting Totals_driver\n');
                [procFname] = CODAR_driver_totals_QC(tCurrent(x),conf, systemType);
                %ingestCodar(systemType, procFname, region,pattSelect)
            catch
                fprintf(1, 'Driver_totals failed because: \n');
                res=lasterror;                
                fprintf(1, '%s\n',res.message);
            end
  
            clear nTime
        end
    end
end
