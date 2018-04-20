function createHourlyTotals(system_type, sites, pattern_type, radial_directory, save_directory, metaData)

    if length(pattern_type) > 1
            pattern_select = 1;
    else
        if strcmpi(pattern_type, 'RDLi')
            pattern_select = 2;
        elseif strcmpi(pattern_type, 'RDLm')
            pattern_select = 3;
        end
    end
        % Build configuration file for use in CODAR_driver_totals
        conf = CODAR_configuration(system_type, sites, pattern_type, radial_directory, save_directory, metaData.region);

        % Hourly Total Creation
        try
            
            fprintf(1,  'Making Hourly Totals. \n');
            [procFname] = CODAR_driver_totals(metaData.current_time(metaData.x), conf, system_type);
            ingestCodar(system_type, procFname, metaData.region, pattern_select) 
            
            if exist(procFname, 'file')            % Check if the file was created. If it was, insert the record.
                fprintf(1, 'The file %s exists. Inserting record into database. \n', procFname);
                fprintf(1, '%s radials were used in the generation of the IOOS region #%s, %s totals. \n',...
                    num2str(metaData.sum_of_recent_radials), num2str(metaData.region), datestr(metaData.current_time(metaData.x)));
                fastinsert(metaData.database_session, 'hfrTotals',...
                    {'systemTypeId', 'radialTimestamp', 'numRadials', 'region', 'pattTypes'},...
                    {systemType, datestr(metaData.current_time(metaData.x), 'yyyy-mm-dd HH:MM:00'),...
                    metaData.sum_of_recent_radials, metaData.region, pattern_select});
            end
            
        catch
            fprintf(1, 'Total Creation failed because: \n');
            res=lasterror;
            fprintf(1, '%s\n',res.message);
        end    

end