function createHourlyTotals(systemType, Sites, pattType, baseDir, saveDir, metaData)
    
    if length(pattType) > 1
            pattSelect = 1;
    else
        if strcmpi(pattType, 'RDLi')
            pattSelect = 2;
        elseif strcmpi(pattType, 'RDLm')
            pattSelect = 3;
        end
    end
        % Build configuration file for use in CODAR_driver_totals
        conf = CODAR_configuration(systemType, Sites, pattType, baseDir, saveDir, metaData.region);
                
        [m, n] = size(conf.Totals.GridFile);
        
        if n == 1
            % Hourly Total Creation
                try

                    fprintf(1,  'Making Hourly Totals. \n');
                    [procFname] = CODAR_driver_totals(metaData.tCurrent(metaData.x), conf, systemType);
                    ncFile = ingestCodar(systemType, procFname, metaData.region, pattSelect);

                    %%
                    % Run Lauras detiding, filtering, diverence, and
                    % vorticity code
                    ind1 = find(ncFile == '/');
                    ind2 = find(ncFile == '.');
                    curName = ncFile(ind1(end)+1:ind2(1)-1);
                    ind3 = find(curName == '_');              
                    curDate = datenum(curName(ind3(2)+1:end), 'yyyy_mm_dd_HHMM');
                    earlyDate = datenum(curDate-6/24);
                    
                    curDate = datestr(curDate, 'yyyy_mm_dd_HHMM');
                    earlyDate = datestr(earlyDate, 'yyyy_mm_dd_HHMM');
                    
                    earlierInputFile = strrep(ncFile, curDate, earlyDate);
                    disp(['Earlier file is ' earlierInputFile])
                    
%                     earlierinputfile=[inputdir 'RU_13MHz_2014_06_10_0300.totals.nc'];
                    coefficientFile = [ncFile(1:ind1(end)) 'detidingcoefficients/detidingcoefficients_current.nc'];
                    % can be done as files come in
                    [~]=detideCurrents(ncFile, coefficientFile, 'addtide', true, 'adddetided', true);

                    % use on earlier file (maybe 6 h previous?)
                    % timestep option will have to switch to .5
                    [~]=filterCurrents(earlierInputFile,'timestep',1,'addlowpass',true,'addhighpass',true,'currenttype','X_detided');

                    % get divergence from raw u and v (can be done as files come in on all but
                    % filtered)
                    [~]=divergence('netcdf','ncfile',ncFile,'adddivergence',true,'addvorticity',true,'currenttype','X');

                    % get 7-day previous raw divergence trend
                    [~]=trendData(ncFile,'div',.1,'addtrend',true,'overwrite',true);

                    % get 7-day previous raw vorticity trend
                    [~]=trendData(ncFile,'vor',.02,'addtrend',true,'overwrite',true);
                    %%
                    
                    if exist(procFname, 'file')            % Check if the file was created. If it was, insert the record.
                        fprintf(1, 'The file %s exists. Inserting record into database. \n', procFname);
                        fprintf(1, '%s radials were used in the generation of the IOOS region #%s, %s totals. \n',...
                            num2str(metaData.newNumRadials), num2str(metaData.region), datestr(metaData.tCurrent(metaData.x)));
                        fastinsert(metaData.dbConn, 'hfrTotals',...
                            {'systemTypeId', 'radialTimestamp', 'numRadials', 'region', 'pattTypes'},...
                            {systemType, datestr(metaData.tCurrent(metaData.x), 'yyyy-mm-dd HH:MM:00'),...
                            metaData.newNumRadials, metaData.region, pattSelect});
                    end

                catch
                    fprintf(1, 'Total Creation failed because: \n');
                    res=lasterror;
                    fprintf(1, '%s\n',res.message);
                end
        else
            for x = 1:length(conf.Totals.GridFile)

                conf_temp = conf;
                conf_temp.Totals.GridFile = conf.Totals.GridFile{x};
                conf_temp.OI.BaseDir = conf.OI.BaseDir{x};
                conf_temp.OI.AsciiDir = conf.OI.AsciiDir{x};

                % Hourly Total Creation
                try

                    fprintf(1,  'Making Hourly Totals. \n');
                    [procFname] = CODAR_driver_totals(metaData.tCurrent(metaData.x), conf_temp, systemType);
                    ncFile = ingestCodar(systemType, procFname, metaData.region, pattSelect);

                                        %%
                    % Run Lauras detiding, filtering, diverence, and
                    % vorticity code
                    ind1 = find(ncFile == '/');
                    ind2 = find(ncFile == '.');
                    curName = ncFile(ind1(end)+1:ind2(1)-1);
                    ind3 = find(curName == '_');              
                    curDate = datenum(curName(ind3(2)+1:end), 'yyyy_mm_dd_HHMM');
                    earlyDate = datenum(curDate-6/24);
                    
                    curDate = datestr(curDate, 'yyyy_mm_dd_HHMM');
                    earlyDate = datestr(earlyDate, 'yyyy_mm_dd_HHMM');
                    
                    earlierInputFile = strrep(ncFile, curDate, earlyDate);
                    disp(['Earlier file is ' earlierInputFile])
                    
%                     earlierinputfile=[inputdir 'RU_13MHz_2014_06_10_0300.totals.nc'];
                    coefficientFile = [ncFile(1:ind1(end)) 'detidingcoefficients/detidingcoefficients_current.nc'];
                    % can be done as files come in
                    [~]=detideCurrents(ncFile, coefficientFile, 'addtide', true, 'adddetided', true);

                    % use on earlier file (maybe 6 h previous?)
                    % timestep option will have to switch to .5
                    [~]=filterCurrents(earlierInputFile,'timestep',1,'addlowpass',true,'addhighpass',true,'currenttype','X_detided');

                    % get divergence from raw u and v (can be done as files come in on all but
                    % filtered)
                    [~]=divergence('netcdf','ncfile',ncFile,'adddivergence',true,'addvorticity',true,'currenttype','X');

                    % get 7-day previous raw divergence trend
                    [~]=trendData(ncFile,'div',.1,'addtrend',true,'overwrite',true);

                    % get 7-day previous raw vorticity trend
                    [~]=trendData(ncFile,'vor',.02,'addtrend',true,'overwrite',true);
                    %%
                    
                    if exist(procFname, 'file')            % Check if the file was created. If it was, insert the record.
                        fprintf(1, 'The file %s exists. Inserting record into database. \n', procFname);
                        fprintf(1, '%s radials were used in the generation of the IOOS region #%s, %s totals. \n',...
                            num2str(metaData.newNumRadials), num2str(metaData.region), datestr(metaData.tCurrent(metaData.x)));
                        fastinsert(metaData.dbConn, 'hfrTotals',...
                            {'systemTypeId', 'radialTimestamp', 'numRadials', 'region', 'pattTypes'},...
                            {systemType, datestr(metaData.tCurrent(metaData.x), 'yyyy-mm-dd HH:MM:00'),...
                            metaData.newNumRadials, metaData.region, pattSelect});
                    end

                catch
                    fprintf(1, 'Total Creation failed because: \n');
                    res=lasterror;
                    fprintf(1, '%s\n',res.message);
                end
            end
        end

end