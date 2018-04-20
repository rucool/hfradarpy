addpath(genpath('/home/codaradm/operational_scripts/totals_toolbox/'))

% Add the snctools jar file
disp( 'Adding NetCDF libraries...' );
lasterr('');
try
    javaaddpath('/home/coolgroup/matlabToolboxes/netcdfAll-4.2.jar');
    addpath(genpath('/home/coolgroup/matlabToolboxes')); 
catch
    disp(['NetCDF-4 Error: ' lasterr]);
end

start_time = datenum(2012, 4, 2, 0, 0, 0);
end_time = datenum(2012,4,3,0,0,0);
%end_time = datenum(2009, 12, 19, 2, 0, 0);
%start_time = datenum(2012,7,9,0,0,0);
%end_time = datenum(2012,8,11,0,0,0)
dtime = [end_time:-1/24:start_time];
strtime = datestr(dtime, 'dd-mmm-yyyy HH:00:00');
dtime = datenum(strtime);

for h = 1:length(dtime)
    file_name = ['/home/michaesm/Codar/reprocessed/pre09LR/totals/maracoos/oi/ascii/5MHz/OI_MARA_' datestr(dtime(h), 'yyyy_mm_dd_HH00')];
    %file_name2 = ['/home/codaradm/data/totals/maracoos/oi/nc/13MHz/RU_13MHz_' datestr(dtime(h), 'yyyy_mm_dd_HH00')];
    %if exist(file_name2, 'file')
    %    continue
    %else
    ingestCodarCustom(2, file_name, 7)
    %end
end
