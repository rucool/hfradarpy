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

start_time = datenum(2013, 2, 28, 0, 0, 0);
end_time = datenum(2013,3,1,14,0,0);
%end_time = datenum(2009, 12, 19, 2, 0, 0);
%start_time = datenum(2012,7,9,0,0,0);
%end_time = datenum(2012,8,11,0,0,0)
dtime = [end_time:-1/24:start_time];
strtime = datestr(dtime, 'dd-mmm-yyyy HH:00:00');
dtime = datenum(strtime);

for h = 1:length(dtime)
    file_name = ['/home/codaradm/data/totals/maracoos/oi/ascii/25MHz/OI_MARASR_' datestr(dtime(h), 'yyyy_mm_dd_HH00')];
    %file_name = ['/home/codaradm/data/totals/maracoos/oi/nc/25MHz/RU_13MHz_' datestr(dtime(h), 'yyyy_mm_dd_HH00')];
    %if exist(file_name2, 'file')
    %    continue
    %else
    ingestCodarEdit(4, file_name, 7)
    %end
end
