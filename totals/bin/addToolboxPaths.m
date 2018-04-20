% Setup paths
javaaddpath('/Users/hroarty/Documents/MATLAB/HFR_RealTime/totals_toolbox/bin/mysql-connector-java-5.1.6-bin.jar');
%addpath(genpath('/Users/hroarty/Documents/MATLAB/HFR_Progs-2_1_3beta/matlab/general/'))
add_subdirectories_to_path('/Users/hroarty/Documents/MATLAB/HFR_Progs-2_1_3beta/matlab/',{'CVS','private','@'});

% Add the snctools jar file
disp( 'Adding NetCDF libraries...' );
lasterr('');
try
    javaaddpath('/Users/hroarty/Documents/MATLAB/netcdfAll-4.2.jar');
    addpath(genpath('/Users/hroarty/Documents/MATLAB/matlabToolboxes'));
catch
    disp(['NetCDF-4 Error: ' lasterr]);
end