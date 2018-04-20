tic
source_dir = '/home/codaradm/data/totals_25MHz/ascii/';
javaaddpath('/home/kerfoot/bin/matlab/toolboxes/toolsUI-2.2.16.jar');
addpath(genpath('/home/kerfoot/bin/matlab/toolboxes/snctools'));


opt = 'new'; % activitate John Wilkin changes for scanning for last NumDaysOld of files
NumDaysOld = 10;

if strcmp(opt,'new')
  source_files = [source_dir 'OI_MARASR_*'];
else
  c=fix(clock);
end

%% This file is used for the MARACOOS extended grid totals /home/lemus/MARACOOS/scripts_nc/template_maracoos_compress.nc
output_base_dir = '/home/codaradm/data/totals_25MHz/netcdf';
template_filename = '/home/codaradm/data/grid_files/nc_templates/template_maracoos25.nc';
%grid_file = '/home/codaradm/data/grid_files/MARACOOS_1km_SRGrid.txt';

%% This file is used for the Josh totals /home/lemus/MARACOOS/scripts_nc/template_macoora_inertial.nc
% Right now, just take the last 60 in the directory.
codar_files = dir(source_files);

% If this variable is non-empty at the end, then we use it to construct a new "latest"
% dataset.
output_netcdf_filename = [];

% Grid to 2D lat/lon
%[P] = maracoos1km_to_grid(grid_file);
load('/home/codaradm/data/grid_files/nc_templates/netcdf_2d_grids/MARACOOS_1km_2dgrid.mat');

lat_size = nc_getdiminfo( template_filename, 'lat' );
lon_size = nc_getdiminfo( template_filename, 'lon' );
lat_size = lat_size.Length;
lon_size = lon_size.Length;

for j = 1:length(codar_files)

	% Don't process a directory, of course.
	if codar_files(j).isdir
		continue;
	end

	% Don't process if the file is empty.
	if codar_files(j).bytes == 0
		continue;
    end
    
    if strcmp(opt,'new')
        % Don't process if the file is more than NumDaysOld
        if datenum(now)-datenum(codar_files(j).date) > NumDaysOld
            continue;
        end
    end
  
	% determine the timestamp
	[d,count] = sscanf ( codar_files(j).name, 'OI_MARASR_%d_%d_%d_%d' );
	if ( count ~= 4 )
		msg = sprintf ( '%s:  could not parse timestamp from %s\n', mfilename, codar_files(j).name );
		error ( msg );
	end

	gmt_time = datenum ( d(1), d(2), d(3), d(4)/100, 0, 0 );
    
	% make sure that the output directory exists.
	[success,message,messageID] = mkdir(output_base_dir);

    output_netcdf_filename = sprintf ( '%s/%s.totals.nc', output_base_dir, codar_files(j).name );
    
	% Does the file already exist?
    if exist ( output_netcdf_filename, 'file' )
        continue;
    else
        [success,message,messageID] = copyfile(template_filename, output_netcdf_filename);
    end
    
    % Add the history attribute
    nc_attput ( output_netcdf_filename, nc_global, 'history', 'Hourly codar combined into one monthly file. See source attribute' );
    
    nc_attput ( output_netcdf_filename, nc_global, 'source', source_files );

	fprintf ( 1, 'Processing %s into %s\n', codar_files(j).name, output_netcdf_filename );

	text_file = sprintf ( '%s%s', source_dir, codar_files(j).name );

	process_codar_totals_into_netcdf ( text_file, output_netcdf_filename, gmt_time, lat_size, lon_size, P);
toc
end
