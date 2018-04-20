function ingestCodarEdit(systemType, asciiName, region)

% Add the snctools jar file
%disp( 'Adding NetCDF libraries...' );
%lasterr('');
%try
%    javaaddpath('/home/coolgroup/matlabToolboxes/netcdfAll-4.2.jar');
%    addpath(genpath('/home/coolgroup/matlabToolboxes'));
%catch
%    disp(['NetCDF-4 Error: ' lasterr]);
%end


switch num2str(systemType)
    
    case '2'
    load('/home/codaradm/data/grid_files/nc_templates/netcdf_2d_grids/OI_6km_ExtendedGrid.mat');
    output_base_dir = '/home/michaesm/Codar/reprocessed/pre09LR/totals/maracoos/oi/nc';
    %     output_base_dir = '/Users/michaesm/Documents/MATLAB/totals_toolbox/5MHz_nc4';
    template_filename = '/home/codaradm/data/grid_files/nc_templates/OI_6km_nc4_template.nc';
    title = 'MARACOOS 5MHz Sea Surface Currents';
    oldHead = 'OI_MARA';
    newHead = 'RU_5MHz';
    summary = ['Optimally Interpolated Total Vectors calculated by HFRProgs toolbox using MATLAB. Mercator lat/lon projection using 6km grid.'];
    keyWd = '6km';

    case '3'
        switch num2str(region)
            case '7'
                load('/home/codaradm/data/grid_files/nc_templates/netcdf_2d_grids/OI_2km_Grid.mat');
                output_base_dir = '/home/michaesm/Codar/reprocessed/totals/ideal/totals/maracoos/oi/nc';
                ncDir = '/home/om/dods-data/thredds/cool/codar/totals/13MHz_2km_realtime';
            %     output_base_dir = '/Users/michaesm/Documents/MATLAB/totals_toolbox/5MHz_nc4';
                template_filename = '/home/codaradm/data/grid_files/nc_templates/OI_2km_nc4_template.nc';
                title = 'MARACOOS 13MHz Sea Surface Currents';
                oldHead = 'OI_BPU';
                newHead = 'RU_13MHz';
                summary = ['Optimally Interpolated Total Vectors calculated by HFRProgs toolbox using MATLAB. Mercator lat/lon projection using 2km grid.'];
                keyWd = '2km';
            case '3'
                load('/home/codaradm/data/grid_files/nc_templates/netcdf_2d_grids/PR_2km_Grid.mat');
                output_base_dir = '/home/codaradm/data/totals/caracoos/oi/nc/13MHz';
%               ncDir = '/home/om/dods-data/thredds/cool/codar/totals/13MHz_2km_realtime';
                ncDir = '/home/om/dods-data/thredds/cool/codar/totals/caracoos/realtime/';
                template_filename = '/home/codaradm/data/grid_files/nc_templates/PR_2km_template.nc';
                title = 'CARACOOS 13MHz Sea Surface Currents';
                oldHead = 'OI_CARA';
                newHead = 'UPRM_13MHz';
                summary = ['Optimally Interpolated Total Vectors calculated by HFRProgs toolbox using MATLAB. Mercator lat/lon projection using 2km grid.'];
                keyWd = '2km';
        end;
    
    case '4'
        load('/home/codaradm/data/grid_files/nc_templates/netcdf_2d_grids/OI_1km_Grid.mat');
        output_base_dir = '/home/codaradm/data/totals/maracoos/oi/nc/25MHz';
        ncDir = '/home/om/dods-data/thredds/cool/codar/totals/25MHz_1km_realtime';
        template_filename = '/home/codaradm/data/grid_files/nc_templates/OI_1km_nc4_template.nc';
        title = 'MARACOOS 25MHz Sea Surface Currents';
        oldHead = 'OI_MARASR';
        newHead = 'RU_25MHz';
        summary = ['Optimally Interpolated Total Vectors calculated by HFRProgs toolbox using MATLAB. Mercator lat/lon projection using 1km grid.'];
        keyWd = '1km';
end

asciiFile = dir(asciiName);

lat_size = nc_getdiminfo( template_filename, 'lat' );
lon_size = nc_getdiminfo( template_filename, 'lon' );
lat_size = lat_size.Length;
lon_size = lon_size.Length;

if isempty(asciiFile);
    return;
end

% Don't process if the file is empty.
if asciiFile.bytes == 0;
    return;
end

[filePath, name, ~] = fileparts(asciiName);
source_files = [filePath '/' oldHead '_*'];

ind = find(name == '_');
yy = str2num(name(ind(2)+1:ind(3)-1));
mm = str2num(name(ind(3)+1:ind(4)-1));
dd = str2num(name(ind(4)+1:ind(5)-1));
hh = str2num(name(ind(5)+1:ind(5)+4));

gmt_time = datenum(yy, mm, dd, hh/100, 0, 0);

new_name = strrep(name, oldHead, newHead);
output_netcdf_filename = sprintf ('%s/%s.totals.nc', output_base_dir, new_name);

[success,message,messageID] = copyfile(template_filename, output_netcdf_filename);

text_file = sprintf('%s/%s', filePath, name);

fid = fopen(text_file, 'r');
header = '';
s = fgetl(fid);

while (strncmp(s, '%', 1));
    header = sprintf('%s%s\n', header, s);
    s = fgetl(fid);
end
fclose(fid);

nc_attput ( output_netcdf_filename, nc_global, 'header', header);
nc_attput ( output_netcdf_filename, nc_global, 'Conventions', 'CF 1.6');
nc_attput ( output_netcdf_filename, nc_global, 'creation_date', [datestr(now, 'ddd mmm dd HH:MM:SS yyyy') ' EDT']);
nc_attput ( output_netcdf_filename, nc_global, 'creator_name', 'John Kerfoot');
nc_attput ( output_netcdf_filename, nc_global, 'creator_email', 'kerfoot@marine.rutgers.edu');
nc_attput ( output_netcdf_filename, nc_global, 'institution', 'Coastal Ocean Observation Lab, Institute of Marine & Coastal Sciences, Rutgers University');
nc_attput ( output_netcdf_filename, nc_global, 'id', new_name);
nc_attput ( output_netcdf_filename, nc_global, 'naming_authority', 'edu.rutgers.marine.rucool');
nc_attput ( output_netcdf_filename, nc_global, 'title', title);
nc_attput ( output_netcdf_filename, nc_global, 'summary', summary);
nc_attput ( output_netcdf_filename, nc_global, 'keywords', ['codar, totals, vectors, currents, optimal, interpolation, ' keyWd]);
nc_attput ( output_netcdf_filename, nc_global, 'source', 'CODAR SeaSonde Current Mapping Device');
nc_attput ( output_netcdf_filename, nc_global, 'geospatial_lat_min', min(lat));
nc_attput ( output_netcdf_filename, nc_global, 'geospatial_lat_max', max(lat));
nc_attput ( output_netcdf_filename, nc_global, 'geospatial_lat_units', 'degrees_north');
nc_attput ( output_netcdf_filename, nc_global, 'geospatial_lon_min', min(lon));
nc_attput ( output_netcdf_filename, nc_global, 'geospatial_lon_max', max(lon));
nc_attput ( output_netcdf_filename, nc_global, 'geospatial_lon_units', 'degrees_east');
nc_attput ( output_netcdf_filename, nc_global, 'history', 'Hourly codar radial data combined into one hourly file containing vectors. See source attribute' );
nc_attput ( output_netcdf_filename, nc_global, 'source', 'surface observation' );

fprintf(1, 'Processing %s into %s\n', name, output_netcdf_filename);

process_codar_totals_into_netcdf(text_file, output_netcdf_filename, gmt_time, lat_size, lon_size, P);

%if num2str(region) == '7';
%    dodsDir = [ncDir '/' new_name '.totals.nc'];
%    disp(dodsDir)
%    copyfile(output_netcdf_filename, dodsDir)
%end
%copyfile(output_netcdf_file_name, '/home/codaradm/data/totals/maracoos/oi/ascii/5MHz/');

