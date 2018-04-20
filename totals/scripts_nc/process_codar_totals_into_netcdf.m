function process_codar_totals_into_netcdf ( text_file, netcdf_filename, gmt_time, lat_size, lon_size, P )

afid = fopen ( text_file, 'r' );

r_time = (24*gmt_time - 24*datenum(2001,1,1))/24;

d = load(text_file);
u = d(:,3);
v = d(:,4);
u_err = d(:,5);
v_err = d(:,6);
num_radials = d(:,7);
site_code = d(:,8);

% The time is in days since 2001-1-1
nc_varput ( netcdf_filename, 'time', r_time, 0, 1 );

Data=-999*ones([1 lat_size lon_size]);

Data(P) = u;
nc_varput ( netcdf_filename, 'u', Data );

Data(P) = v;
nc_varput ( netcdf_filename, 'v', Data );

Data(P) = u_err;
nc_varput ( netcdf_filename, 'u_err', Data );

Data(P) = v_err;
nc_varput ( netcdf_filename, 'v_err', Data );

Data(P) = num_radials;
nc_varput ( netcdf_filename, 'num_radials', Data );

Data(P) = site_code;
nc_varput ( netcdf_filename, 'site_code', Data );
fclose(afid)
return

function define_netcdf_file ( netcdf_filename, num_pts )

[ncid, status] = mexnc ( 'create', netcdf_filename, nc_clobber_mode );
if status, error ( mexnc('STRERROR', status) ), end

[point_dimid] = mexnc ( 'DEF_DIM', ncid, 'point', num_pts );
if status, error ( mexnc('STRERROR', status) ), end


[time_dimid] = mexnc ( 'DEF_DIM', ncid, 'time', 1 );
if status, error ( mexnc('STRERROR', status) ), end

% time
[varid, status] = mexnc ( 'DEF_VAR', ncid, 'time', nc_double, 1, [time_dimid] );
if status, error ( mexnc('STRERROR', status) ), end

% Longitude
[varid, status] = mexnc ( 'DEF_VAR', ncid, 'lon', nc_float, 1, [point_dimid] );
if status, error ( mexnc('STRERROR', status) ), end

[varid, status] = mexnc ( 'DEF_VAR', ncid, 'lat', nc_float, 1, [point_dimid] );
if status, error ( mexnc('STRERROR', status) ), end

[varid, status] = mexnc ( 'DEF_VAR', ncid, 'u', nc_float, 1, [point_dimid] );
if status, error ( mexnc('STRERROR', status) ), end

[varid, status] = mexnc ( 'DEF_VAR', ncid, 'v', nc_float, 1, [point_dimid] );
if status, error ( mexnc('STRERROR', status) ), end

[varid, status] = mexnc ( 'DEF_VAR', ncid, 'site', nc_int, 1, [point_dimid] );
if status, error ( mexnc('STRERROR', status) ), end

[varid, status] = mexnc ( 'DEF_VAR', ncid, 'num_radials', nc_int, 1, [point_dimid] );
if status, error ( mexnc('STRERROR', status) ), end

[varid, status] = mexnc ( 'DEF_VAR', ncid, 'maperr_speed', nc_float, 1, [point_dimid] );
if status, error ( mexnc('STRERROR', status) ), end

[varid, status] = mexnc ( 'DEF_VAR', ncid, 'fitdiff', nc_float, 1, [point_dimid] );
if status, error ( mexnc('STRERROR', status) ), end

status = mexnc ( 'CLOSE', ncid );
if status, error ( mexnc('STRERROR', status) ), end

% Easier to define attributes here.
nc_attput ( netcdf_filename, 'time', 'longname', 'Time' );
nc_attput ( netcdf_filename, 'time', 'units', 'days since 2001-01-01 00:00:00' );

nc_attput ( netcdf_filename, 'lon', 'longname', 'Longitude' );
nc_attput ( netcdf_filename, 'lon', 'shortname', 'lon' );
nc_attput ( netcdf_filename, 'lon', 'units', 'degrees_east' );
nc_attput ( netcdf_filename, 'lon', '_FillValue', single(-999) );

nc_attput ( netcdf_filename, 'lat', 'longname', 'Latitude' );
nc_attput ( netcdf_filename, 'lat', 'shortname', 'lat' );
nc_attput ( netcdf_filename, 'lat', 'units', 'degrees_north' );
nc_attput ( netcdf_filename, 'lat', '_FillValue', single(-999) );

nc_attput ( netcdf_filename, 'u', 'longname', 'Eastward Velocity' );
nc_attput ( netcdf_filename, 'u', 'shortname', 'u' );
nc_attput ( netcdf_filename, 'u', 'units', 'cm/s' );
nc_attput ( netcdf_filename, 'u', '_FillValue', single(-999) );

nc_attput ( netcdf_filename, 'v', 'longname', 'Northward Velocity' );
nc_attput ( netcdf_filename, 'v', 'shortname', 'v' );
nc_attput ( netcdf_filename, 'v', 'units', 'cm/s' );
nc_attput ( netcdf_filename, 'v', '_FillValue', single(-999) );

nc_attput ( netcdf_filename, 'site', '_FillValue', int32(-999) );
nc_attput ( netcdf_filename, 'num_radials', '_FillValue', int32(-999) );

nc_attput ( netcdf_filename, 'maperr_speed', '_FillValue', single(-999) );
nc_attput ( netcdf_filename, 'fitdiff', '_FillValue', single(-999) );

return
