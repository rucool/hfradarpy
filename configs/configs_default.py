import datetime as dt
from collections import OrderedDict


def db_configs():
    database = {
        'drivername': 'mysql',
        'username': 'admin',
        'password': 'root',
        'host': 'localhost',
        'database': 'coolops'}
    return database


def netcdf_global_attributes(user_attributes, time_string):
    """
    define global attributes
    :param lon: longitude data
    :param lat: latitude data
    :param time_string: time string of hourly file
    :return:
    """
    created = dt.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ') # Timestamp to add to global attributes for creation time

    global_attrs = [('ncei_template_version', 'NCEI_NetCDF_Grid_Template_v2.0'),
                    ('title', user_attributes['title']),
                    ('summary', 'Optimally Interpolated Total Vectors calculated by HFRProgs toolbox using MATLAB. Mercator lat/lon projection'),
                    ('keywords', 'Environmental Advisories > Marine Advisories > Marine Weather/Forecast, Oceans > Coastal Processes, Oceans > Ocean Circulation, Oceans > Ocean Waves, Oceans > Ocean Winds, Oceans > Ocean Tides, Spectral/Engineering > Radar'),
                    ('Conventions', 'CF-1.6, ACDD-1.3'),
                    ('naming_authority', 'edu.rutgers.marine.rucool'),
                    ('history', 'Hourly codar radial data combined into one hourly file containing vectors.'),
                    ('source', 'CODAR SeaSonde Surface Current Mapping Device'),
                    ('processing_level', 'Level 3'),
                    ('comment', user_attributes['comment']),
                    ('acknowledgment', user_attributes['acknowledgment']),
                    ('standard_name_vocabulary', 'CF Standard Name Table v41'),
                    ('date_created', created),
                    ('creator_name', user_attributes['creator_name']),
                    ('creator_email', user_attributes['creator_email']),
                    ('creator_url', user_attributes['creator_url']),
                    ('institution', user_attributes['institution']),
                    ('project', user_attributes['project']),
                    ('publisher_name', 'NOAA National Centers for Environmental Information'),
                    ('publisher_email', 'ncei.info@noaa.gov'),
                    ('publisher_url', 'www.ncei.noaa.gov'),
                    ('geospatial_lat_min', -90),
                    ('geospatial_lat_max', 90),
                    ('geospatial_lon_min', -180),
                    ('geospatial_lon_max', 180),
                    ('geospatial_vertical_min', 0.0),
                    ('geospatial_vertical_max', 0.0),
                    ('geospatial_vertical_positive', 'down'),
                    ('time_coverage_start', time_string),
                    ('time_coverage_end', time_string),
                    ('sea_name', user_attributes['sea_name']),
                    ('creator_type', 'person'),
                    ('creator_institution', user_attributes['creator_institution']),
                    ('contributor_name', user_attributes['contributor_name']),
                    ('contributor_role', user_attributes['contributor_role']),
                    ('geospatial_lat_units', 'degrees_north'),
                    ('geospatial_lon_units', 'degrees_east'),
                    ('date_modified', created),
                    ('date_issued', created),
                    ('date_metadata_modified', created),
                    ('keywords_vocabulary', 'GCMD Science Keywords'),
                    ('platform', user_attributes['platform']),
                    ('instrument', user_attributes['instrument']),
                    ('cdm_data_type', 'Grid'),
                    ('references', user_attributes['references'])]

    global_attrs = OrderedDict(global_attrs)
    return global_attrs