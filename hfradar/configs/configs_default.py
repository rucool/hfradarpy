import datetime as dt
import pandas as pd
from collections import OrderedDict


def mysql_configs():
    database = {
        'drivername': 'mysql+mysqlconnector',
        'username': 'admin',
        'password': 'root',
        'host': 'localhost',
        'database': 'coolops'}
    return database


def mongodb_configs():
    database = {
        'uri': 'mongodb://127.0.0.1:27017'}
    return database


def netcdf_global_attributes(required_attributes, time_start, time_end):
    """

    :param required_attributes:
    :param time_start: time string of earliest data point in file
    :param time_end: time string of latest data point in file
    :return:
    """
    datetime_format = '%Y%m%dT%H%M%SZ'
    created = pd.Timestamp(dt.datetime.utcnow()).strftime(datetime_format)  # creation time Timestamp

    time_start = pd.Timestamp(str(time_start)).strftime(datetime_format)
    time_end = pd.Timestamp(str(time_end)).strftime(datetime_format)

    # Required global attributes
    global_attrs = [('publisher_name', required_attributes['publisher_name']),
                    ('publisher_email', required_attributes['publisher_email']),
                    ('publisher_url', required_attributes['publisher_url']),
                    ('title', required_attributes['title']),
                    ('summary', required_attributes['summary']),
                    ('keywords', required_attributes['keywords']),
                    ('Conventions', 'CF-1.6, ACDD-1.3'),
                    ('naming_authority', required_attributes['naming_authority']),
                    ('history', required_attributes['history']),
                    ('source', required_attributes['source']),
                    ('processing_level', required_attributes['processing_level']),
                    ('comment', required_attributes['comment']),
                    ('acknowledgment', required_attributes['acknowledgment']),
                    ('standard_name_vocabulary', 'CF Standard Name Table v41'),
                    ('date_created', created),
                    ('creator_name', required_attributes['creator_name']),
                    ('creator_email', required_attributes['creator_email']),
                    ('creator_url', required_attributes['creator_url']),
                    ('institution', required_attributes['institution']),
                    ('project', required_attributes['project']),
                    ('geospatial_lat_min', -90),
                    ('geospatial_lat_max', 90),
                    ('geospatial_lon_min', -180),
                    ('geospatial_lon_max', 180),
                    ('geospatial_vertical_min', 0.0),
                    ('geospatial_vertical_max', 0.0),
                    ('geospatial_vertical_positive', 'down'),
                    ('time_coverage_start', time_start),
                    ('time_coverage_end', time_end),
                    ('sea_name', required_attributes['sea_name']),
                    ('creator_type', 'person'),
                    ('creator_institution', required_attributes['creator_institution']),
                    ('contributor_name', required_attributes['contributor_name']),
                    ('contributor_role', required_attributes['contributor_role']),
                    ('geospatial_lat_units', 'degrees_north'),
                    ('geospatial_lon_units', 'degrees_east'),
                    ('date_modified', created),
                    ('date_issued', created),
                    ('date_metadata_modified', created),
                    ('keywords_vocabulary', 'GCMD Science Keywords'),
                    ('platform', required_attributes['platform']),
                    ('instrument', required_attributes['instrument']),
                    ('cdm_data_type', required_attributes['cdm_data_type']),
                    ('references', required_attributes['references'])]

    if 'ncei_template_version' in required_attributes.keys():
        global_attrs = [('ncei_template_version', required_attributes['ncei_template_version'])] + global_attrs

    global_attrs = OrderedDict(global_attrs)
    return global_attrs