import logging
import os
import pandas as pd
import re
from hfradarpy.configs import database_tables
from sqlalchemy.engine.url import URL
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

try:
    from hfradarpy.configs.configs import mysql_configs
except ModuleNotFoundError:
    from hfradarpy.configs.configs_default import mysql_configs


logger = logging.getLogger(__name__)


def check_freq(session, freq):
    """
    Check the frequency of the radar with the table, hfrSystemTypes. Return the frequency id
    :param session: SQLAlchemy database session instance
    :param freq: Frequency floating point number
    :return: ID of the appropriate system type in hfrSystemTypes table
    """
    if freq:
        result = session.query(database_tables.SystemTypes).filter(database_tables.SystemTypes.frequency_min < freq, freq < database_tables.SystemTypes.frequency_max).first()
    else:
        logger.warning('Center Frequency not present. Filling with null value. Please check file header')
        return None

    try:
        id = result.id
    except AttributeError:
        logger.warning('Invalid system type. Please update hfrSystemTypes with appropriate system type')
        id = None
    return id


def check_file_upload(session, f_name, tables):
    """
    This function checks if a file, f_name, has been uploaded to the database in the table defined in 'tables'
    :param session: SQLAlchemy database session instance
    :param f_name: File name to be checked
    :param tables: SQLAlchemy table object. Defined in ~/configs/database_tables.py
    :return: True - file has already been uploaded. False - if it has not been uploaded
    """

    data_dict = dict(filename=os.path.splitext(os.path.basename(f_name))[0])
    uploaded = session.query(tables).filter_by(**data_dict).first()
    if uploaded:
        logging.debug('{} - Database record found - File ID: {}'.format(f_name, uploaded))
    else:
        logging.debug('{} - No database record found. Uploading to database'.format(f_name))
    return uploaded


def db_session():
    """
    Allows the script to interact with the database via the engine
    :return: sqlalchemy session instance: http://docs.sqlalchemy.org/en/latest/orm/session_basics.html#what-does-the-session-do
    """
    return sessionmaker(bind=db_engine())()


def db_engine():
    """
    Performs database connection using database settings from settings.py.
    Returns sqlalchemy engine instance
    """
    return create_engine(URL(**mysql_configs()))


def get_file_type_id(session, file_type):
    """
    From database, get the id of the file type specified in the hfrFileTypes table. If one doesn't exist, create it.
    :param session: SQLAlchemy database session instance
    :param file_type: Four letter string that prepends CODAR CTF file name
    :return: ID of CODAR file type in hfrFileTypes table
    """
    type_dict = dict(type=file_type)
    result = session.query(database_tables.FileTypes).filter_by(**type_dict).first()

    if not result:
        logger.info('{} - New file type detected. Adding to hfrFileTypes'.format(file_type))
        new_site = dict(type=file_type)
        result = database_tables.FileTypes(**new_site)
        session.add(result)
        session.commit()
        session.flush()
    return result.id


def get_pattern_types(session):
    """
    Get a list of all pattern types in hfrPatternTypes table
    :param session: SQLAlchemy database session instance
    :return: dataframe of hfrPatternTypes
    """
    df = pd.read_sql(session.query(database_tables.PatternTypes).statement, session.bind)
    return df


def get_sites(session):
    """
    Get list of all sites in hfrSites table.

    :param session: SQLAlchemy database session instance
    :return: dataframe of hfrSites
    """
    df = pd.read_sql(session.query(database_tables.Sites).statement, session.bind)
    return df


def update_site_table(session, site, freq=None, origin=None):
    """
    Check if site exists in table, hfrSites. If it doesn't, add it to the table.
    :param session: SQLAlchemy database session instance
    :param site: Four letter code of CODAR site
    :param freq: Center Frequency (float) of CODAR site
    :param origin: Origin (lon, lat) of CODAR Receive Antenna
    :return: ID of CODAR site in hfrSites table
    """
    type_id = check_freq(session, freq)
    try:
        lonlat = re.findall(r"[-+]?\d*\.\d+|\d+", origin)
    except TypeError:
        lonlat = (0, 0)

    logger.info('{} - New HFR site detected. Adding to table `hfrSites`'.format(site))
    new_site = dict(site=site, transmitCenterFrequency=freq, lat=lonlat[0], lon=lonlat[1], type=type_id)
    ref = database_tables.Sites(**new_site)
    session.add(ref)
    session.commit()
    session.flush()
    return ref


def site_check(session, site, freq=None, origin=None):
    """
    Check if site exists in table, hfrSites. If it doesn't, add it to the table.
    :param session: SQLAlchemy database session instance
    :param site: Four letter code of CODAR site
    :param freq: Center Frequency (float) of CODAR site
    :param origin: Origin (lon, lat) of CODAR Receive Antenna
    :return: ID of CODAR site in hfrSites table
    """

    site_dict = dict(site=site)
    result = session.query(database_tables.Sites).filter_by(**site_dict).first()

    if not result:
        type_id = check_freq(session, freq)
        try:
            lonlat = re.findall(r"[-+]?\d*\.\d+|\d+", origin)
        except TypeError:
            lonlat = (0, 0)

        logger.info('{} - New HFR site detected. Adding to table `hfrSites`'.format(site))
        new_site = dict(site=site, transmitCenterFrequency=freq, lat=lonlat[0], lon=lonlat[1], type=type_id)
        ref = database_tables.Sites(**new_site)
        session.add(ref)
        session.commit()
        session.flush()
        return ref
    else:
        return result
