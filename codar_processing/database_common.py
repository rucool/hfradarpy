import logging
import os
import re
from configs import database_tables
from configs.configs import db_configs
from sqlalchemy.engine.url import URL
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

logger = logging.getLogger(__name__)


def check_freq(session, freq):
    """
    Grab the id from the table, hfrSystemTypes, where the frequency lies between
    :param session: database session connector
    :param freq: frequency floating point number
    :return: id of the appropriate system type
    """
    result = session.query(database_tables.hfrSystemTypes).filter(database_tables.hfrSystemTypes.frequency_min < freq, freq < database_tables.hfrSystemTypes.frequency_max).first()

    try:
        id = result.id
    except AttributeError:
        logger.warn('Invalid system type. Please update hfrSystemTypes with appropriate system type')
        id = None

    return id


def check_file_upload(session, f_name, tables):
    """
    :param f_name: file name of the spreadsheet to be checked
    :return: returns true or false if the spreadsheet has already been uploaded before
    """
    data_dict = dict(filename=os.path.splitext(os.path.basename(f_name))[0])
    uploaded = session.query(tables).filter_by(**data_dict).first()
    logging.debug('{} - Database record found - File ID: {}'.format(f_name, uploaded))
    return uploaded


def db_session():
    """
    Allows the script to interact with the database via the engine
    :return: sqlalchemy session instance
    """
    return sessionmaker(bind=db_engine())()


def db_engine():
    """
    Performs database connection using database settings from settings.py.
    Returns sqlalchemy engine instance
    """
    return create_engine(URL(**db_configs()))


def get_file_type_id(session, file_type):
    """

    :param file_type: Four letter string that prepends CODAR CTF file name
    :return: hfrFileTypes id
    """
    type_dict = dict(type=file_type)
    result = session.query(database_tables.hfrFileTypes).filter_by(**type_dict).first()

    if not result:
        logger.info('{} - New file type detected. Adding to hfrFileTypes'.format(file_type))
        new_site = dict(type=file_type)
        result = database_tables.hfrFileTypes(**new_site)
        session.add(result)
        session.commit()
        session.flush()
    return result.id


def get_pattern_type_id(session, pattern_type):
    """

    :param file_type: Four letter string that prepends CODAR CTF file name
    :return: hfrFileTypes id
    """
    type_dict = dict(type=pattern_type)
    result = session.query(database_tables.hfrPatternTypes).filter_by(**type_dict).first()

    return int(result.id)


def site_check(session, site, freq=None, origin=None):
    """
    Check if site exists in database. If it doesn't, create it.
    :param session:
    :param site:
    :param freq:
    :param origin:
    :return:
    """

    site_dict = dict(site=site)
    result = session.query(database_tables.Sites).filter_by(**site_dict).first()

    if freq and origin:
        lonlat = re.findall(r"[-+]?\d*\.\d+|\d+", origin)
        if not result:
            logger.info('{} - New HFR site detected. Adding to table `hfrSites`'.format(site))
            type_id = check_freq(session, freq)
            new_site = dict(site=site, transmitCenterFrequency=freq, lat=lonlat[0], lon=lonlat[1], type=type_id)
            ref = database_tables.Sites(**new_site)
            session.add(ref)
            session.commit()
            session.flush()
            ref_id = ref.id
            return ref_id
        else:
            return result
    else:
        return result


def upload_radial_header(session, header):
    result = database_tables.hfrRadialFilesMetadata(**header)
    session.add(result)
    session.commit()
    session.flush()
    logger.debug('Header uploaded')
    return result.id