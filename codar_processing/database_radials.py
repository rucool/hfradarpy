import logging
from configs.database_tables import *
from sqlalchemy import exists

logger = logging.getLogger(__name__)


def get_pattern_type_id(session, pattern_type):
    """
    From database, get the id of the pattern type specified in the hfrPatternTypes table.
    :param session: SQLAlchemy database session instance
    :param file_type: 'Ideal' or 'Measured'
    :return: ID of pattern type in hfrPatternTypes table
    """
    type_dict = dict(type=pattern_type)
    result = session.query(PatternTypes).filter_by(**type_dict).first()

    return int(result.id)


def update_latest_radials(session, filename, timestamp, site_id, radial_id):
    radial_info = dict(siteId=site_id, radialId=radial_id, TimeStamp=timestamp, filename=filename)

    result = session.query(RadialsLatest).filter_by(siteId=radial_info['siteId']).first()

    if not result:
        result = RadialsLatest(**radial_info)
        session.add(result)
        session.commit()
        session.flush()
        logging.info('{} - Table `hfrLatestRadials` - Added siteId: {} and latest timestamp: {}'.format(filename, site_id, timestamp))
    elif result:
        if result.TimeStamp < timestamp:
            result.TimeStamp = timestamp
            result.radialId = radial_id
            result.filename = filename
            session.commit()
            logging.info('{} - Table `hfrLatestRadials` - Updated siteId: {} to latest timestamp: {}'.format(filename, site_id, timestamp))


def upload_diagnostics(session, table_object, data, id):
    data['datetime'] = data['datetime'].dt.strftime('%Y-%m-%d %H:%M:%S')
    data = data.where(data.notnull(), None)
    bulk_list = []
    for row in data.itertuples():
        (ret,), = session.query(exists().where(table_object.datetime == row.datetime).
            where(table_object.id_site == id))
        if ret:
            continue
        line_dict = row._asdict()
        line_dict.pop('Index', None)
        line_ref = table_object(**line_dict)
        bulk_list.append(line_ref)

    session.bulk_save_objects(bulk_list)
    session.commit()
    session.flush()

    # data.to_sql(table_name, con=engine, if_exists='append', index=False)
    return


def upload_radial_header(session, header):
    """

    :param session:
    :param header:
    :return:
    """
    result = RadialMetadata(**header)
    session.add(result)
    session.commit()
    session.flush()
    logging.info('{} - Table `hfrRadialFilesMetadata` - Radial header upload complete.'.format(header['filename']))
    return result.id