import logging
from hfradarpy.configs import database_tables
from sqlalchemy import func

logger = logging.getLogger(__name__)


def iterate_through_data(session, df, site_info, file_id):
    """
    Loop through the data contained in the pandas dataframe
    :param session: SQLAlchemy database session instance
    :param df: Pandas dataframe containing codar data
    :param site_info:  ID of codar site in hfrSites table
    :param file_id: ID for the file metadata that was updated to hfrWaveFilesMetadata
    :return: Number of rows of data committed to hfrWaveData table
    """
    round_dict = dict(MWHT=2, MWPD=2, WAVB=2, WNDB=2, PMWH=2, DIST=4, WHSD=2)

    bulk_list = []

    # Convert pandas timestamp to a string because mysql-connector doesn't recognize 'Timestamp' objects
    df['datetime'] = df['datetime'].dt.strftime('%Y-%m-%d %H:%M:%S')

    for row in df.itertuples():
        line_dict = row._asdict()
        line_dict.pop('Index', None)
        for k, v in line_dict.items():
            if k in round_dict.keys():
                line_dict[k] = round(line_dict[k], round_dict[k])

        line_dict['site_id'] = int(site_info)
        line_dict['file_id'] = int(file_id)
        line_ref = database_tables.WaveData(**line_dict)
        bulk_list.append(line_ref)

    session.bulk_save_objects(bulk_list)
    session.commit()
    session.flush()
    return len(bulk_list)


def data_route(session, df, site_info, file_id, fname, initial_upload=True):
    """

    :param session: SQLAlchemy database session instance
    :param df: Pandas dataframe containing CODAR data
    :param site_id: ID of CODAR site in hfrSites table
    :param file_id: ID for the file metadata that was updated to hfrWaveFilesMetadata
    :param initial_upload: True for initial upload of Wave File to database. False for recurring update
    :return:
    """

    if not initial_upload:
        max_time = session.query(func.max(database_tables.WaveData.datetime)).filter_by(file_id=file_id).one()
        df = df[df.datetime > max_time[0]]

    inserted = iterate_through_data(session, df, site_info, file_id)

    if not inserted == 0:
        logger.info('{} - Inserted {} rows'.format(fname, inserted))
    else:
        logger.info('{} - Database up to date. No rows inserted'.format(fname))
