import logging
from configs import database_tables
from sqlalchemy import func

logger = logging.getLogger(__name__)


def iterate_through_data(session, df, site_info, file_id):
    """
    Loop through the data contained in the text file
    :param data: .readlines() data from an open text file
    :param site_id:  id for the hfr site
    :param file_id: id for the file metadata that was updated
    :param bulk: boolean whether to upload the entire file in bulk (completed wave files) or each line separately (for files still being written).
    :return: nothing
    """
    round_dict = dict(MWHT=2, MWPD=2, WAVB=2, WNDB=2, PMWH=2, DIST=4, WHSD=2)

    bulk_list = []

    # make sure that we only grab valid keys from the datasets that exist in the mysql database
    # valid_keys = [x for x in database_tables.WaveData.__table__.columns.keys() if
    #               x not in ('id', 'site_id', 'file_id', 'mwht_flag')]

    # Need to convert datetime from a pandas timestamp to a string because mysql-connector doesn't recognize 'Timestamp' objects
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

    :param df:
    :param site_id:
    :param file_id:
    :param initial_upload:
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
