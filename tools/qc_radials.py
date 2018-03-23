import functions.common as cf
import functions.qc_functions as qc
import os
from functions.radials import create_ruv


def main(radial, database=False):
    """

    :param radial:
    :param database:
    :return:
    """

    # Get directory and filename of radial file
    pathname, filename = os.path.split(radial)
    pathname_qc = os.path.join(pathname.split('radials')[0], 'radials_qc', filename[5:9])
    filename_qc = os.path.join(pathname_qc, filename)

    # Create new dir from pathname if it doesn't already exist
    cf.create_dir(pathname_qc)

    # Parse the radial file with the generic codar CTF/LLUV parser located in functions/common
    radial_file_data = cf.parse_lluv(radial)

    # Create radial_data variable with alias for main datatable in radial_file_data
    radial_data = radial_file_data['tables']['1']['data']

    # QARTOD QC Tests
    radial_data = qc.qc_location(radial_data)
    radial_data = qc.qc_speed(radial_data, 250)
    radial_file_data['header']['qc_qartod_radial_count'] = str(qc.qc_radial_count(radial_data, 150, 300))

    # Modify any tableheader information that needs to be updated
    radial_file_data['tables']['1']['TableColumns'] = radial_data.shape[1]
    header_list = radial_data.columns.tolist()
    header_list.remove('%%')
    radial_file_data['tables']['1']['TableColumnTypes'] = ' '.join(header_list)

    # Generate quality controlled radial file
    create_ruv(radial_file_data, filename_qc)

    # # TODO: Finish optional database insertion
    # if database:  # if database == True
    #
    #     # Create a global session variable
    #     global Session
    #     Session = cf.db_session()
    #
    #     try:
    #         mod_time = os.path.getmtime(radial)
    #     except OSError:
    #         mod_time = 0
    #
    #     # Create a dictionary
    #     header_dict = cf.clean_radial_header(radial_file_data['header'])
    #     header_dict['filename'] = filename
    #     header_dict['fileModTime'] = datetime.fromtimestamp(mod_time)
    #
    #     # Add database ids for the following variables (for future database relations)
    #     header_dict['processingMethod'] = int(cf.get_file_type_id(Session, filename[:4]))
    #     header_dict['PatternType'] = int(cf.get_pattern_type_id(Session, header_dict['PatternType']))
    #     header_dict['Site'] = cf.check_site(Session, header_dict['Site'],
    #                                         freq=header_dict['TransmitCenterFreqMHz'],
    #                                         origin=header_dict['Origin'])
    #
    #     # Check if this radial has been uploaded already
    #     cf.check_radial_upload(Session, header_dict['UUID'])
    #
    #     # upload radial file header information to the database
    #     id_radial = cf.upload_radial_header(Session, header_dict)
    #
    #     print


if __name__ == '__main__':
    radial = cf.path_within_module('data/radials/SEAB/2018_03/RDLi_SEAB_2018_03_01_0000.ruv')
    main(radial, database=False)