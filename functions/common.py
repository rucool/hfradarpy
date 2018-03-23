import io
import logging
import os
import pandas as pd
from collections import OrderedDict

logger = logging.getLogger(__name__)


def create_dir(new_dir):
    # Check if dir exists.. if it doesn't... create it.
    if not os.path.isdir(new_dir):
        try:
            os.makedirs(new_dir)
        except OSError:
            if os.path.exists(new_dir):
                pass
            else:
                raise


def dateparse(yr, mo, da, hr, mn, ss):
    dt = '{} {} {} {} {} {}'.format(yr, mo, da, hr, mn, ss)
    return pd.datetime.strptime(dt, '%Y %m %d %H %M %S')


def parse_header_line(line):
    """
    Parse a line into a key, value
    :param line: a line from a text file
    :type line: string
    :return: a tuple containing the key, value for the line
    :rtype: tuple
    """

    line = line.replace('%', '') # Strip the % sign from the line
    line = line.replace('\n', '') # Strip the new line character from the end of the line
    key = line[0:line.find(':')] # Grab the key, which is everything to the left of the colon
    value = line[line.find(':') + 2:] # Grab the value, which is everything to the right of the colon
    return key, value


def parse_metadata(codar_file):
    metadata_list = []

    for i, line in enumerate(codar_file):
        if line.startswith('%'):
            metadata_list.append([i, line.strip('\n')])
    return metadata_list


def parse_lluv(lluv):
    # Load file into Generic LLUV container
    lluv_container = dict(header=OrderedDict(), tables=OrderedDict(), footer=OrderedDict())

    with open(lluv, 'r') as open_file:
        open_lluv = open_file.readlines()

    # Parse header and footer metadata
    table_count = 0
    table = False  # Set table to False. Once a table is found, switch to True.
    processing_info = []
    for i, line in enumerate(open_lluv):
        if not table:  # If we are not looking at a table
            if line.startswith('%%'):
                continue
            elif line.startswith('%'):  # Parse the single commented header lines
                key, value = parse_header_line(line)
                if 'TableType' in line:  # Save this data as global header information
                    table = True  # we found a table
                    table_count = table_count + 1  # this is the nth table
                    data_header = []  # initialize an empty list for the data header information
                    table_data = u''
                    lluv_container['tables'][str(table_count)] = OrderedDict()
                    lluv_container['tables'][str(table_count)][key] = value
                elif table_count > 0:
                    if key == 'ProcessingTool':
                        processing_info.append(value)
                    elif key == 'End':
                        lluv_container['footer']['ProcessingTool'] = processing_info
                        lluv_container['footer'][key] = value
                    else:
                        lluv_container['footer'][key] = value
                else:
                    lluv_container['header'][key] = value
        elif table:
            if line.startswith('%%'):  # table header information
                line = line.replace('%%', '')
                temp = line.split()
                if 'comp' in temp:
                    temp = [x for x in temp if x not in ('comp', 'Distance')]
                data_header.append(tuple(temp))
            elif line.startswith('%'):
                if len(line.split(':')) == 1:
                    line = line.replace('%', '')
                    table_data += '{}'.format(line)
                    # table_data.append(line.split())
                else:
                    key, value = parse_header_line(line)
                    if 'TableEnd' in line:
                        # use pandas read_csv rather than loading directly into dataframe because it automatically
                        # interprets the datatype for each column of the csv
                        tdf = pd.read_csv(io.StringIO(unicode(table_data)), sep=' ', header=None,
                                          names=lluv_container['tables'][str(table_count)]['TableColumnTypes'].split(),
                                          skipinitialspace=True)

                        # tdf = pd.DataFrame(table_data, columns=lluv_container['tables'][str(table_count)]['TableColumnTypes'].split())
                        if table_count > 1:
                            tdf.insert(0, '%%', '%')
                        else:
                            tdf.insert(0, '%%', '')
                        lluv_container['tables'][str(table_count)]['data'] = tdf
                        # lluv_container['tables'][str(table_count)]['TableHeader'] = data_header
                        table = False
                    else:
                        lluv_container['tables'][str(table_count)][key] = value
            else:  # Uncommented lines are the main data table.
                table_data += '{}'.format(line)
                # table_data.append(line.split())
    return lluv_container


def parse_header(codar_file):
    """
    Parse the generalized seasonde lluv format header data into a dictionary
    :param data: .readlines() data from an open text file
    :return header_dict: dictionary containing all the important header information for each file
    :rtype: dictionary
    """
    header_dict = {} # initialize an empty dictionary
    dist_dict = {}
    # headers = None
    # col_headers = None
    units = None
    n = 0
    for i, line in enumerate(codar_file):
        if '%' in line:
            if 'Distance:' in line:
                key, value = parse_header_line(line)
                key = '{}_{}'.format(key.lower(), str(n))
                dist_dict[key] = float(value.split(' ')[0])
                n = n + 1
            elif 'TableStart' in line:
                headers = i + 1
                units = headers + 1
            # elif i == headers:
            #     col_headers = [x.strip() for x in line.split('  ')[1:] if x]
            #     # headers += ['', '']  # Site Contributors
            elif i == units:
                # units = [x.strip().replace('(', '').replace(')', '') for x in line.split('  ')[1:] if x]
                break
            else:
                key, value = parse_header_line(line)
                header_dict[key] = value
    return header_dict, dist_dict #, col_headers


def path_within_module(file_path):
    complete_path = os.path.join(os.path.dirname(__file__), '..', file_path)
    return complete_path