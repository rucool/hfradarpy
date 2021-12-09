import datetime as dt
import io
import os
import re
from abc import ABCMeta, abstractmethod
from collections import OrderedDict
import numpy as np
import pandas as pd

import logging
logger = logging.getLogger(__name__)


class CTFParser(object):
    """
    A generic parser for the CODAR CTF file format.
    url: http://support.codar.com/Technicians_Information_Page_for_SeaSondes/Manuals_Documentation_Release_8/File_Formats/File_CTF.pdf
    """

    __metaclass__ = ABCMeta

    def __init__(self, fname):
        """
        Return an LLUVParser object whose
        """
        split_path = os.path.split(fname)
        self.file_path = split_path[0]
        self.file_name = split_path[1]
        self.full_file = os.path.realpath(fname)
        self.metadata = OrderedDict()
        self._tables = OrderedDict()

        # Load the LLUV Data with this generic LLUV parsing routine below
        table_count = 0
        table = False  # Set table to False. Once a table is found, switch to True.
        self.is_wera = False  # Default false. If 'WERA' is detected in the Manufacturer flag. It is set to True
        processing_info = []

        with open(self.full_file, 'r', encoding='ISO-8859-1') as open_file:
            open_lluv = open_file.readlines()
            if any('%End:' in s or s.strip() == '%End' for s in open_lluv):  # if there is no %End: the file is corrupt!
                # Parse header and footer metadata
                for line in open_lluv:

                    # Fix for older WERA files
                    # Add a colon to the end of '%End'
                    if line.strip() == '%End':
                        line += ':'

                    if not table:  # If we are not looking at a table or a tables header information
                        if line.startswith('%%'):
                            continue
                        elif line.startswith('%'):  # Parse the single commented header lines
                            key, value = self._parse_header_line(line)
                            if 'TableType' in line:  # Save this data as global header information
                                table = True  # we found a table
                                table_count = table_count + 1  # this is the nth table
                                table_data = u''
                                # self._data_header[table_count] = []
                                self._tables[str(table_count)] = OrderedDict()
                                self._tables[str(table_count)][key] = value
                                self._tables[str(table_count)]['_TableHeader'] = []
                            elif 'Manufacturer' in line:
                                if 'WERA' in value:
                                    self.is_wera = True
                                self.metadata[key] = value
                            elif table_count > 0:
                                if key == 'ProcessingTool':
                                    processing_info.append(value)
                                else:
                                    self.metadata[key] = value
                            else:
                                self.metadata[key] = value
                    elif table:
                        if line.startswith(('%', ' %')):
                            if line.startswith(('%%', ' %%')):  # table header information
                                rep = {' comp': '_comp', ' Distance': '_Distance',' Ratio': '_Ratio',' (dB)': '_(dB)',' Width': '_Width', ' Resp': '_Resp', 'Value ': 'Value_','FOL ':'FOL_' }
                                rep = dict((re.escape(k), v) for k, v in rep.items())
                                pattern = re.compile('|'.join(rep.keys()))
                                temp = pattern.sub(lambda m: rep[re.escape(m.group(0))], line).strip('% \n')
                                temp = [x.replace('_', ' ') for x in re.sub(' +', ' ', temp).split(' ')]  # Get rid of underscores
                                # temp[0] = '%%   {}'.format(temp[0])

                                self._tables[str(table_count)]['_TableHeader'].append(temp)
                            else:  # Table metadata and diagnostic data are prepended by at least 1 % sign
                                if len(line.split(':')) == 1:  # Diagnostic Data
                                    line = line.replace('%', '').strip()
                                    table_data += '{}\n'.format(line)
                                else:  # Table data
                                    key, value = self._parse_header_line(line)
                                    # if 'TableColumnTypes' not in self._tables[str(table_count)]:
                                    #     raise ValueError("TableColumnTypes not defined")
                                    if 'TableEnd' in line:
                                        if 'TableColumnTypes' in self._tables[str(table_count)]:
                                            # use pandas read_csv because it interprets the datatype for each column of the csv
                                            tdf = pd.read_csv(
                                                io.StringIO(table_data),
                                                sep=' ',
                                                header=None,
                                                names=self._tables[str(table_count)]['TableColumnTypes'].split(),
                                                skipinitialspace=True
                                            )
                                        else:
                                            tdf = pd.DataFrame()

                                        self._tables[str(table_count)]['data'] = tdf
                                        table = False
                                    else:
                                        key, value = self._parse_header_line(line)
                                        self._tables[str(table_count)][key] = value
                        else:  # Uncommented lines are the main data table.
                            table_data += '{}'.format(line)
                self.metadata['ProcessingTool'] = processing_info
                # if self.is_wera:
                #     self._tables['1']['data'] = pd.read_csv(io.StringIO(table_data),
                #                                             sep=' ',
                #                                             header=None,
                #                                             names=['LOND', 'LATD', 'VELU', 'VELV', 'VFLG', 'EACC', 'RNGE', 'BEAR', 'VELO', 'HEAD'],  # WERA has incorrect TableColumnTypes in their files.....
                #                                             skipinitialspace=True, )
                self._iscorrupt = False
            else:
                logging.error('{}: File corrupt. Skipping to next file.'.format(self.full_file))
                self._iscorrupt = True
        try:
            self.time = dt.datetime(*[int(s) for s in self.metadata['TimeStamp'].split()])
        except KeyError:
            pass

    def is_valid(self, table='1'):
        """
        Check if the data table for the file contains data
        :param table: string containing the table number to validate. Defaults to the primary data table '1'
        :return: True or False
        """
        try:
            return not self._tables[table]['data'].empty
        except:
            return False

    @staticmethod
    def _parse_header_line(line):
        """
        Parse a line into a key, value
        :param line: a line from a text file
        :type line: string
        :return: a tuple containing the key, value for the line
        :rtype: tuple
        """

        line = line.replace('%', '')  # Strip the % sign from the line
        line = line.replace('\n', '')  # Strip the new line character from the end of the line
        line_split = line.split(':')
        key = line_split[0]  # save key variable
        value = line_split[1].strip()  # save value variable and strip whitespace from value
        return key, value

    @abstractmethod
    def file_type(self):
        """Return a string representing the type of file this is."""
        pass

    def replace_invalid_values(self, values=[999.00, 1080.0]):
        """
        Convert invalid CODAR values to NaN
        :param df: dataframe
        :param values: list of CODAR fill values that reflect non calculable values
        :return: dataframe with invalid values set to NaN
        """
        logging.info('Replacing invalid values {} with NaN'.format(values))
        self.data.replace(values, np.nan, inplace=True)
