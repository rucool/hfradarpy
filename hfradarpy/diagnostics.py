import os
import re
import numpy as np
import datetime as dt
from hfradarpy.ctf import CTFParser
from pathlib import Path
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def concat(flist):
    """
    This function takes a list of Wave objects or wave file paths and
    combines them along the time dimension using xarrays built-in concatenation
    routines.

    Args:
        wave_list (list): list of wave files or Wave objects that you want to concatenate
        enhance (bool, optional): Changes variable names to something meaningful and adds attributes. Defaults to False.

    Returns:
        xarray.Dataset: diagnostics concatenated into an xarray dataset by (time) or (time, dist)
    """
    diag_list = []
    for f in sorted(flist):
        if not isinstance(f, Diags):
            d = Diags(f)
        diag_list.append(d.data)

    df = pd.concat(diag_list)
    return df


class Diags(CTFParser):
    """
    Diagnostics Subclass.

    This class should be used when loading a CODAR diagnostic (.hdt/.rdt/.xdt) file. 
    This class inherits the generic LLUV class from hfradarpy/ctf.py 
    """

    def __init__(self, fname, ):
        """
        Initalize a diagnostic object from a diagnostic  file.

        Args:
            fname (str or path.Path): Filename to be loaded
        """
        logging.info("Loading diagnostic file: {}".format(fname))
        super().__init__(fname)

        if self._iscorrupt:
            return

        self.data = self._tables[1]["data"]
        # self.df_index = "time"  

        # Use separate date and time columns to create atetime column and drop those columns.
        self.data["time"] = self.data[["TYRS", "TMON", "TDAY", "THRS", "TMIN", "TSEC"]].apply(
            lambda s: dt.datetime(*s), axis=1
        )
        self.data = self.data.set_index('time')

        # if not self.data.empty:
        #     if replace_invalid:
        #         self.replace_invalid_values()

    def __repr__(self):
        """
        String representation of Diagnostic object

        Returns:
            str: string representation of Diagnostic object
        """
        return "<Diagnostic: {}>".format(self.file_name)


    def to_file(self, filename, validate=True, overwrite=False):
        """
        Create a diagnostic file (e.g. *.rdt, *.hdt)

        Args:
            filename (str or Path): User defined filename of diagnostic you want to save
            validate (boolean): If False, no validation check will be performed before creating the file.
            overwrite (bool): If True, an exported file can overwrite an existing file with the same name. Defaults to False.

        """
        # Make sure filename is converted into a Path object
        filename = Path(filename)

        if validate:
                if not self.is_valid():
                    raise ValueError("Could not export ASCII data, the input file was invalid.")

        # Ensure that the filename passed into the export function is not the same as the filename that we read in.
        # # We do not want to overwrite the original wave file by accident.
        if not overwrite:
            if self.full_file == str(filename):
                suffix = f'{filename.suffix}.mod'
                filename = filename.with_suffix(suffix)

        if os.path.isfile(filename):
            os.remove(filename)

        os.makedirs(os.path.dirname(filename), exist_ok=True)

        with open(filename, "w") as f:
            # Write header
            for metadata_key, metadata_value in self.metadata.items():
                if "ProcessingTool" in metadata_key:
                    break
                elif "End" in metadata_key:
                    break
                else:
                    #print(metadata_key)
                    f.write("%{}: {}\n".format(metadata_key, metadata_value))

            # Write data tables. Anything beyond the first table is commented out.
            for table in self._tables.keys():

                if "datetime" in self._tables[table]["data"].keys():
                    self._tables[table]["data"] = self._tables[table]["data"].drop(["datetime"], axis=1)

                for table_key, table_value in self._tables[table].items():
                    if table_key != 'data':
                        #if (table_key == 'TableType') & (table == 1):
                        #    f.write("%{}: {}\n".format(table_key, table_value))
                        #elif table_key == "TableColumns":
                            #f.write("%TableColumns: {}\n".format(len(self._tables[table]["data"].columns)))
                        if table_key == "TableRows":
                            f.write("%TableRows: {}\n".format(self.data.shape[0]))
                        elif table_key == "TableColumnTypes":
                            f.write("%TableColumnTypes: {}\n".format(" ".join(self.data.columns.to_list())))
                        elif table_key == "TableStart":
                            f.write("%{}: {}\n".format(table_key, table_value))
                        elif table_key == "_TableHeader":
                            pass
                        else:
                            f.write("%{}: {}\n".format(table_key, table_value))

                if table == 1:
                    # Fill NaN with 999.000 which is the standard fill value for codar lluv files
                    self.data = self.data.fillna(999.000)

                    # Convert _TableHeader to a new dataframe and concatenate to dataframe containing radial data
                    # This allows for the output format to follow CODARS CTF specifications
                    # The below block of code adds the weird header and units format that codar uses in their files

                    row_df = pd.DataFrame([self._tables[1]["_TableHeader"][0]], columns=self._tables[1]["_TableHeader"][0])
                    row_df2 = pd.DataFrame([self._tables[1]["_TableHeader"][1]], columns=self._tables[1]["_TableHeader"][0])
                    self.data.columns = self._tables[1]["_TableHeader"][0]
                    self.data = pd.concat([row_df, row_df2, self.data], ignore_index=True)
                    self.data.insert(0, "%%", np.nan)  # Insert column at the beginning of dataframe of NaNs
                    self.data.iloc[0, self.data.columns.get_loc("%%")] = "%%"  # make the first row in the first column a '%%'
                    self.data.iloc[1, self.data.columns.get_loc("%%")] = "%%"  # make the second row in the first column a '%%'

                    # Output data table to string
                    # self.data.to_string(f, index=False, justify='center', header=True, na_rep=' ')
                    self.data.temp = re.sub(
                        " %%", "%%", self.data.to_string(index=False, justify="right", header=False, na_rep=" ")
                    )

                    f.write(self.data.temp)

                f.write("\n%TableEnd: \n")
                f.write("%%\n")

            # Write footer containing processing information
            #f.write("%ProcessedTimeStamp: {}\n".format(self.metadata["ProcessedTimeStamp"]))
            for tool in self.metadata["ProcessingTool"]:
                f.write("%ProcessingTool: {}\n".format(tool))
                #f.write('%{}: {}\n'.format(footer_key, footer_value))
            f.write("%End:")

if __name__ == "__main__":
    from pathlib import Path
    
    data_path = (Path(__file__).parent.with_name("examples") / "data").resolve()
    
    hdt_file = data_path / "diagnostics" /  "SEAB" / "STAT_SEAB_2018_01_01.hdt"
    rdt_file = data_path / "diagnostics" /  "SEAB" / "STAT_SEAB_2018_01_01.rdt"    
    xdt_file = data_path / "diagnostics" /  "SEAB" / "STAT_SEAB_2018_01_01.xdt"
        
    d = Diags(hdt_file)
    print(d.file_type())
    
    d = Diags(rdt_file)
    print(d.file_type())
    
    d = Diags(xdt_file)
    print(d.file_type())
    
    from glob import glob
    
    files = glob(str(data_path / "diagnostics" /  "SEAB" / "*.hdt"))

    diag_data = concat(files)
    print(f"Concatenated Data Shape: {diag_data.shape}")
