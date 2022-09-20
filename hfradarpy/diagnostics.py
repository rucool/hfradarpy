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
