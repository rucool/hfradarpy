# codar_processing
Rutgers Center for Ocean Observing Leadership High Frequency Radar (CODAR) Processing toolbox

## Install Miniconda
Download and follow installation instructions for the appropriate Miniconda installer from http://conda.pydata.org/miniconda.html. 

Make sure to add the channel, conda-forge, to your .condarc. You can find out more about conda-forge from their website: https://conda-forge.org/
 
You can do this with the following command:

`conda config --add channels conda-forge`

## Clone codar_processing repository
Either use git to clone:

`git clone https://github.com/rucool/codar_processing.git`

or download the zip from file from, https://github.com/rucool/codar_processing/archive/master.zip and extract to the directory of your choice


## Create  environment
Change your current working directory to the location that you downloaded codar_processing to. 

`cd /Users/mikesmith/Documents/git/codar_processing/`

Create conda environment from the included environment.yml file:

`conda env create -f environment.yml`

Once the environment is done building, you can activate the environment by typing:

    conda activate codar_processing # OSX/Unix
    
## Install toolbox to environment

Now we need to install the toolbox to the conda environment. We can do this as follows from the root directory of the codar_processing toolbox:

`pip install .`

The toolbox should now be installed to your conda environment.

If you are developing new code in the toolbox you should installed this library as "editable":

`pip install --no-deps --force-reinstall --ignore-installed -e .`

## Using the toolbox
### Open your python interpreter
    (codar_processing)
    mikesmith@nbp-56-180: ~
    $ python
    Python 3.6.5 | packaged by conda-forge | (default, Apr  6 2018, 13:44:09)
    [GCC 4.2.1 Compatible Apple LLVM 6.1.0 (clang-602.0.53)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>>
    
### Import CODAR file classes
The toolbox is based on the LLUV file format developed by CODAR Ocean Sensors. There is a generic, inheritable class called 'LLUV' that is written to dynamically read in all CODAR LLUV format files.
There are currently two subclasses, Radial and Waves, that inherit from this LLUV class. Other LLUV files are currently being implemented. In order to load one of these files, you will need to import
one (or all) of the the classes into your python environment before you can utilize the toolbox.
  
* For loading radial files: `from codar_processing.radials import Radial`
* For loading wave files: `from codar_processing.waves import Waves`

### Loading the CODAR files using the imported classes
    In [1]: radial_file = '/Users/mikesmith/Documents/git/rucool/codar_processing/data/radials/SEAB/2018_03/RDLi_SEAB_2018_03_01_0000.ruv'
    In [2]: r = Radial(radial_file)
    
#### This command lists methods (functions) of the built-in class. Don't bother with any methods starting with a _   
    In [3]: r.__dir__() 
    ['file_path', 'file_name', 'header', 'tables', 'footer', 'data', 'diags_radial', 'diags_hardware', '__module__', '__doc__', '__init__', 'file_type', 'validate_key', 'validate_header', 'create_ruv', 'export', 'qc_qartod_location', 'qc_qartod_radial_count', 'qc_qartod_speed', 'reset', '__metaclass__', 'is_valid', '_parse_header_line', '__dict__', '__weakref__', '__repr__', '__hash__', '__str__', '__getattribute__', '__setattr__', '__delattr__', '__lt__', '__le__', '__eq__', '__ne__', '__gt__', '__ge__', '__new__', '__reduce_ex__', '__reduce__', '__subclasshook__', '__init_subclass__', '__format__', '__sizeof__', '__dir__', '__class__']

#### List the file name    
    In [4]: r.file_name
    'RDLi_SEAB_2018_03_01_0000.ruv'
    
#### List the file path    
    In [5]: r.file_path
    '/Users/mikesmith/Documents/git/rucool/codar_processing/data/radials/SEAB/2018_03/RDLi_SEAB_2018_03_01_0000.ruv'

#### List the file type
    In [6]: r.file_type()
    'radial'

#### List header information
    In [6]: import pprint
    In [7]: pprint.pprint(r.header)
    OrderedDict([('CTF', '1.00'),
        ('FileType', 'LLUV rdls "RadialMap"'),
        ('LLUVSpec', '1.26  2016 10 07'),
        ('UUID', '0D8D79F8-6BD0-4A82-8DD5-B7AD5F924C78'),
        ('Manufacturer', 'CODAR Ocean Sensors. SeaSonde'),
        ('Site', 'SEAB ""'),
        ('TimeStamp', '2018 03 01  00 00 00'),
        ('TimeZone', '"UTC" +0.000 0 "Atlantic/Reykjavik"'),
        ('TimeCoverage', '75.000 Minutes'),
        ('Origin', '40.3668167  -73.9735333'),
        ('GreatCircle', '"WGS84" 6378137.000  298.257223562997'),
        ('GeodVersion', '"CGEO" 1.70  2014 09 09'),
        ('LLUVTrustData', 'all  all lluv xyuv rbvd'),
        ('RangeStart', '2'),
        ('RangeEnd', '21'),
        ('RangeResolutionKMeters', '3.020300'),
        ('RangeCells', '31'),
        ('DopplerCells', '512'),
        ('DopplerInterpolation', '2'),
        ('AntennaBearing', '152.0 True'),
        ('ReferenceBearing', '0 True'),
        ('AngularResolution', '5 Deg'),
        ('SpatialResolution', '5 Deg'),
        ('PatternType', 'Ideal'),
        ('PatternDate', '2110 05 10  14 59 41'),
        ('PatternResolution', '1.0 deg'),
        ('TransmitCenterFreqMHz', '13.450000'),
        ('TransmitBandwidthKHz', '-49.629688'),
        ('TransmitSweepRateHz', '2.000000'),
        ('DopplerResolutionHzPerBin', '0.001953125'),
        ('FirstOrderMethod', '0'),
        ('BraggSmoothingPoints', '4'),
        ('CurrentVelocityLimit', '100.0'),
        ('BraggHasSecondOrder', '1'),
        ('RadialBraggPeakDropOff', '100.000'),
        ('RadialBraggPeakNull', '6.310'),
        ('RadialBraggNoiseThreshold', '5.000'),
        ('PatternAmplitudeCorrections', '1.0000  1.0000'),
        ('PatternPhaseCorrections', '50.00  62.00'),
        ('PatternAmplitudeCalculations', '0.2313  0.4446'),
        ('PatternPhaseCalculations', '50.90  61.40'),
        ('RadialMusicParameters', '40.000 20.000 2.000'),
        ('RadialMinimumMergePoints', '2'),
        ('FirstOrderCalc', '1'),
        ('MergeMethod', '1 MedianVectors'),
        ('PatternMethod', '1 PatternVectors'),
        ('MergedCount', '7')])
        
#### List footer information
    In [8]: pprint.pprint(r.footer)
    OrderedDict([('ProcessedTimeStamp', '2018 03 01  00 38 56'),
        ('ProcessingTool',
         ['"RadialMerger" 11.5.0',
         '"SpectraToRadial" 11.5.1',
         '"RadialSlider" 12.0.0',
         '"RadialArchiver" 12.0.0',
         '"AnalyzeSpectra" 10.9.6']),
        ('End', '')])
        
#### Look at radial data
    In [9]: r.data
                 LOND       LATD    VELU    VELV  VFLG     ESPC     ETMP    MAXV    MINV  ERSC  ERTC     XDST     YDST     RNGE   BEAR    VELO   HEAD  SPRC
    0      -73.971049  40.421183   0.413  11.818   128    1.089    2.614  -9.647 -11.825     2     5   0.2108   6.0369   6.0406    2.0 -11.825  182.0     2
    1      -73.964859  40.420810   1.974  16.060   128  999.000  999.000 -16.181 -16.181     1     2   0.7362   5.9956   6.0406    7.0 -16.181  187.0     2
    2      -73.946871  40.417252  15.865  39.234   128  999.000   13.349 -42.320 -42.320     1     3   2.2628   5.6007   6.0406   22.0 -42.320  202.0     2
    3      -73.941222  40.415282  15.763  30.909   128   15.940   11.929  -9.647 -48.855     3     5   2.7424   5.3822   6.0406   27.0 -34.696  207.0     2
    4      -73.935820  40.412944  22.441  35.880   128  999.000   15.095 -42.320 -42.320     1     5   3.2010   5.1227   6.0406   32.0 -42.320  212.0     2
    [758 rows x 19 columns] # Cut down table for brevity


#### Look at radial diagnostic data
    In [5]: r.diags_radial
    %%  TIME   AMP1   AMP2  PH13  PH23  CPH1  CPH2   SNF1   SNF2         ...           RABA  RTYP  STYP  TYRS  TMON  TDAY  THRS  TMIN  TSEC            datetime
    0  % -1800  0.239  0.498  49.1  61.6  50.0  62.0 -144.0 -144.0         ...          111.3     1    68  2018     2    28    23    30     0 2018-02-28 23:30:00
    1  % -1200  0.241  0.473  49.3  60.3  50.0  62.0 -143.0 -144.0         ...          112.6     1    68  2018     2    28    23    40     0 2018-02-28 23:40:00
    2  %  -600  0.240  0.532  48.6  59.1  50.0  62.0 -144.0 -144.0         ...          114.9     1    68  2018     2    28    23    50     0 2018-02-28 23:50:00
    3  %     0  0.231  0.493  51.0  61.3  50.0  62.0 -144.0 -145.0         ...          114.3     1    68  2018     3     1     0     0     0 2018-03-01 00:00:00
    4  %   600  0.236  0.471  50.0  62.4  50.0  62.0 -144.0 -145.0         ...          115.7     1    68  2018     3     1     0    10     0 2018-03-01 00:10:00
    5  %  1200  0.242  0.484  49.4  60.3  50.0  62.0 -144.0 -144.0         ...          121.0     1    68  2018     3     1     0    20     0 2018-03-01 00:20:00
    6  %  1800  0.231  0.445  50.9  61.4  50.0  62.0 -143.0 -145.0         ...          117.8     1    68  2018     3     1     0    30     0 2018-03-01 00:30:00
    
#### Look at hardware diagnostic data 
    r.diags_hardware
    %%  TIME  RTMP  MTMP  XTRP     RUNT  SP24  SP05  SN05   SP12         ...          EXTA  EXTB     CRUN  TYRS  TMON  TDAY  THRS  TMIN  TSEC            datetime
    % -35.0    29    41     0  5330369   0.0  5.08 -5.03  12.16         ...             0     0  1406.61  2018     2    28    23    25     0 2018-02-28 23:25:00
    % -30.0    29    41     0  5330668   0.0  5.08 -5.05  12.16         ...             0     0  1411.78  2018     2    28    23    30     0 2018-02-28 23:30:00
    % -25.0    29    41     0  5330967   0.0  5.08 -5.03  12.16         ...             0     0  1416.62  2018     2    28    23    35     0 2018-02-28 23:35:00
    % -20.0    29    41     0  5331267   0.0  5.08 -5.05  12.16         ...             0     0  1421.64  2018     2    28    23    40     0 2018-02-28 23:40:00
    % -15.0    29    41     0  5331567   0.0  5.08 -5.05  12.16         ...             0     0  1426.65  2018     2    28    23    45     0 2018-02-28 23:45:00
    % -10.0    29    41     0  5331866   0.0  5.08 -5.05  12.16         ...             0     0  1431.66  2018     2    28    23    50     0 2018-02-28 23:50:00
    %  -5.0    29    41     0  5332166   0.0  5.08 -5.03  12.16         ...             0     0  1436.67  2018     2    28    23    55     0 2018-02-28 23:55:00
    %   0.0    28    41     0  5332456   0.0  5.08 -5.05  12.16         ...             0     0     1.47  2018     3     1     0     0     0 2018-03-01 00:00:00
    %   5.0    29    41     0  5332753   0.0  5.08 -5.05  12.16         ...             0     0     6.47  2018     3     1     0     5     0 2018-03-01 00:05:00
    %  10.0    29    41     0  5333053   0.0  5.08 -5.05  12.16         ...             0     0    11.48  2018     3     1     0    10     0 2018-03-01 00:10:00
    %  15.0    29    41     0  5333352   0.0  5.08 -5.05  12.16         ...             0     0    16.49  2018     3     1     0    15     0 2018-03-01 00:15:00
    %  20.0    29    41     0  5333652   0.0  5.08 -5.03  12.16         ...             0     0    21.50  2018     3     1     0    20     0 2018-03-01 00:20:00
    %  25.0    29    41     0  5333951   0.0  5.08 -5.03  12.16         ...             0     0    26.70  2018     3     1     0    25     0 2018-03-01 00:25:00
    %  30.0    29    41     0  5334251   0.0  5.08 -5.05  12.16         ...             0     0    31.52  2018     3     1     0    30     0 2018-03-01 00:30:00
    %  35.0    29    41     0  5334549   0.0  5.08 -5.05  12.16         ...             0     0    36.53  2018     3     1     0    35     0 2018-03-01 00:35:00
    [15 rows x 35 columns]


## Running tests

After setting up your environment you can run all of the tests:

`pytest`
