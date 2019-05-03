from sqlalchemy import Column, ForeignKey
from sqlalchemy.dialects.mysql import DOUBLE, TINYINT, TIMESTAMP, VARCHAR, FLOAT, CHAR, INTEGER, SMALLINT, DATETIME
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class FileTypes(Base):
    __tablename__ = 'hfrFileTypes'

    id = Column(INTEGER, primary_key=True)
    type = Column(CHAR)
    description = Column(VARCHAR)


class HardwareDiagnostics(Base):
    __tablename__ = 'hfrHardwareDiagnostics'

    id = Column(INTEGER, primary_key=True)
    id_site = Column(INTEGER)
    id_radial = Column(INTEGER)
    RTMP = Column(TINYINT)
    MTMP = Column(TINYINT)
    XTRP = Column(VARCHAR)
    RUNT = Column(INTEGER)
    SP24 = Column(DOUBLE)
    SP05 = Column(DOUBLE)
    SN05 = Column(DOUBLE)
    SP12 = Column(DOUBLE)
    XPHT = Column(TINYINT)
    XAHT = Column(TINYINT)
    XAFW = Column(TINYINT)
    XARW = Column(TINYINT)
    X1VR = Column(DOUBLE)
    XP28 = Column(DOUBLE)
    XP05 = Column(DOUBLE)
    GRMD = Column(TINYINT)
    GDMD = Column(TINYINT)
    GSLK = Column(TINYINT)
    GSUL = Column(TINYINT)
    PLLL = Column(TINYINT)
    HTMP = Column(DOUBLE)
    HUMI = Column(TINYINT)
    RBIA = Column(DOUBLE)
    EXTA = Column(TINYINT)
    EXTB = Column(TINYINT)
    CRUN = Column(DOUBLE)
    datetime = Column(DATETIME)


class PatternTypes(Base):
    __tablename__ = 'hfrPatternTypes'

    id = Column(INTEGER, primary_key=True)
    type = Column(CHAR)


class QCValues(Base):
    __tablename__ = 'hfrQCValues'

    id = Column(INTEGER, ForeignKey('hfrSites.id'), primary_key=True)
    radial_min_count = Column(SMALLINT)
    radial_low_count = Column(SMALLINT)
    radial_max_speed = Column(SMALLINT)
    radial_smed_range_cell_limit = Column(FLOAT)
    radial_smed_angular_limit = Column(SMALLINT)
    radial_smed_current_difference = Column(SMALLINT)


class RadialDiagnostics(Base):
    __tablename__ = 'hfrRadialDiagnostics'

    id = Column(INTEGER, primary_key=True)
    id_site = Column(INTEGER)
    id_radial = Column(INTEGER)
    AMP1 = Column(DOUBLE)
    AMP2 = Column(DOUBLE)
    PH13 = Column(DOUBLE)
    PH23 = Column(DOUBLE)
    CPH1 = Column(DOUBLE)
    CPH2 = Column(DOUBLE)
    SNF1 = Column(DOUBLE)
    SNF2 = Column(DOUBLE)
    SNF3 = Column(DOUBLE)
    SSN1 = Column(DOUBLE)
    SSN2 = Column(DOUBLE)
    SSN3 = Column(DOUBLE)
    DGRC = Column(TINYINT)
    DOPV = Column(SMALLINT)
    DDAP = Column(TINYINT)
    RADV = Column(SMALLINT)
    RAPR = Column(TINYINT)
    RARC = Column(TINYINT)
    RADR = Column(DOUBLE)
    RMCV = Column(DOUBLE)
    RACV = Column(DOUBLE)
    RABA = Column(DOUBLE)
    RTYP = Column(TINYINT)
    STYP = Column(TINYINT)
    datetime = Column(DATETIME)
    upload_timestamp = Column(TIMESTAMP)


class RadialsLatest(Base):
    __tablename__ = 'hfrLatestRadials'
    siteId = Column(INTEGER, primary_key=True)
    radialId = Column(INTEGER)
    TimeStamp = Column(DATETIME)
    filename = Column(VARCHAR)
    dateAddedEasternTime = Column(TIMESTAMP)


class RadialMetadata(Base):
    __tablename__ = 'hfrRadialFilesMetadata'

    id = Column(INTEGER, primary_key=True)
    filename = Column(VARCHAR)
    CTF = Column(DOUBLE)
    FileType = Column(VARCHAR)
    LLUVSpec = Column(VARCHAR)
    UUID = Column(VARCHAR)
    Manufacturer = Column(VARCHAR)
    Site = Column(INTEGER)
    TimeStamp = Column(DATETIME)
    TimeZone = Column(CHAR)
    TimeCoverage = Column(SMALLINT)
    Origin = Column(VARCHAR)
    GreatCircle = Column(VARCHAR)
    GeodVersion = Column(VARCHAR)
    LLUVTrustData = Column(VARCHAR)
    RangeStart = Column(TINYINT)
    RangeEnd = Column(TINYINT)
    RangeResolutionKMeters = Column(FLOAT)
    AntennaBearing = Column(SMALLINT)
    ReferenceBearing = Column(SMALLINT)
    AngularResolution = Column(SMALLINT)
    SpatialResolution = Column(SMALLINT)
    PatternType = Column(TINYINT)
    PatternDate = Column(DATETIME)
    PatternResolution = Column(SMALLINT)
    TransmitCenterFreqMHz = Column(DOUBLE)
    DopplerResolutionHzPerBin = Column(DOUBLE)
    FirstOrderMethod = Column(TINYINT)
    BraggSmoothingPoints = Column(TINYINT)
    CurrentVelocityLimit = Column(SMALLINT)
    BraggHasSecondOrder = Column(TINYINT)
    RadialBraggPeakDropOff = Column(DOUBLE)
    RadialBraggPeakNull = Column(DOUBLE)
    RadialBraggNoiseThreshold = Column(DOUBLE)
    PatternAmplitudeCorrections = Column(VARCHAR)
    PatternPhaseCorrections = Column(VARCHAR)
    PatternAmplitudeCalculations = Column(VARCHAR)
    PatternPhaseCalculations = Column(VARCHAR)
    RadialMusicParameters = Column(VARCHAR)
    MergedCount = Column(TINYINT)
    RadialMinimumMergePoints = Column(TINYINT)
    FirstOrderCalc = Column(TINYINT)
    MergeMethod = Column(VARCHAR)
    PatternMethod = Column(VARCHAR)
    TransmitSweepRateHz = Column(FLOAT)
    TransmitBandwidthKHz = Column(FLOAT)
    SpectraRangeCells = Column(TINYINT)
    SpectraDopplerCells = Column(SMALLINT)
    TableType = Column(VARCHAR)
    TableColumns = Column(TINYINT)
    TableColumnTypes = Column(VARCHAR)
    TableRows = Column(SMALLINT)
    fileModTime = Column(DATETIME)
    # processingMethod = Column(TINYINT)
    dateAdded = Column(TIMESTAMP)


class Sites(Base):
    __tablename__ = "hfrSites"

    id = Column(INTEGER, primary_key=True)
    site = Column(VARCHAR)
    type = Column(TINYINT)
    transmitCenterFrequency = Column(DOUBLE)
    description = Column(VARCHAR)
    lat = Column(DOUBLE)
    lon = Column(DOUBLE)
    # active = Column(TINYINT)
    # numSites = Column(TINYINT)
    # radialProcessingType = Column(INTEGER)
    # useRadials = Column(TINYINT)


class SystemTypes(Base):
    __tablename__ = 'hfrSystemTypes'

    id = Column(INTEGER, primary_key=True)
    frequency_min = Column(FLOAT)
    frequency_max = Column(FLOAT)


class WaveData(Base):
    __tablename__ = "hfrWaveData"

    id = Column(INTEGER, primary_key=True)
    site_id = Column(SMALLINT)
    TIME = Column(INTEGER)
    datetime = Column(DATETIME)
    MWHT = Column(FLOAT)
    MWPD = Column(FLOAT)
    WAVB = Column(FLOAT)
    WNDB = Column(FLOAT)
    PMWH = Column(FLOAT)
    ACNT = Column(TINYINT)
    DIST = Column(FLOAT)
    RCLL = Column(TINYINT)
    WDPT = Column(TINYINT)
    MTHD = Column(TINYINT)
    FLAG = Column(TINYINT)
    WHNM = Column(SMALLINT)
    WHSD = Column(FLOAT)
    TYRS = Column(SMALLINT)
    TMON = Column(TINYINT)
    TDAY = Column(TINYINT)
    THRS = Column(TINYINT)
    TMIN = Column(TINYINT)
    TSEC = Column(TINYINT)
    file_id = Column(SMALLINT)
    mwht_flag = Column(TINYINT)


class WaveFileHeader(Base):
    __tablename__ = "hfrWaveFilesMetadata"

    id = Column(INTEGER, primary_key=True)
    filename = Column(VARCHAR)
    CTF = Column(FLOAT)
    FileType = Column(VARCHAR)
    UUID = Column(VARCHAR)
    Manufacturer = Column(VARCHAR)
    Site = Column(INTEGER)
    TimeStamp = Column(TIMESTAMP)
    TimeZone = Column(CHAR)
    TimeCoverage = Column(SMALLINT)
    Origin = Column(VARCHAR)
    RangeCells = Column(SMALLINT)
    RangeResolutionKMeters = Column(FLOAT)
    AntennaBearing = Column(SMALLINT)
    TransmitCenterFreqMHz = Column(FLOAT)
    BraggSmoothingPoints = Column(SMALLINT)
    CurrentVelocityLimit = Column(SMALLINT)
    TransmitSweepRateHz = Column(FLOAT)
    TransmitBandwidthKHz = Column(FLOAT)
    CoastlineSector = Column(VARCHAR)
    DopplerCells = Column(FLOAT)
    MaximumWavePeriod = Column(FLOAT)
    WaveBearingLimits = Column(VARCHAR)
    WaveBraggNoiseThreshold = Column(FLOAT)
    WaveBraggPeakDropOff = Column(FLOAT)
    WaveBraggPeakNull = Column(FLOAT)
    WaveMergeMethod = Column(SMALLINT)
    WaveMinDopplerPoints = Column(SMALLINT)
    WaveUseInnerBragg = Column(SMALLINT)
    WavesFollowTheWind = Column(SMALLINT)
    BraggHasSecondOrder = Column(TINYINT)
    TableWaveMode = Column(TINYINT)
