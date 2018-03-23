from sqlalchemy import Column
from sqlalchemy.dialects.mysql import DOUBLE, TINYINT, TIMESTAMP, VARCHAR, FLOAT, CHAR, INTEGER, SMALLINT, DATETIME
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class WaveFile(Base):
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
    TableColumnTypes = Column(VARCHAR)
    TableRows = Column(SMALLINT)
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


class Sites(Base):
    __tablename__ = "hfrSites"

    id = Column(INTEGER, primary_key=True)
    site = Column(VARCHAR)
    lat = Column(DOUBLE)
    lon = Column(DOUBLE)
    transmitCenterFrequency = Column(DOUBLE)
    type = Column(TINYINT)


class hfrSystemTypes(Base):
    __tablename__ = 'hfrSystemTypes'

    id = Column(INTEGER, primary_key=True)
    frequency_min = Column(FLOAT)
    frequency_max = Column(FLOAT)


class hfrRadialFilesMetadata(Base):
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
    processingMethod = Column(TINYINT)
    dateAdded = Column(TIMESTAMP)


class hfrFileTypes(Base):
    __tablename__ = 'hfrFileTypes'

    id = Column(INTEGER, primary_key=True)
    type = Column(CHAR)
    description = Column(VARCHAR)


class hfrPatternTypes(Base):
    __tablename__ = 'hfrPatternTypes'

    id = Column(INTEGER, primary_key=True)
    type = Column(CHAR)