# ************************************************************
# Sequel Pro SQL dump
# Version 4541
#
# http://www.sequelpro.com/
# https://github.com/sequelpro/sequelpro
#
# Host: 127.0.0.1 (MySQL 5.7.20)
# Database: coolops
# Generation Time: 2018-03-23 15:18:17 +0000
# ************************************************************


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;


# Dump of table hfrFileTypes
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrFileTypes`;

CREATE TABLE `hfrFileTypes` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` char(4) NOT NULL,
  `description` varchar(256) DEFAULT NULL,
  `dateAdded` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `type` (`type`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

LOCK TABLES `hfrFileTypes` WRITE;
/*!40000 ALTER TABLE `hfrFileTypes` DISABLE KEYS */;

INSERT INTO `hfrFileTypes` (`id`, `type`, `description`, `dateAdded`)
VALUES
	(2,'RDLi','Radials from Ideal Antenna Pattern','2011-12-14 15:24:18'),
	(3,'RDLm','Radials from Measured Antenna Pattern','2011-12-14 15:24:18'),
	(4,'RDLx','Short-time Radials from Ideal Antenna Pattern','2011-12-14 15:24:18'),
	(5,'RDLy','Short-time Radials from Measured Antenna Pattern','2011-12-14 15:24:18'),
	(6,'ELTi','Ellipiticals from Ideal Antenna Pattern','2011-12-14 15:24:18'),
	(7,'ELTm','Ellipticals from Measured Antenna Pattern','2011-12-14 15:24:18'),
	(8,'TOTL','Total Current Vectors','2011-12-14 15:24:18');

/*!40000 ALTER TABLE `hfrFileTypes` ENABLE KEYS */;
UNLOCK TABLES;


# Dump of table hfrLatestRadials
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrLatestRadials`;

CREATE TABLE `hfrLatestRadials` (
  `siteId` smallint(5) unsigned NOT NULL,
  `radialId` int(10) unsigned NOT NULL,
  `TimeStamp` datetime NOT NULL,
  `filename` varchar(64) NOT NULL,
  `dateAddedEasternTime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`siteId`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table hfrPatternTypes
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrPatternTypes`;

CREATE TABLE `hfrPatternTypes` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` varchar(32) NOT NULL,
  `dateAdded` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `type` (`type`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

LOCK TABLES `hfrPatternTypes` WRITE;
/*!40000 ALTER TABLE `hfrPatternTypes` DISABLE KEYS */;

INSERT INTO `hfrPatternTypes` (`id`, `type`, `dateAdded`)
VALUES
	(2,'Ideal','2011-12-14 12:32:13'),
	(3,'Measured','2011-12-14 12:32:23');

/*!40000 ALTER TABLE `hfrPatternTypes` ENABLE KEYS */;
UNLOCK TABLES;


# Dump of table hfrRadialFilesMetadata
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrRadialFilesMetadata`;

CREATE TABLE `hfrRadialFilesMetadata` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `filename` varchar(64) NOT NULL,
  `CTF` double DEFAULT NULL,
  `FileType` varchar(64) DEFAULT NULL,
  `LLUVSpec` varchar(64) DEFAULT NULL,
  `UUID` varchar(64) DEFAULT NULL,
  `Manufacturer` varchar(64) DEFAULT NULL,
  `Site` int(10) NOT NULL,
  `TimeStamp` datetime NOT NULL,
  `TimeZone` char(3) NOT NULL,
  `TimeCoverage` smallint(6) NOT NULL,
  `Origin` varchar(64) DEFAULT NULL,
  `GreatCircle` varchar(64) DEFAULT NULL,
  `GeodVersion` varchar(64) DEFAULT NULL,
  `LLUVTrustData` varchar(64) DEFAULT NULL,
  `RangeStart` tinyint(3) unsigned DEFAULT NULL,
  `RangeEnd` tinyint(3) unsigned DEFAULT NULL,
  `RangeResolutionKMeters` float NOT NULL,
  `AntennaBearing` smallint(6) NOT NULL,
  `ReferenceBearing` smallint(6) NOT NULL,
  `AngularResolution` smallint(6) NOT NULL,
  `SpatialResolution` smallint(6) NOT NULL,
  `PatternType` tinyint(4) DEFAULT '0',
  `PatternDate` datetime DEFAULT NULL,
  `PatternResolution` smallint(6) NOT NULL,
  `TransmitCenterFreqMHz` double DEFAULT NULL,
  `DopplerResolutionHzPerBin` double DEFAULT NULL,
  `FirstOrderMethod` tinyint(3) unsigned DEFAULT NULL,
  `BraggSmoothingPoints` tinyint(3) unsigned DEFAULT NULL,
  `CurrentVelocityLimit` smallint(5) unsigned DEFAULT NULL,
  `BraggHasSecondOrder` tinyint(3) unsigned DEFAULT NULL,
  `RadialBraggPeakDropOff` double DEFAULT NULL,
  `RadialBraggPeakNull` double DEFAULT NULL,
  `RadialBraggNoiseThreshold` double DEFAULT NULL,
  `PatternAmplitudeCorrections` varchar(64) DEFAULT NULL,
  `PatternPhaseCorrections` varchar(64) DEFAULT NULL,
  `PatternAmplitudeCalculations` varchar(64) DEFAULT NULL,
  `PatternPhaseCalculations` varchar(64) DEFAULT NULL,
  `RadialMusicParameters` varchar(64) DEFAULT NULL,
  `MergedCount` tinyint(3) unsigned DEFAULT NULL,
  `RadialMinimumMergePoints` tinyint(3) unsigned DEFAULT NULL,
  `FirstOrderCalc` tinyint(3) unsigned DEFAULT NULL,
  `MergeMethod` varchar(64) DEFAULT NULL,
  `PatternMethod` varchar(64) DEFAULT NULL,
  `TransmitSweepRateHz` float DEFAULT NULL,
  `TransmitBandwidthKHz` float DEFAULT NULL,
  `SpectraRangeCells` tinyint(3) unsigned DEFAULT NULL,
  `SpectraDopplerCells` smallint(5) unsigned DEFAULT NULL,
  `TableType` varchar(64) DEFAULT NULL,
  `TableColumns` tinyint(4) DEFAULT NULL,
  `TableColumnTypes` varchar(255) DEFAULT NULL,
  `TableRows` smallint(5) unsigned DEFAULT '0',
  `fileModTime` datetime NOT NULL,
  `processingMethod` tinyint(3) unsigned NOT NULL,
  `dateAdded` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `radialFileType` (`TimeStamp`,`Site`,`processingMethod`),
  KEY `ts` (`TimeStamp`),
  KEY `sid` (`Site`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table hfrRegions
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrRegions`;

CREATE TABLE `hfrRegions` (
  `id` tinyint(3) unsigned NOT NULL AUTO_INCREMENT,
  `ra` varchar(16) DEFAULT NULL,
  `url` varchar(128) DEFAULT NULL,
  `date_added` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

LOCK TABLES `hfrRegions` WRITE;
/*!40000 ALTER TABLE `hfrRegions` DISABLE KEYS */;

INSERT INTO `hfrRegions` (`id`, `ra`, `url`, `date_added`)
VALUES
	(1,'Unspecified','','2012-05-21 16:18:31'),
	(2,'AOOS','http://www.aoos.org','2012-05-20 22:56:30'),
	(3,'CaraCOOS','http://cara.uprm.edu','2012-05-20 22:56:30'),
	(4,'CenCOOS','http://www.cencoos.org','2012-05-20 22:56:30'),
	(5,'GCOOS','http://gcoos.org','2012-05-20 22:56:30'),
	(6,'GLOS','http://glos.us','2012-05-20 22:56:30'),
	(7,'MARACOOS','http://maracoos.org','2012-05-20 22:56:30'),
	(8,'NANOOS','http://www.nanoos.org','2012-05-20 22:56:30'),
	(9,'NERACOOS','http://www.neracoos.org','2012-05-20 22:56:30'),
	(10,'PacIOOS','http://www.pacioos.org','2012-05-20 22:56:30'),
	(11,'SCCOOS','http://www.scoos.org','2012-05-20 22:56:30'),
	(12,'SECOORA','http://secoora.org','2012-05-20 22:56:30'),
	(13,'Palmer Deep',NULL,'2014-11-18 15:43:04');

/*!40000 ALTER TABLE `hfrRegions` ENABLE KEYS */;
UNLOCK TABLES;


# Dump of table hfrSites
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrSites`;

CREATE TABLE `hfrSites` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `site` varchar(8) NOT NULL,
  `type` tinyint(4) DEFAULT NULL,
  `transmitCenterFrequency` double DEFAULT NULL,
  `description` varchar(256) DEFAULT NULL,
  `lat` double DEFAULT NULL,
  `lon` double DEFAULT NULL,
  `dateAdded` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `active` tinyint(3) unsigned NOT NULL DEFAULT '1',
  `webServer` varchar(128) DEFAULT NULL,
  `numSites` tinyint(3) unsigned DEFAULT '1',
  `region` tinyint(3) unsigned DEFAULT '1',
  `radialProcessingType` int(10) unsigned NOT NULL DEFAULT '2',
  `useRadials` tinyint(3) unsigned NOT NULL DEFAULT '0',
  `test_radialProcessingType` int(10) unsigned NOT NULL DEFAULT '2' COMMENT 'for testing only (added by Teresa Updyke)',
  PRIMARY KEY (`id`),
  UNIQUE KEY `site` (`site`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table hfrSiteStatus
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrSiteStatus`;

CREATE TABLE `hfrSiteStatus` (
  `statusId` int(1) NOT NULL,
  `statusType` varchar(10) NOT NULL,
  UNIQUE KEY `statusId` (`statusId`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

LOCK TABLES `hfrSiteStatus` WRITE;
/*!40000 ALTER TABLE `hfrSiteStatus` DISABLE KEYS */;

INSERT INTO `hfrSiteStatus` (`statusId`, `statusType`)
VALUES
	(0,'Retired'),
	(1,'Active'),
	(2,'Inactive');

/*!40000 ALTER TABLE `hfrSiteStatus` ENABLE KEYS */;
UNLOCK TABLES;


# Dump of table hfrSystemTypes
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrSystemTypes`;

CREATE TABLE `hfrSystemTypes` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` char(2) DEFAULT '',
  `description` varchar(128) DEFAULT NULL,
  `dateAdded` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `avg_dist` smallint(2) DEFAULT NULL,
  `frequency_min` float DEFAULT NULL,
  `frequency_max` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `type` (`type`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

LOCK TABLES `hfrSystemTypes` WRITE;
/*!40000 ALTER TABLE `hfrSystemTypes` DISABLE KEYS */;

INSERT INTO `hfrSystemTypes` (`id`, `type`, `description`, `dateAdded`, `avg_dist`, `frequency_min`, `frequency_max`)
VALUES
	(2,'lr','Long Range','2017-06-01 16:24:17',12,4.3,5.4),
	(3,'mr','Medium Range','2017-06-01 16:24:32',6,11.5,14),
	(4,'sr','Standard Range','2017-06-01 16:24:41',2,24,27);

/*!40000 ALTER TABLE `hfrSystemTypes` ENABLE KEYS */;
UNLOCK TABLES;


# Dump of table hfrTotals
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrTotals`;

CREATE TABLE `hfrTotals` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `systemTypeId` int(10) unsigned NOT NULL,
  `radialTimestamp` datetime NOT NULL,
  `numRadials` int(1) unsigned NOT NULL,
  `dateProcessed` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `region` tinyint(2) NOT NULL,
  `pattTypes` tinyint(4) NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table hfrWaveAffirmation
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrWaveAffirmation`;

CREATE TABLE `hfrWaveAffirmation` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `affirmation` tinytext,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

LOCK TABLES `hfrWaveAffirmation` WRITE;
/*!40000 ALTER TABLE `hfrWaveAffirmation` DISABLE KEYS */;

INSERT INTO `hfrWaveAffirmation` (`id`, `affirmation`)
VALUES
	(0,'No'),
	(1,'Yes');

/*!40000 ALTER TABLE `hfrWaveAffirmation` ENABLE KEYS */;
UNLOCK TABLES;


# Dump of table hfrWaveData
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrWaveData`;

CREATE TABLE `hfrWaveData` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `site_id` smallint(5) unsigned NOT NULL,
  `datetime` datetime DEFAULT NULL,
  `upload_date` timestamp NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `TIME` int(7) DEFAULT NULL,
  `MWHT` float(5,2) DEFAULT NULL,
  `MWPD` float(5,2) DEFAULT NULL,
  `WAVB` float(6,2) DEFAULT NULL,
  `WNDB` float(6,2) DEFAULT NULL,
  `PMWH` float(6,2) DEFAULT NULL,
  `ACNT` tinyint(2) DEFAULT NULL,
  `DIST` float(6,4) DEFAULT NULL,
  `RCLL` tinyint(2) DEFAULT NULL,
  `WDPT` tinyint(3) DEFAULT NULL,
  `MTHD` tinyint(1) DEFAULT NULL,
  `FLAG` int(6) DEFAULT NULL,
  `WHNM` smallint(2) DEFAULT NULL,
  `WHSD` float(3,2) DEFAULT NULL,
  `TYRS` smallint(4) DEFAULT NULL,
  `TMON` tinyint(2) DEFAULT NULL,
  `TDAY` tinyint(2) DEFAULT NULL,
  `THRS` tinyint(2) DEFAULT NULL,
  `TMIN` tinyint(2) DEFAULT NULL,
  `TSEC` tinyint(2) DEFAULT NULL,
  `file_id` int(10) DEFAULT NULL,
  `mwht_flag` tinyint(1) DEFAULT '1',
  PRIMARY KEY (`id`),
  KEY `file_id` (`file_id`),
  KEY `site_id` (`site_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table hfrWaveFileMode
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrWaveFileMode`;

CREATE TABLE `hfrWaveFileMode` (
  `id` tinyint(11) unsigned NOT NULL AUTO_INCREMENT,
  `mode` varchar(20) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

LOCK TABLES `hfrWaveFileMode` WRITE;
/*!40000 ALTER TABLE `hfrWaveFileMode` DISABLE KEYS */;

INSERT INTO `hfrWaveFileMode` (`id`, `mode`)
VALUES
	(1,'Classic'),
	(2,'Averaged');

/*!40000 ALTER TABLE `hfrWaveFileMode` ENABLE KEYS */;
UNLOCK TABLES;


# Dump of table hfrWaveFilesMetadata
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrWaveFilesMetadata`;

CREATE TABLE `hfrWaveFilesMetadata` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `filename` varchar(64) NOT NULL,
  `CTF` float DEFAULT NULL,
  `FileType` varchar(64) DEFAULT NULL,
  `UUID` varchar(64) DEFAULT NULL,
  `Manufacturer` varchar(64) DEFAULT NULL,
  `Site` smallint(10) NOT NULL,
  `TimeStamp` datetime NOT NULL,
  `TimeZone` char(3) NOT NULL,
  `TimeCoverage` smallint(6) NOT NULL,
  `Origin` varchar(64) DEFAULT NULL,
  `RangeCells` tinyint(3) unsigned DEFAULT NULL,
  `RangeResolutionKMeters` float DEFAULT NULL,
  `AntennaBearing` smallint(6) NOT NULL,
  `TransmitCenterFreqMHz` float DEFAULT NULL,
  `BraggSmoothingPoints` tinyint(3) unsigned DEFAULT NULL,
  `CurrentVelocityLimit` smallint(5) unsigned DEFAULT NULL,
  `TransmitSweepRateHz` float DEFAULT NULL,
  `TransmitBandwidthKHz` float DEFAULT NULL,
  `TableRows` smallint(5) unsigned DEFAULT '0',
  `date_updated` timestamp NULL DEFAULT NULL ON UPDATE CURRENT_TIMESTAMP,
  `date_added` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `CoastlineSector` varchar(12) DEFAULT NULL,
  `DopplerCells` float DEFAULT NULL,
  `MaximumWavePeriod` float DEFAULT NULL,
  `WaveBearingLimits` varchar(15) DEFAULT NULL,
  `WaveBraggNoiseThreshold` float DEFAULT NULL,
  `WaveBraggPeakDropOff` float DEFAULT NULL,
  `WaveBraggPeakNull` float DEFAULT NULL,
  `BraggHasSecondOrder` tinyint(1) unsigned DEFAULT NULL,
  `WaveMergeMethod` tinyint(1) unsigned DEFAULT NULL,
  `WaveMinDopplerPoints` tinyint(3) unsigned DEFAULT NULL,
  `WaveUseInnerBragg` tinyint(1) unsigned DEFAULT NULL,
  `WavesFollowTheWind` tinyint(1) unsigned DEFAULT NULL,
  `TableColumnTypes` varchar(255) DEFAULT NULL,
  `TableWaveMode` tinyint(11) unsigned NOT NULL DEFAULT '1',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table hfrWaveMergeMethod
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrWaveMergeMethod`;

CREATE TABLE `hfrWaveMergeMethod` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `method` tinytext,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

LOCK TABLES `hfrWaveMergeMethod` WRITE;
/*!40000 ALTER TABLE `hfrWaveMergeMethod` DISABLE KEYS */;

INSERT INTO `hfrWaveMergeMethod` (`id`, `method`)
VALUES
	(0,'None'),
	(1,'Average'),
	(2,'Median');

/*!40000 ALTER TABLE `hfrWaveMergeMethod` ENABLE KEYS */;
UNLOCK TABLES;


# Dump of table hfrWaveRange
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrWaveRange`;

CREATE TABLE `hfrWaveRange` (
  `site_id` int(10) DEFAULT NULL,
  `min_dist` double DEFAULT NULL,
  `max_dist` double DEFAULT NULL,
  `avg_dist` double DEFAULT NULL,
  `dist_units` varchar(255) NOT NULL DEFAULT 'Km'
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

LOCK TABLES `hfrWaveRange` WRITE;
/*!40000 ALTER TABLE `hfrWaveRange` DISABLE KEYS */;

INSERT INTO `hfrWaveRange` (`site_id`, `min_dist`, `max_dist`, `avg_dist`, `dist_units`)
VALUES
	(4,4.5,7.5,6,'Km'),
	(4,7.5,10.5,9,'Km'),
	(4,10.5,13.5,12,'Km'),
	(4,13.5,16.5,15,'Km'),
	(4,16.5,19.5,18,'Km'),
	(4,19.5,22.5,21,'Km'),
	(7,4.5,7.5,6,'Km'),
	(7,7.5,10.5,9,'Km'),
	(7,10.5,13.5,12,'Km'),
	(7,13.5,16.5,15,'Km'),
	(7,16.5,19.5,18,'Km'),
	(7,19.5,22.5,21,'Km'),
	(7,22.5,25.5,24,'Km'),
	(8,1.5,4.5,3,'Km'),
	(8,4.5,7.5,6,'Km'),
	(8,7.5,10.5,9,'Km'),
	(8,10.5,13.5,12,'Km'),
	(8,13.5,16.5,15,'Km'),
	(8,16.5,19.5,18,'Km'),
	(8,19.5,22.5,21,'Km'),
	(8,22.5,25.5,24,'Km'),
	(8,28.5,31.5,30,'Km'),
	(8,34.5,37.5,36,'Km'),
	(8,40.5,43.5,42,'Km'),
	(8,49.5,52.5,51,'Km'),
	(26,1.5,4.5,3,'Km'),
	(26,4.5,7.5,6,'Km'),
	(26,7.5,10.5,9,'Km'),
	(26,10.5,13.5,12,'Km'),
	(26,13.5,16.5,15,'Km'),
	(26,16.5,19.5,18,'Km'),
	(26,19.5,22.5,21,'Km'),
	(26,22.5,25.5,24,'Km'),
	(27,0.5,1.5,1,'Km'),
	(27,1.5,4.5,3,'Km'),
	(27,4.5,7.5,6,'Km'),
	(27,7.5,10.5,9,'Km'),
	(27,10.5,13.5,12,'Km'),
	(27,13.5,16.5,15,'Km'),
	(27,16.5,19.5,18,'Km'),
	(27,19.5,22.5,21,'Km'),
	(27,22.5,25.5,24,'Km'),
	(27,28.5,31.5,30,'Km'),
	(27,34.5,37.5,36,'Km'),
	(27,40.5,43.5,42,'Km'),
	(27,49.5,52.5,51,'Km'),
	(28,1.5,4.5,3,'Km'),
	(28,4.5,7.5,6,'Km'),
	(28,7.5,10.5,9,'Km'),
	(28,10.5,13.5,12,'Km'),
	(28,13.5,16.5,15,'Km'),
	(28,16.5,19.5,18,'Km'),
	(28,19.5,22.5,21,'Km'),
	(28,22.5,25.5,24,'Km'),
	(30,1.5,4.5,3,'Km'),
	(30,4.5,7.5,6,'Km'),
	(30,7.5,10.5,9,'Km'),
	(30,10.5,13.5,12,'Km'),
	(30,13.5,16.5,15,'Km'),
	(30,16.5,19.5,18,'Km'),
	(30,19.5,22.5,21,'Km'),
	(30,22.5,25.5,24,'Km'),
	(30,25.5,28.5,27,'Km'),
	(30,28.5,31.5,30,'Km'),
	(30,31.5,34.5,33,'Km'),
	(30,34.5,37.5,36,'Km'),
	(30,40.5,43.5,42,'Km'),
	(30,46.5,49.5,48,'Km'),
	(30,49.5,52.5,51,'Km'),
	(30,52.5,55.5,54,'Km'),
	(30,55.5,58.5,57,'Km'),
	(30,61.5,64.5,63,'Km'),
	(30,67.5,70.5,69,'Km'),
	(31,1.5,4.5,3,'Km'),
	(31,4.5,7.5,6,'Km'),
	(31,7.5,10.5,9,'Km'),
	(31,10.5,13.5,12,'Km'),
	(31,13.5,16.5,15,'Km'),
	(31,16.5,19.5,18,'Km'),
	(31,19.5,22.5,21,'Km'),
	(31,22.5,25.5,24,'Km'),
	(31,25.5,28.5,27,'Km'),
	(31,28.5,31.5,30,'Km'),
	(31,31.5,34.5,33,'Km'),
	(31,34.5,37.5,36,'Km'),
	(31,37.5,40.5,39,'Km'),
	(31,40.5,43.5,42,'Km'),
	(31,46.5,49.5,48,'Km'),
	(42,0.5,1.5,1,'Km'),
	(42,1.5,4.5,3,'Km'),
	(42,4.5,7.5,6,'Km'),
	(42,7.5,10.5,9,'Km'),
	(83,0.5,1.5,1,'Km'),
	(83,1.5,4.5,3,'Km'),
	(83,4.5,7.5,6,'Km'),
	(83,7.5,10.5,9,'Km');

/*!40000 ALTER TABLE `hfrWaveRange` ENABLE KEYS */;
UNLOCK TABLES;



/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;
/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
