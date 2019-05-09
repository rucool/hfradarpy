# ************************************************************
# Sequel Pro SQL dump
# Version 4541
#
# http://www.sequelpro.com/
# https://github.com/sequelpro/sequelpro
#
# Host: 127.0.0.1 (MySQL 5.7.20)
# Database: coolops
# Generation Time: 2018-06-08 19:14:53 +0000
# ************************************************************


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;


# Dump of table hfrChecks
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrChecks`;

CREATE TABLE `hfrChecks` (
  `diag_id` mediumint(8) unsigned NOT NULL AUTO_INCREMENT,
  `diag_site_id` smallint(6) unsigned NOT NULL,
  `freq` double unsigned NOT NULL,
  `time_chk` binary(1) NOT NULL DEFAULT '0',
  `rx_temp` tinyint(2) unsigned NOT NULL,
  `fwd_power` tinyint(3) unsigned NOT NULL,
  `refl_power` tinyint(2) NOT NULL,
  `disk_space_avail` smallint(4) unsigned NOT NULL,
  `plot_amps` binary(1) NOT NULL DEFAULT '0',
  `seaphase_loop1` smallint(4) NOT NULL,
  `seaphase_loop2` smallint(4) NOT NULL,
  `plot_signal` binary(1) NOT NULL DEFAULT '0',
  `plot_rx_temp` binary(1) NOT NULL DEFAULT '0',
  `plot_power` binary(1) NOT NULL DEFAULT '0',
  `archive` binary(1) NOT NULL DEFAULT '0',
  `diag_note` varchar(500) CHARACTER SET latin1 COLLATE latin1_general_ci NOT NULL,
  `diag_op_id` smallint(6) NOT NULL,
  `diag_creation_date` datetime NOT NULL,
  `diag_edit_date` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00' ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`diag_id`),
  KEY `diag_op_id` (`diag_op_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



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


# Dump of table hfrHardwareDiagnostics
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrHardwareDiagnostics`;

CREATE TABLE `hfrHardwareDiagnostics` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `id_site` int(11) DEFAULT NULL,
  `id_radial` int(11) DEFAULT NULL,
  `RTMP` smallint(3) DEFAULT NULL,
  `MTMP` smallint(3) DEFAULT NULL,
  `XTRP` varchar(2) DEFAULT NULL,
  `RUNT` int(8) DEFAULT NULL,
  `SP24` double(3,2) DEFAULT NULL,
  `SP05` double(3,2) DEFAULT NULL,
  `SN05` double(3,2) DEFAULT NULL,
  `SP12` double(4,2) DEFAULT NULL,
  `XPHT` tinyint(2) DEFAULT NULL,
  `XAHT` tinyint(2) DEFAULT NULL,
  `XAFW` tinyint(2) DEFAULT NULL,
  `XARW` tinyint(2) DEFAULT NULL,
  `X1VR` double(2,1) DEFAULT NULL,
  `XP28` double(4,2) DEFAULT NULL,
  `XP05` double(4,2) DEFAULT NULL,
  `GRMD` tinyint(1) DEFAULT NULL,
  `GDMD` tinyint(1) DEFAULT NULL,
  `GSLK` tinyint(1) DEFAULT NULL,
  `GSUL` tinyint(1) DEFAULT NULL,
  `PLLL` tinyint(1) DEFAULT NULL,
  `HTMP` double(4,1) DEFAULT NULL,
  `HUMI` tinyint(3) DEFAULT NULL,
  `RBIA` double(3,2) DEFAULT NULL,
  `EXTA` tinyint(1) DEFAULT NULL,
  `EXTB` tinyint(1) DEFAULT NULL,
  `CRUN` double(7,2) DEFAULT NULL,
  `datetime` datetime DEFAULT NULL,
  `upload_timestamp` timestamp NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table hfrInstitutions
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrInstitutions`;

CREATE TABLE `hfrInstitutions` (
  `id` tinyint(3) unsigned NOT NULL AUTO_INCREMENT,
  `name` tinytext NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



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



# Dump of table hfrNotify
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrNotify`;

CREATE TABLE `hfrNotify` (
  `nt_id` mediumint(10) unsigned NOT NULL AUTO_INCREMENT,
  `nt_op` smallint(3) unsigned NOT NULL,
  `nt_region` smallint(3) unsigned NOT NULL,
  `nt_inst` smallint(3) unsigned NOT NULL,
  PRIMARY KEY (`nt_id`),
  KEY `nt_op_id` (`nt_op`),
  KEY `nt_site_id` (`nt_inst`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table hfrOperators
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrOperators`;

CREATE TABLE `hfrOperators` (
  `op_id` smallint(5) unsigned NOT NULL AUTO_INCREMENT,
  `op_username` tinytext NOT NULL,
  `op_pass` char(60) NOT NULL,
  `op_salt` char(22) NOT NULL,
  `op_email` varchar(50) NOT NULL,
  `op_firstname` varchar(25) NOT NULL,
  `op_lastname` varchar(25) NOT NULL,
  `op_initials` tinytext NOT NULL,
  `op_inst` tinyint(3) unsigned NOT NULL,
  `op_region_id` tinyint(3) unsigned NOT NULL,
  `op_role` varchar(15) NOT NULL,
  `op_status` varchar(15) NOT NULL,
  PRIMARY KEY (`op_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table hfrOutages
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrOutages`;

CREATE TABLE `hfrOutages` (
  `out_id` smallint(5) unsigned NOT NULL AUTO_INCREMENT,
  `out_site_id` smallint(5) unsigned NOT NULL DEFAULT '0',
  `out_code` smallint(3) unsigned NOT NULL DEFAULT '999',
  `out_cause` varchar(250) CHARACTER SET latin1 COLLATE latin1_general_ci NOT NULL,
  `out_desc` text CHARACTER SET latin1 COLLATE latin1_general_ci,
  `out_est_repair_date` date DEFAULT '0000-00-00',
  `out_status` binary(1) NOT NULL DEFAULT '1' COMMENT '1 if current outage, 0 if outage has been resolved',
  `out_op_id` tinyint(3) unsigned NOT NULL,
  `out_start_date` datetime DEFAULT NULL,
  `out_resolved_date` datetime DEFAULT NULL,
  `out_data_lost` binary(1) NOT NULL DEFAULT '1',
  `out_edit_date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`out_id`),
  KEY `out_op_id` (`out_op_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table hfrOutageTypes
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrOutageTypes`;

CREATE TABLE `hfrOutageTypes` (
  `outtype_id` tinyint(3) unsigned NOT NULL AUTO_INCREMENT,
  `outtype_desc` varchar(50) CHARACTER SET latin1 COLLATE latin1_general_ci NOT NULL,
  `outtype_code` smallint(3) unsigned NOT NULL,
  PRIMARY KEY (`outtype_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

LOCK TABLES `hfrOutageTypes` WRITE;
/*!40000 ALTER TABLE `hfrOutageTypes` DISABLE KEYS */;

INSERT INTO `hfrOutageTypes` (`outtype_id`, `outtype_desc`, `outtype_code`)
VALUES
	(1,'unknown',999),
	(2,'hardware',100),
	(3,'software',200),
	(4,'communications',300),
	(5,'power',500),
	(6,'operations',400);

/*!40000 ALTER TABLE `hfrOutageTypes` ENABLE KEYS */;
UNLOCK TABLES;


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
	(1,'Ideal','2018-05-15 10:54:58'),
	(2,'Measured','2018-05-15 10:55:02');

/*!40000 ALTER TABLE `hfrPatternTypes` ENABLE KEYS */;
UNLOCK TABLES;


# Dump of table hfrQCValues
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrQCValues`;

CREATE TABLE `hfrQCValues` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `radial_min_count` smallint(4) unsigned NOT NULL DEFAULT '25' COMMENT 'failure threshold for low radial count test',
  `radial_low_count` smallint(4) unsigned NOT NULL DEFAULT '150',
  `radial_max_speed` smallint(3) unsigned NOT NULL DEFAULT '300',
  PRIMARY KEY (`id`),
  CONSTRAINT `site_id` FOREIGN KEY (`id`) REFERENCES `hfrSites` (`id`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table hfrRadialDiagnostics
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrRadialDiagnostics`;

CREATE TABLE `hfrRadialDiagnostics` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `id_site` int(11) DEFAULT NULL,
  `id_radial` int(11) DEFAULT NULL,
  `AMP1` double(5,4) DEFAULT NULL,
  `AMP2` double(5,4) DEFAULT NULL,
  `PH13` double(4,1) DEFAULT NULL,
  `PH23` double(4,1) DEFAULT NULL,
  `CPH1` double(4,1) DEFAULT NULL,
  `CPH2` double(4,1) DEFAULT NULL,
  `SNF1` double(4,0) DEFAULT NULL,
  `SNF2` double(4,0) DEFAULT NULL,
  `SNF3` double(4,0) DEFAULT NULL,
  `SSN1` double(4,0) DEFAULT NULL,
  `SSN2` double(4,0) DEFAULT NULL,
  `SSN3` double(4,0) DEFAULT NULL,
  `DGRC` tinyint(1) DEFAULT NULL,
  `DOPV` smallint(4) DEFAULT NULL,
  `DDAP` tinyint(2) DEFAULT NULL,
  `RADV` smallint(4) DEFAULT NULL,
  `RAPR` tinyint(2) DEFAULT NULL,
  `RARC` tinyint(2) DEFAULT NULL,
  `RADR` double(4,1) DEFAULT NULL,
  `RMCV` double(4,1) DEFAULT NULL,
  `RACV` double(4,1) DEFAULT NULL,
  `RABA` double(4,1) DEFAULT NULL,
  `RTYP` tinyint(1) DEFAULT NULL,
  `STYP` tinyint(2) DEFAULT NULL,
  `datetime` datetime DEFAULT NULL,
  `upload_timestamp` timestamp NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



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
  `AngularResolution` smallint(6) DEFAULT NULL,
  `SpatialResolution` smallint(6) DEFAULT NULL,
  `PatternType` tinyint(4) DEFAULT '0',
  `PatternDate` datetime DEFAULT NULL,
  `PatternResolution` smallint(6) DEFAULT NULL,
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
  `dateAdded` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
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



# Dump of table hfrSites
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrSites`;

CREATE TABLE `hfrSites` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `site` varchar(8) NOT NULL,
  `type` tinyint(4) DEFAULT NULL COMMENT 'see hfrSystemTypes for meaning of codes',
  `transmitCenterFrequency` double DEFAULT NULL,
  `description` varchar(256) DEFAULT NULL,
  `lat` double DEFAULT NULL,
  `lon` double DEFAULT NULL,
  `dateAdded` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `active` tinyint(3) unsigned NOT NULL DEFAULT '1' COMMENT 'see hfrSiteStatus for codes',
  `webServer` varchar(128) DEFAULT NULL,
  `institution` tinyint(3) unsigned NOT NULL,
  `region` tinyint(3) unsigned DEFAULT '1',
  `radialProcessingType` int(10) unsigned NOT NULL DEFAULT '2' COMMENT 'see hfrPatternTypes table for meaning of codes',
  `useRadials` tinyint(3) unsigned NOT NULL DEFAULT '0' COMMENT 'this field is not currently in use',
  `low_count` smallint(4) unsigned NOT NULL DEFAULT '25' COMMENT 'failure threshold for low radial count test',
  `operator` smallint(5) unsigned DEFAULT NULL,
  `application` smallint(1) DEFAULT '0',
  PRIMARY KEY (`id`),
  UNIQUE KEY `site` (`site`),
  CONSTRAINT `id` FOREIGN KEY (`id`) REFERENCES `hfrQCValues` (`id`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


DELIMITER ;;
/*!50003 SET SESSION SQL_MODE="NO_AUTO_CREATE_USER,NO_ENGINE_SUBSTITUTION" */;;
/*!50003 CREATE */ /*!50017 DEFINER=`coolops`@`%` */ /*!50003 TRIGGER `after_hfrSites_update` AFTER UPDATE ON `hfrSites` FOR EACH ROW BEGIN
   IF NEW.active <> OLD.active THEN
   INSERT INTO hfrSitesAudit
   SET site_id = OLD.id,
    field = 'active',
    old_value = OLD.active,
    new_value = NEW.active,
    modified_date = NOW();
   END IF;   
   IF NEW.radialProcessingType <> OLD.radialProcessingType THEN
   INSERT INTO hfrSitesAudit
   SET site_id = OLD.id,
    field = 'radialProcessingType',
    old_value = OLD.radialProcessingType,
    new_value = NEW.radialProcessingType,
    modified_date = NOW();
    END IF;
 END */;;
DELIMITER ;
/*!50003 SET SESSION SQL_MODE=@OLD_SQL_MODE */;


# Dump of table hfrSitesAudit
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrSitesAudit`;

CREATE TABLE `hfrSitesAudit` (
  `id` int(8) unsigned NOT NULL AUTO_INCREMENT,
  `site_id` mediumint(8) unsigned NOT NULL,
  `field` tinytext NOT NULL,
  `old_value` tinyint(3) unsigned NOT NULL,
  `new_value` tinyint(3) unsigned NOT NULL,
  `modified_date` datetime DEFAULT NULL,
  `comment` varchar(200) DEFAULT NULL,
  `operator_id` tinyint(3) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`)
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
	(0,'Inactive'),
	(1,'Active'),
	(2,'Retired');

/*!40000 ALTER TABLE `hfrSiteStatus` ENABLE KEYS */;
UNLOCK TABLES;


# Dump of table hfrSystemTypes
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrSystemTypes`;

CREATE TABLE `hfrSystemTypes` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` char(2) NOT NULL,
  `description` varchar(128) DEFAULT NULL,
  `dateAdded` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `avg_dist` smallint(2) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `type` (`type`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

LOCK TABLES `hfrSystemTypes` WRITE;
/*!40000 ALTER TABLE `hfrSystemTypes` DISABLE KEYS */;

INSERT INTO `hfrSystemTypes` (`id`, `type`, `description`, `dateAdded`, `avg_dist`)
VALUES
	(2,'lr','Long Range','2017-05-19 14:26:33',12),
	(3,'mr','Medium Range','2017-05-19 14:26:37',6),
	(4,'sr','Standard Range','2017-05-19 14:26:39',2);

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



# Dump of table hfrTotalsConfig
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrTotalsConfig`;

CREATE TABLE `hfrTotalsConfig` (
  `totals_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `setup_id` smallint(5) unsigned NOT NULL,
  `site_id` smallint(5) unsigned NOT NULL,
  `use_radials` tinyint(3) unsigned NOT NULL,
  `radial_type` tinyint(3) unsigned NOT NULL,
  PRIMARY KEY (`totals_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table hfrTotalsSetup
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrTotalsSetup`;

CREATE TABLE `hfrTotalsSetup` (
  `setup_id` smallint(5) unsigned NOT NULL AUTO_INCREMENT,
  `system_type` smallint(5) unsigned NOT NULL,
  `region_id` tinyint(3) unsigned NOT NULL,
  `effective_date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`setup_id`)
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
	(1,'No'),
	(2,'Yes');

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
  `ACNT` tinyint(1) DEFAULT NULL,
  `DIST` float(6,4) DEFAULT NULL,
  `RCLL` tinyint(2) DEFAULT NULL,
  `WDPT` tinyint(3) DEFAULT NULL,
  `MTHD` tinyint(1) DEFAULT NULL,
  `FLAG` smallint(3) DEFAULT NULL,
  `WHNM` smallint(3) DEFAULT NULL,
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
  KEY `site_id` (`site_id`),
  KEY `TYRS` (`TYRS`),
  KEY `TMON` (`TMON`),
  KEY `TDAY` (`TDAY`),
  KEY `THRS` (`THRS`),
  KEY `TMIN` (`TMIN`),
  KEY `TSEC` (`TSEC`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table hfrWaveFilesMetadata
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hfrWaveFilesMetadata`;

CREATE TABLE `hfrWaveFilesMetadata` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `filename` varchar(64) NOT NULL,
  `CTF` double DEFAULT NULL,
  `FileType` varchar(64) DEFAULT NULL,
  `UUID` varchar(64) DEFAULT NULL,
  `Manufacturer` varchar(64) DEFAULT NULL,
  `Site` int(10) NOT NULL,
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
  `TableRows` smallint(5) unsigned DEFAULT NULL,
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
  `date_updated` timestamp NULL DEFAULT NULL ON UPDATE CURRENT_TIMESTAMP,
  `date_added` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
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
	(1,'Average'),
	(2,'Median'),
	(3,'None');

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
  `dist_units` varchar(255) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;




/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;
/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
