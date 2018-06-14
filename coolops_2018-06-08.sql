# ************************************************************
# Sequel Pro SQL dump
# Version 4541
#
# http://www.sequelpro.com/
# https://github.com/sequelpro/sequelpro
#
# Host: 127.0.0.1 (MySQL 5.7.20)
# Database: coolops
# Generation Time: 2018-06-08 19:59:31 +0000
# ************************************************************


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;


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




/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;
/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
