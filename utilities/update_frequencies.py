#!/usr/bin/env python
import codar_processing.common as cf
import codar_processing.database_common as dbc
import os.path
from decimal import *


def update_frequency(session, site, frequency):
    """
    Check frequency of a site. It it changes, update it.
    :param session: sqlalchemy session binding
    :param site: four digit site code
    :return:
    """
    result = dbc.site_check(session, site)
    site_freq_db = result.transmitCenterFrequency
    center_freq = Decimal(frequency[1])

    if not site_freq_db == center_freq:
        print('{}: Updating frequency from {} MHz to {} MHz'.format(site, site_freq_db,  frequency[1]))
        result.transmitCenterFrequency = center_freq
        session.commit()
    else:
        print('{}: No frequency change required'.format(site))


global session
session = dbc.db_session()

sourceDir = '/home/codaradm/data/radials/'
site_codes = ['ASSA', 'BELM', 'BLCK', 'BRIG', 'BRMR', 'BRNT', 'CBBT', 'CDDO', 'CEDR', 'CMPT', 'DUCK', 'FURA', 'GCAP',
              'HATY', 'HLPN', 'HEMP', 'HOOK', 'LISL', 'LOVE', 'MISQ', 'MNTK', 'MRCH', 'MVCO', 'NANT', 'NAUS', 'PORT',
              'RATH', 'SILD', 'SLTR', 'STLI', 'SUNS', 'VIEW', 'WILD', 'WOOD']


for site in site_codes:
    siteDir = sourceDir + site + '/'
    logfiles = sorted([f for f in os.listdir(siteDir) if f.startswith('RDL')])

    if not logfiles:
        continue
    else:
        latestRadial = siteDir + logfiles[-1]
        radialFile = open(latestRadial)

        for line in radialFile:
            if line.startswith("%"):
                if 'TransmitCenterFreqMHz' in line:
                    freq = cf.parse_header_line(line)
                    update_frequency(session, site, freq)