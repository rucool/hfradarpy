import datetime as dt
import re


def clean_wave_header(head_dict):
    """
    Cleans the header data from the wave data for proper input into MySQL database
    :param head_dict: dictionary containing the header data
    :return: dictionary containing the cleaned header data
    :rtype:
    """

    keep = ['TimeCoverage', 'WaveMinDopplerPoints', 'AntennaBearing','DopplerCells','TransmitCenterFreqMHz',
            'CTF', 'TableColumnTypes', 'TimeZone', 'WaveBraggPeakDropOff','RangeResolutionKMeters',
            'CoastlineSector', 'WaveMergeMethod', 'RangeCells', 'WaveBraggPeakNull', 'WaveUseInnerBragg',
            'BraggSmoothingPoints', 'Manufacturer', 'TimeStamp', 'FileType', 'TableRows', 'BraggHasSecondOrder',
            'Origin', 'MaximumWavePeriod', 'UUID', 'WaveBraggNoiseThreshold', 'TransmitBandwidthKHz', 'Site',
            'TransmitSweepRateHz', 'WaveBearingLimits', 'WavesFollowTheWind','CurrentVelocityLimit']

    for key in head_dict.keys():
        if not key in keep:
            del head_dict[key]

    for k, v in head_dict.iteritems():
        if 'Site' in k:
            head_dict[k] = ''.join(e for e in v if e.isalnum())
        elif 'TimeStamp' in k:
            t_list = v.split()
            t_list = [int(s) for s in t_list]
            head_dict[k] = dt.datetime(t_list[0], t_list[1], t_list[2], t_list[3], t_list[4], t_list[5]).strftime('%Y-%m-%d %H:%M:%S')
        elif k in ('TimeCoverage', 'RangeResolutionKMeters'):
            head_dict[k] = re.findall("\d+\.\d+", v)[0]
        elif k in ('WaveMergeMethod', 'WaveUseInnerBragg', 'WavesFollowTheWind'):
            head_dict[k] = re.search(r'\d+', v).group()
        elif 'TimeZone' in k:
            head_dict[k] = re.search('"(.*)"', v).group(1)
        elif k in ('WaveBearingLimits', 'CoastlineSector'):
            bearings = re.findall(r"[-+]?\d*\.\d+|\d+", v)
            head_dict[k] = ', '.join(e for e in bearings)
        else:
            continue
    return head_dict


def flag_bad_data(wave_height, wave_min=0.2, wave_max=5):
    if wave_min <= wave_height <= wave_max:
        flag = 1
    else:
        flag = 4
    return flag