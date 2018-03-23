import datetime as dt


def clean_radial_header(head_dict):
    keep = ['CTF', 'FileType', 'LLUVSpec', 'UUID', 'Manufacturer', 'Site', 'TimeStamp', 'TimeZone', 'TimeCoverage',
            'Origin', 'GreatCircle', 'GeodVersion', 'LLUVTrustData', 'RangeStart', 'RangeEnd', 'RangeResolutionKMeters',
            'AntennaBearing', 'ReferenceBearing', 'AngularResolution', 'SpatialResolution', 'PatternType', 'PatternDate',
            'PatternResolution', 'TransmitCenterFreqMHz', 'DopplerResolutionHzPerBin', 'FirstOrderMethod',
            'BraggSmoothingPoints', 'CurrentVelocityLimits', 'BraggHasSecondOrder', 'RadialBraggPeakDropOff',
            'RadialBraggPeakNull', 'RadialBraggNoiseThreshold', 'PatternAmplitudeCorrections', 'PatternPhaseCorrections',
            'PatternAmplitudeCalculations', 'PatternPhaseCalculations', 'RadialMusicParameters', 'MergedCount',
            'RadialMinimumMergePoints', 'FirstOrderCalc', 'MergeMethod', 'PatternMethod', 'TransmitSweepRateHz',
            'TransmitBandwidthKHz', 'SpectraRangeCells', 'SpectraDopplerCells', 'TableType', 'TableColumns', 'TableColumnTypes',
            'TableRows', 'TableStart', 'CurrentVelocityLimit']
    # TableColumnTypes
    for key in head_dict.keys():
        if not key in keep:
            del head_dict[key]

    for k, v in head_dict.iteritems():
        if 'Site' in k:
            head_dict[k] = ''.join(e for e in v if e.isalnum())
        elif k in ('TimeStamp', 'PatternDate'):
            t_list = [int(s) for s in v.split()]
            head_dict[k] = dt.datetime(*t_list)
        elif 'TimeZone' in k:
            head_dict[k] = v.split('"')[1]
        elif 'TableColumnTypes' in k:
            head_dict[k] = ' '.join([x.strip() for x in v.strip().split(' ')])
        elif 'Origin' in k:
            head_dict[k] = v.lstrip()
        elif k in ('RangeStart', 'RangeEnd', 'AntennaBearing', 'ReferenceBearing', 'AngularResolution', 'SpatialResolution',
                   'FirstOrderMethod', 'BraggSmoothingPoints', 'BraggHasSecondOrder', 'MergedCount',
                   'RadialMinimumMergePoints', 'FirstOrderCalc', 'SpectraRangeCells', 'SpectraDopplerCells',
                   'TableColumns', 'TableRows',  'PatternResolution', 'CurrentVelocityLimit', 'TimeCoverage'):
            try:
                head_dict[k] = int(v)
            except ValueError:
                temp = v.split(' ')[0]
                try:
                    head_dict[k] = int(temp)
                except ValueError:
                    head_dict[k] = int(temp.split('.')[0])
        elif k in ('RangeResolutionKMeters', 'CTF', 'TransmitCenterFreqMHz', 'DopplerResolutionHzPerBin',
                   'RadialBraggPeakDropOff', 'RadialBraggPeakNull', 'RadialBraggNoiseThreshold', 'TransmitSweepRateHz',
                   'TransmitBandwidthKHz'):
            try:
                head_dict[k] = float(v)
            except ValueError:
                head_dict[k] = float(v.split(' ')[0])
        else:
            continue
    return head_dict


def create_ruv(radial_file_data, filename):
    # Create new radial file
    f = open(filename, 'w')

    # Write header
    for header_key, header_value in radial_file_data['header'].iteritems():
        f.write('%{}: {}\n'.format(header_key, header_value))

    for table in radial_file_data['tables'].keys():
        for table_key, table_value in radial_file_data['tables'][table].iteritems():
            if table_key is not 'data':
                f.write('%{}: {}\n'.format(table_key, table_value))

        radial_file_data['tables'][table]['data'].to_string(f, index=False, justify='center')

        if int(table) > 1:
            f.write('\n%TableEnd: {}\n'.format(table))
        else:
            f.write('\n%TableEnd: \n')
        f.write('%%\n')

    for footer_key, footer_value in radial_file_data['footer'].iteritems():
        if footer_key == 'ProcessingTool':
            for tool in radial_file_data['footer']['ProcessingTool']:
                f.write('%ProcessingTool: {}\n'.format(tool))
            f.write('%ProcessingTool: {}\n'.format('"hfr_processing/qc_radials.py"'))
        else:
            f.write('%{}: {}\n'.format(footer_key, footer_value))

    f.close()