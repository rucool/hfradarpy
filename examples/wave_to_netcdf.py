from hfradarpy.waves import Waves
import logging
import sys
from pathlib import Path

# Set up the parse_wave_files logger
logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)


def main(filename, save_dir, min=0.2, max=5):
    # Convert save_dir to a Path object if its just a string.
    save_dir = Path(save_dir)

    try:
        w = Waves(filename)
    except Exception:
            return

    # Flag heights less than min and greater than max
    w.flag_wave_heights(min=min, max=max, remove=True)

    sname = save_dir / w.file_name
    w.export(sname, file_type='netcdf', prepend_ext=True)


if __name__ == '__main__':
    data_root = (Path(__file__).parent.with_name('examples') / 'data').resolve()
    output_path = (Path(__file__).parent.with_name('examples') / 'output').resolve()

    wave_path = data_root / 'waves' / 'wls' / 'SEAB'
    wave_output = output_path / 'waves' / 'nc' / 'SEAB'

    for f in sorted(wave_path.glob('*.wls')):
        try:
            print(str(f))
            main(f, wave_output)
        except Exception:
            logger.exception('Exception in main(): ')
            exit(1)
