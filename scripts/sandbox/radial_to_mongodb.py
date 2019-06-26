from mongoengine import *
import datetime as dt
from codar_processing.src.radials import Radial

connect('codar')


class Site(Document):
    name = StringField(required=True, max_length=4, unique=True)
    center_frequency = DecimalField(required=True)
    date_added = DateTimeField(default=dt.datetime.now)
    meta = {
        'ordering': ['-name']
    }


class RadialMetadata(DynamicEmbeddedDocument):
    date_added = DateTimeField(default=dt.datetime.now)


class RadialDiagnostics(DynamicEmbeddedDocument):
    date_added = DateTimeField(default=dt.datetime.now)


class HardwareDiagnostics(DynamicEmbeddedDocument):
    date_added = DateTimeField(default=dt.datetime.now)


class RadialFile(Document):
    filename = StringField(max_length=120, required=True, unique=True)
    site_code = ReferenceField(Site, required=True, primary_key=True)
    radial_metadata = EmbeddedDocumentField(RadialMetadata)
    hardware_diagnostics = EmbeddedDocumentField(HardwareDiagnostics)
    radial_diagnostics = EmbeddedDocumentField(RadialDiagnostics)
    date_added = DateTimeField(default=dt.datetime.now)
    # meta = {'allow_inheritance': True}


def main(file):
    r = Radial(file)
    r.clean_header()
    print(r.file_name)
    r.metadata['filename'] = r.file_name

    # Upload site information to database
    try:
        site = Site(name=r.metadata['Site'], center_frequency=r.metadata['TransmitCenterFreqMHz'])
        site.save()
    except NotUniqueError:  # except if its already uploaded
        site = Site.objects(name=r.metadata['Site'])[0]

    r.metadata['site_code'] = site.id

    hardware_diagnostics = HardwareDiagnostics(**r.diagnostics_hardware.to_dict(orient='list'))
    radial_diagnostics = RadialDiagnostics(**r.diagnostics_radial.to_dict(orient='list'))
    radial_metadata = RadialMetadata(**r.metadata)

    object_info = {}
    object_info['filename'] = r.file_name
    object_info['site_code'] = site.id
    object_info['radial_metadata'] = radial_metadata
    object_info['hardware_diagnostics'] = hardware_diagnostics
    object_info['radial_diagnostics'] = radial_diagnostics

    # Upload the radial file header information
    RadialFile(**object_info).save()


if __name__ == '__main__':
    from glob import glob
    import os

    radial_dir = '/Volumes/home/codaradm/data/radials/BRMR/2018_03/'
    file_list = sorted(glob(os.path.join(radial_dir, '*.ruv')))
    for radial in file_list:
        main(radial)



