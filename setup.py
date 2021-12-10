#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
import versioneer

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Click>=7.0', 
    'xarray==0.16.1',
    'netcdf4==1.5.4',
    'numpy==1.20',
    'pyparsing==2.4.7',
    'scipy',
    'geopandas==0.10.2',
    'oceans',
    'jupyter',
    'numba',
    'pytest',
    'cmocean',
    'pymongo',
    'dask',
    'python-dateutil',
    'rtree'
    ]

test_requirements = ['pytest>=3', ]

setup(
    author="Michael Smith",
    author_email='michaesm@marine.rutgers.edu',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Toolbox to read in High Frequency Radar (HFR) files written in CODAR Tabular Format (CTF).",
    entry_points={
        'console_scripts': [
            'hfradarpy=hfradarpy.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='hfradarpy',
    name='hfradarpy',
    packages=find_packages(include=['hfradarpy', 'hfradarpy.*'], exclude="tests"),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/rucool/hfradarpy',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
    extras_require={
        "plot": ['cartopy']
    }
)
