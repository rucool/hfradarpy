#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
import versioneer

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = ["xarray", "netcdf4", "pyproj", "click"]

test_requirements = [
    "pytest>=3",
]

setup(
    author="Michael Smith",
    author_email="michaesm@marine.rutgers.edu",
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    description="Toolbox to read in High Frequency Radar (HFR) files written in CODAR Tabular Format (CTF).",
    entry_points={
        "console_scripts": [
            "hfradarpy=hfradarpy.cli:main",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description_content_type="text/x-rst",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="hfradarpy",
    name="hfradarpy",
    packages=find_packages(include=["hfradarpy", "hfradarpy.*"]),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/rucool/hfradarpy",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
    extras_require={},
)
