from setuptools import setup

setup(
    name='hfradar',
    version='0.3.3',
    packages=['scripts', 'scripts.sandbox', 'hfradar.utilities', 'hfradar', 'hfradar.src', 'hfradar.configs', 'hfradar.methods', 'hfradar.methods.waves', 'hfradar.methods.totals', 'hfradar.methods.radials', 'hfradar.plotting'],
    url='https://github.com/rucool/hfradarpy',
    license='MIT',
    author='mikesmith',
    author_email='michaesm@marine.rutgers.edu',
    description='Rutgers Center for Ocean Observing Leadership High Frequency Radar Processing toolbox'
)
