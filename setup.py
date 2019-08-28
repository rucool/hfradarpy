from setuptools import setup

setup(
    name='codar_processing',
    version='0.3.1',
    packages=['scripts', 'scripts.sandbox', 'codar_processing.utilities', 'codar_processing', 'codar_processing.src', 'codar_processing.configs', 'codar_processing.methods', 'codar_processing.methods.waves', 'codar_processing.methods.totals', 'codar_processing.methods.radials', 'codar_processing.plotting'],
    url='https://github.com/rucool/codar_processing',
    license='MIT',
    author='mikesmith',
    author_email='michaesm@marine.rutgers.edu',
    description='Rutgers Center for Ocean Observing Leadership High Frequency Radar (CODAR) Processing toolbox'
)
