# codar_processing
Rutgers Center for Ocean Observing Leadership High Frequency Radar (CODAR) Processing toolbox

# Install Miniconda
Download and follow installation instructions for the appropriate Miniconda installer from http://conda.pydata.org/miniconda.html. 

Make sure to add the channel, conda-forge, to your .condarc. You can find out more about conda-forge from their website: https://conda-forge.org/
 
You can do this with the following command:

`conda config --add channels conda-forge`

# Clone codar_processing repository
Either use git to clone:

`git clone https://github.com/rucool/codar_processing.git`

or download the zip from file from, https://github.com/rucool/codar_processing/archive/master.zip and extract to the directory of your choice


# Create codar_processing environment
Change your current working directory to the location that you downloaded codar_processing to. 

`cd /Users/mikesmith/Documents/git/codar_processing/`

Create conda environment from the included environment.yml file:

`conda create env -f environment.yml`

Once the environment is done building, you can activate the environment by typing:

    source activate codar_processing # OSX/Unix
    
# Install codar_processing toolbox to environment

Now we need to install the toolbox to the conda environment. We can do this as follows from the root directory of the codar_processing toolbox:

`pip install .`

The toolbox should now be installed to your conda environment.
