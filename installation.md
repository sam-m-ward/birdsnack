## Begin Bird-Snack installation by first setting up a snpy conda environment
#### Navigate to https://csp.obs.carnegiescience.edu/data/snpy/installing_snoopy2
#### Follow instructions, i.e. download snpy install script, and create new conda environment with choice of environment name, e.g. birdsnack, via:

`bash install-snpy birdsnack`

# Next, activate new conda environment and pip install additional packages
`conda activate birdsnack`

`pip install pandas==1.3.5 arviz==0.12.1 george==0.4.0 extinction==0.4.6 h5py==3.7.0 sncosmo==2.8.0 pystan==3.3.0 jupyter==1.0.0`

# Create a new folder to clone the birdsnack github repo into
`cd ../`

`mkdir birdsnack-clone`

`cd birdsnack-clone`

`git clone git@github.com:sam-m-ward/birdsnack.git`

# Next, update snoopy filters using complete filter set in Bird-Snack repo
`cd birdsnack/data/FILTERS/`

`python update_conda_env_filters.py`

# Bird-snack installation is now good to go, try running
`cd ../../analysis`

`python analyse_fiducial_sample.py`

#### Some scripts require BayeSN for simulating data.
#### To setup BirdSnack-BayeSN interaction see .txt instructions at sbc/BAYESN_GITHUB/README_BIRDSNACK_BAYESN_INSTRUCTIONS.txt
