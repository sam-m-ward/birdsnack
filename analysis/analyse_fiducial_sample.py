#Append model_files to path directory
import sys
sys.path.append('../model_files/')
from load_raw_data import *
from birdsnack_model import BIRDSNACK

#Load up light curves of fiducial sample of SNe, and metadata
dataloader = LOAD_DATA()
SNSsnpy_fiducial = dataloader.load_SNSsnpy('SNSsnpy_fiducial.pkl')
dfmeta           = dataloader.dfmeta
edit_dict        = {}

#Load into Bird-Snack
bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_fiducial}, configname='loader_config.yaml', dfmeta=dfmeta, edit_dict=edit_dict)

#Get peak magnitudes
bs.get_peak_mags()

#Perform additional sample cuts
bs.additional_cuts()

#Plot Mag Deviations
bs.plot_mag_deviations()

#Plot Colours
bs.plot_colour_corner()

#Fit HBM to data
bs.fit_stan_model()

#Plot posterior samples
bs.plot_posterior_samples()
