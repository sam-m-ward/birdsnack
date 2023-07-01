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
#edit_dict = {  'preproc_parameters'    :   {'DF_savekey':'uBVriJH','pblist':[s for s in 'uBVriJH']},
#                'analysis_parameters'   :   {'HBM_savekey':'uBVriJH_CensoredCut1.0','lam_choice':'central','CensoredData':True,'CensoredCut':1.0}}
#edit_dict = {'analysis_parameters':{'HBM_savekey':'Fiducial_AVGamma_CensoredCut1.0','CensoredData':True,'CensoredCut':1.0,'AVprior':'Gamma','n_sampling':10000}}
#edit_dict = {'analysis_parameters'      :   {'HBM_savekey':'LowRVs_CensoredCut1.0','CensoredData':True,'CensoredCut':1.0},
#            'additional_cut_parameters' :   {'extra_drop_SNe':{sn:'Low RVs' for sn in ['2009ds','16abc']}}}

#Load into Bird-Snack
bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_fiducial}, configname='loader_config.yaml', dfmeta=dfmeta, edit_dict=edit_dict)

#Get peak magnitudes
bs.get_peak_mags()

#Perform additional sample cuts
bs.additional_cuts()

#Plot Mag Deviations
#bs.plot_mag_deviations()

#Plot Colours
#bs.plot_colour_corner()

#bs.plot_lcs()

#Fit HBM to data
bs.fit_stan_model()

#Plot posterior samples
bs.plot_posterior_samples()
