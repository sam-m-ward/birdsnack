#Append model_files to path directory
import sys
sys.path.append('../model_files/')
from load_raw_data import *
from birdsnack_model import BIRDSNACK

#Load up light curves of fiducial sample of SNe, and metadata
dataloader = LOAD_DATA()
SNSsnpy_fiducial = dataloader.load_SNSsnpy('SNSsnpy_fiducial.pkl')
dfmeta           = dataloader.dfmeta
#Grab x65 low-to-moderate reddening SNe
bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_fiducial}, configname='loader_config.yaml', dfmeta=dfmeta, edit_dict={'additional_cut_parameters':{'BVcut':True,'BVcutval':1.0}})
bs.get_peak_mags()
bs.additional_cuts()
N65SNe = bs.sns
extra_drop_SNe = {sn:'BV less than 1.0' for sn in N65SNe}
#New fit to keep only x4 high reddening SNe, and fix intrinsic hyperparameters at the fiducial fit values
edit_dict        = {'additional_cut_parameters':{'extra_drop_SNe':extra_drop_SNe},
                    'analysis_parameters':{
                        'fixint':'AVExp_BVcut1.0',
                        'HBM_savekey':'HighRedSNe_AVGauss_FixedInt',
                        'disp_sigmaAV':1,
                        'AVprior':'Gauss',
                        'n_warmup':500,
                        'n_sampling':500
                        }
                    }
#Load into Bird-Snack
bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_fiducial}, configname='loader_config.yaml', dfmeta=dfmeta, edit_dict=edit_dict)
#Get peak magnitudes
bs.get_peak_mags()
#Perform additional sample cuts
bs.additional_cuts()
#Fit HBM to data
bs.fit_stan_model()#Rhat_threshold=1.05)
#Plot posterior samples
bs.plot_posterior_samples()
