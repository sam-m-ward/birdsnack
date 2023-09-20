#Append model_files to path directory
import sys
sys.path.append('../model_files/')
from load_raw_data import *
from birdsnack_model import BIRDSNACK

#Load up light curves of fiducial sample of SNe, and metadata
dataloader = LOAD_DATA()
SNSsnpy_fiducial = dataloader.load_SNSsnpy('SNSsnpy_fiducial.pkl')
dfmeta           = dataloader.dfmeta
edit_dict        = {'additional_cut_parameters':{'BVcut':True,'BVcutval':1.0},
                    'analysis_parameters':{
                        'fixint':'AVExp_BVcut1.0',
                        'HBM_savekey':'AVExp_BVcut1.0_FixedInt'
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
