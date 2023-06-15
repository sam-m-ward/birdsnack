#Append model_files to path directory
import sys
sys.path.append('../model_files/')
from load_raw_data import *

#Load up fiducial sample
dataloader = LOAD_DATA()
SNSsnpy_fiducial = dataloader.load_SNSsnpy('SNSsnpy_fiducial.pkl')

from birdsnack_model import BIRDSNACK

bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_fiducial}, configname='loader_config.yaml', dfmeta=dataloader.dfmeta)

bs.get_peak_mags()
#bs.plot_lcs()
bs.additional_cuts()
bs.fit_stan_model()
