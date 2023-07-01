#Append model_files to path directory
import sys
sys.path.append('../model_files/')
from load_raw_data import *
from CYCLE_DICTS import *
from birdsnack_model import BIRDSNACK, get_edit_dict
import argparse
parser = argparse.ArgumentParser(description="BirdSnack Analysis Cycle Script")
parser.add_argument("--mode",default='Science',help='Options are Science,CensoredData')
MODE   = parser.parse_args().mode
CYCLE_DICT = dict(zip(['Science','CensoredData'],[CYCLE_DICT_Science,CYCLE_DICT_CensoredData]))[MODE]

#Load up light curves of fiducial sample of SNe, and metadata
dataloader = LOAD_DATA()
SNSsnpy_fiducial = dataloader.load_SNSsnpy('SNSsnpy_fiducial.pkl')
dfmeta           = dataloader.dfmeta

for HBM_savekey in CYCLE_DICT['RUNS']:
    edit_dict = get_edit_dict(dataloader.choices,CYCLE_DICT,HBM_savekey)
    print ('###'*10)
    print (f"MODE{MODE}; Performing analysis on {edit_dict['analysis_parameters']['HBM_savekey']}")
    #Load into Bird-Snack
    bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_fiducial}, configname='loader_config.yaml', dfmeta=dfmeta, edit_dict=edit_dict)
    #Get peak magnitudes
    bs.get_peak_mags()
    #Perform additional sample cuts; set cutter=False (so SNe with magerrs>0.3mag are retained, ensures common sample fitted for consistency across pre-proc analysis variants)
    bs.additional_cuts(cutter=False)
    #Fit HBM to data
    if not os.path.exists(f"{bs.FITSpath}FIT{bs.choices['analysis_parameters']['HBM_savekey']}.pkl"):
        bs.fit_stan_model()
